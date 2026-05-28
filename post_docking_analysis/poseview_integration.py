"""
PoseView integration module for post-docking analysis pipeline.

Generates 2D interaction diagrams using the proteins.plus PoseView API.
Supports both PDB database codes and local PDB file uploads.

Reference: https://proteins.plus/help/poseview_rest
"""
import requests
import time
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
import re
from collections import Counter

logger = logging.getLogger(__name__)


@dataclass
class PoseViewResult:
    """Results from PoseView API analysis."""
    complex_name: str
    pdb_file: Optional[Path] = None
    pdb_code: Optional[str] = None
    ligand: Optional[str] = None
    endpoint: Optional[str] = None
    success: bool = False
    job_id: Optional[str] = None
    png_file: Optional[Path] = None
    svg_file: Optional[Path] = None
    pdf_file: Optional[Path] = None
    error_message: Optional[str] = None


class PoseViewClient:
    """
    REST API client for proteins.plus PoseView service.
    
    This client handles:
    - Uploading local PDB files
    - Creating PoseView jobs
    - Polling for job completion
    - Downloading result images (PNG, SVG, PDF)
    """
    
    BASE_URL = "https://proteins.plus/api"
    
    def __init__(self, config: Dict = None):
        """
        Initialize PoseView client.
        
        Parameters
        ----------
        config : Dict, optional
            Configuration dictionary with settings
        """
        self.config = config or {}
        self.poll_interval = self.config.get('poll_interval_seconds', 5)
        self.timeout = self.config.get('timeout_seconds', 300)
        self.max_retries = self.config.get('max_retries', 3)
        self.upload_poll_interval = int(self.config.get('upload_poll_interval_seconds', 2))
        self.upload_timeout = int(self.config.get('upload_timeout_seconds', 180))
        self.last_error_message: str = ""
        self._job_status_urls: Dict[str, str] = {}
        self._job_endpoint_labels: Dict[str, str] = {}
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        })
        self.job_endpoints = self._build_job_endpoints()

    def _build_job_endpoints(self) -> List[Dict[str, str]]:
        """
        Build ordered endpoint candidates for PoseView job creation/status.
        """
        configured = self.config.get("job_endpoints")
        if isinstance(configured, list) and configured:
            endpoints: List[Dict[str, str]] = []
            for item in configured:
                if not isinstance(item, dict):
                    continue
                base = str(item.get("base_url") or self.BASE_URL).rstrip("/")
                path = str(item.get("path") or "").strip("/")
                payload_key = str(item.get("payload_key") or "").strip()
                label = str(item.get("label") or f"{base}/{path}").strip()
                if not path or not payload_key:
                    continue
                endpoints.append(
                    {
                        "base_url": base,
                        "path": path,
                        "payload_key": payload_key,
                        "label": label,
                    }
                )
            if endpoints:
                return endpoints

        # Default: try classic endpoint first, then PoseEdit variants.
        return [
            {
                "base_url": "https://proteins.plus/api",
                "path": "poseview_rest",
                "payload_key": "poseview",
                "label": "proteins.plus/poseview_rest",
            },
            {
                "base_url": "https://proteins.plus/api",
                "path": "poseview2_rest",
                "payload_key": "poseview2",
                "label": "proteins.plus/poseview2_rest",
            },
            {
                "base_url": "https://poseedit.proteins.plus/api",
                "path": "poseview2_rest",
                "payload_key": "poseview2",
                "label": "poseedit.proteins.plus/poseview2_rest",
            },
        ]

    @staticmethod
    def _safe_json(response: requests.Response) -> Dict[str, Any]:
        try:
            payload = response.json()
            return payload if isinstance(payload, dict) else {}
        except Exception:
            return {}

    @staticmethod
    def _extract_location(data: Dict[str, Any]) -> str:
        location = str(data.get("location") or "").strip()
        return location

    @staticmethod
    def _extract_file_id(data: Dict[str, Any], fallback_location: str = "") -> Optional[str]:
        """
        Extract uploaded file ID from upload/status payload.
        """
        if not isinstance(data, dict):
            return None
        file_id = data.get('id') or data.get('hashid')
        if file_id:
            return str(file_id)
        location = str(data.get('location') or fallback_location or "")
        if location:
            return location.rstrip("/").split("/")[-1]
        return None

    @staticmethod
    def _safe_int(value: Any) -> Optional[int]:
        """Parse status-like values without raising on textual API payloads."""
        if isinstance(value, bool):
            return None
        if isinstance(value, int):
            return int(value)
        if isinstance(value, float):
            try:
                return int(value)
            except Exception:
                return None
        text = str(value or "").strip()
        if not text:
            return None
        if text.isdigit() or (text.startswith("-") and text[1:].isdigit()):
            try:
                return int(text)
            except Exception:
                return None
        return None

    @classmethod
    def _normalize_status(
        cls,
        data: Dict[str, Any],
        http_status: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Normalize mixed numeric/text API status into canonical states:
        `completed`, `pending`, `rate_limited`, `failed`.
        """
        payload = data if isinstance(data, dict) else {}
        code = cls._safe_int(payload.get("status_code"))
        if code is None:
            code = cls._safe_int(payload.get("code"))
        if code is None:
            code = cls._safe_int(http_status)

        token_candidates = [
            payload.get("status"),
            payload.get("state"),
            payload.get("job_status"),
            payload.get("phase"),
            payload.get("message"),
            payload.get("detail"),
            payload.get("status_code"),
        ]
        status_text = " ".join(str(value or "").strip().lower() for value in token_candidates).strip()

        success_tokens = {
            "done", "complete", "completed", "success", "succeeded", "finished", "ready", "ok",
        }
        pending_tokens = {
            "accepted", "queued", "queue", "processing", "running", "pending", "submitted", "in progress",
            "in_progress", "working", "started",
        }
        rate_tokens = {"rate limit", "too many requests", "throttle"}

        if any(token in status_text for token in rate_tokens):
            return {"state": "rate_limited", "status_code": code or 429, "status_text": status_text}
        if any(token in status_text for token in success_tokens):
            return {"state": "completed", "status_code": code or 200, "status_text": status_text}
        if any(token in status_text for token in pending_tokens):
            return {"state": "pending", "status_code": code or 202, "status_text": status_text}

        if code == 429:
            return {"state": "rate_limited", "status_code": 429, "status_text": status_text}
        if code in (200, 201):
            return {"state": "completed", "status_code": code, "status_text": status_text}
        if code == 202:
            return {"state": "pending", "status_code": 202, "status_text": status_text}

        if code is not None and code >= 400:
            return {"state": "failed", "status_code": code, "status_text": status_text}
        return {"state": "failed", "status_code": code or 500, "status_text": status_text}

    def _wait_for_upload_ready(self, location_url: str) -> Dict[str, Any]:
        """
        Poll upload location until structure preprocessing is complete.
        """
        start = time.time()
        last_payload: Dict[str, Any] = {}
        while True:
            if time.time() - start > self.upload_timeout:
                return last_payload
            try:
                response = self.session.get(location_url, timeout=30)
                payload = self._safe_json(response)
            except requests.exceptions.RequestException:
                return last_payload
            except Exception:
                return last_payload

            last_payload = payload if isinstance(payload, dict) else {}
            status = self._normalize_status(last_payload, response.status_code)
            state = status.get("state")
            if state == "completed":
                return last_payload
            if state == "rate_limited":
                time.sleep(max(10, self.upload_poll_interval))
                continue
            if state != "pending":
                return last_payload
            time.sleep(self.upload_poll_interval)
    
    def upload_pdb_file(self, pdb_file: Path) -> Optional[str]:
        """
        Upload a local PDB file to proteins.plus.
        
        Parameters
        ----------
        pdb_file : Path
            Path to the local PDB file
            
        Returns
        -------
        str or None
            The uploaded file ID (to use as pdbCode), or None if failed
        """
        pdb_file = Path(pdb_file)
        
        if not pdb_file.exists():
            logger.error(f"PDB file not found: {pdb_file}")
            return None
        
        logger.info(f"📤 Uploading PDB file: {pdb_file.name}")
        
        try:
            # Use multipart form data for file upload
            with open(pdb_file, 'rb') as f:
                files = {
                    'pdb_file[pathvar]': (pdb_file.name, f, 'chemical/x-pdb')
                }
                headers = {'Accept': 'application/json'}
                
                response = requests.post(
                    f"{self.BASE_URL}/pdb_files_rest",
                    files=files,
                    headers=headers,
                    timeout=60
                )
            
            if response.status_code in [200, 201, 202]:
                data = self._safe_json(response)
                location = str(data.get("location") or "")
                status = self._normalize_status(data, response.status_code)

                # Poll asynchronous uploads until server-side parsing is completed.
                if status.get("state") == "pending" and location:
                    data = self._wait_for_upload_ready(location) or data

                # The response should contain the file ID or location.
                file_id = self._extract_file_id(data, fallback_location=location)
                if file_id:
                    logger.info(f"✅ Uploaded successfully: {file_id}")
                    return file_id
                else:
                    logger.warning(f"⚠️ Upload response missing ID: {data}")
                    return None
            else:
                logger.error(f"❌ Upload failed: {response.status_code} - {response.text}")
                return None
                
        except requests.exceptions.RequestException as e:
            logger.error(f"❌ Upload request failed: {e}")
            return None
    
    def create_poseview_job(
        self,
        pdb_code: str,
        ligand: str = ""
    ) -> Optional[str]:
        """
        Create a new PoseView job.
        
        Parameters
        ----------
        pdb_code : str
            PDB code or uploaded file ID
        ligand : str
            Ligand identifier (e.g., "JE2_A_701")
            
        Returns
        -------
        str or None
            Job ID/location, or None if failed
        """
        logger.info(f"🔬 Creating PoseView job: {pdb_code}, ligand: {ligand}")
        self.last_error_message = ""
        
        ligand_text = str(ligand or "").strip()
        errors: List[str] = []

        for endpoint in self.job_endpoints:
            poseview_payload: Dict[str, Any] = {"pdbCode": pdb_code}
            # For auto-detection, omit ligand entirely instead of sending empty text.
            if ligand_text:
                poseview_payload["ligand"] = ligand_text
            payload = {endpoint["payload_key"]: poseview_payload}
            endpoint_url = f"{endpoint['base_url'].rstrip('/')}/{endpoint['path'].strip('/')}"

            try:
                response = self.session.post(
                    endpoint_url,
                    json=payload,
                    timeout=30
                )
            except requests.exceptions.RequestException as e:
                errors.append(f"{endpoint['label']}: request error ({e})")
                continue

            data = self._safe_json(response)
            status = self._normalize_status(data, response.status_code)
            status_code = int(status.get("status_code", response.status_code))
            location = self._extract_location(data)

            if status.get("state") in {"completed", "pending"}:
                job_id = location.split('/')[-1] if location else self._extract_file_id(data)
                if not job_id:
                    # Some responses provide message but no explicit location.
                    errors.append(f"{endpoint['label']}: missing location in success payload")
                    continue
                status_url = location if location.startswith("http") else f"{endpoint_url.rstrip('/')}/{job_id}"
                self._job_status_urls[job_id] = status_url
                self._job_endpoint_labels[job_id] = endpoint["label"]
                logger.info(f"✅ Job created: {job_id} via {endpoint['label']}")
                return job_id

            message = str(data.get('message', response.text or 'Unknown error'))
            errors.append(f"{endpoint['label']}: {status_code} {message}")

        self.last_error_message = " | ".join(errors) if errors else "Unknown error"
        logger.error(f"❌ Job creation failed: {self.last_error_message}")
        return None
    
    def get_job_status(self, job_id: str) -> Dict:
        """
        Get the status of a PoseView job.
        
        Parameters
        ----------
        job_id : str
            The job ID
            
        Returns
        -------
        Dict
            Job status and results
        """
        try:
            status_url = self._job_status_urls.get(job_id, f"{self.BASE_URL}/poseview_rest/{job_id}")
            response = self.session.get(
                status_url,
                timeout=30
            )
            data = self._safe_json(response)
            normalized = self._normalize_status(data, response.status_code)
            data["status_code"] = normalized["status_code"]
            data["status_state"] = normalized["state"]
            data["status_text"] = normalized["status_text"]
            data["_http_status"] = response.status_code
            return data
        except requests.exceptions.RequestException as e:
            logger.error(f"❌ Status request failed: {e}")
            return {'status_code': 500, 'error': str(e)}
    
    def wait_for_completion(self, job_id: str) -> Dict:
        """
        Wait for a job to complete, polling periodically.
        
        Parameters
        ----------
        job_id : str
            The job ID
            
        Returns
        -------
        Dict
            Final job result
        """
        start_time = time.time()
        
        while True:
            elapsed = time.time() - start_time
            if elapsed > self.timeout:
                return {'status_code': 408, 'error': 'Timeout waiting for job'}
            
            result = self.get_job_status(job_id)
            normalized = self._normalize_status(
                result,
                self._safe_int(result.get("_http_status")) or self._safe_int(result.get("status_code")) or 500
            )
            status_code = int(normalized.get("status_code", 500))
            status_state = str(normalized.get("state") or "")
            result["status_code"] = status_code
            result["status_state"] = status_state
            result["status_text"] = normalized.get("status_text", "")
            
            if status_state == "completed":
                # Job completed
                return result
            elif status_state == "pending":
                # Still processing
                logger.debug(f"Job {job_id} still processing...")
                time.sleep(self.poll_interval)
            elif status_state == "rate_limited":
                # Rate limited
                logger.warning("Rate limited, waiting 30 seconds...")
                time.sleep(30)
            else:
                # Error
                return result
    
    def download_result_image(
        self,
        url: str,
        output_path: Path,
        image_type: str = "png"
    ) -> bool:
        """
        Download a result image from the API.
        
        Parameters
        ----------
        url : str
            URL of the image
        output_path : Path
            Local path to save the image
        image_type : str
            Type of image (png, svg, pdf)
            
        Returns
        -------
        bool
            True if download successful
        """
        try:
            response = requests.get(url, timeout=60)
            
            if response.status_code == 200:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                logger.info(f"✅ Downloaded {image_type}: {output_path.name}")
                return True
            else:
                logger.warning(f"⚠️ Download failed: {response.status_code}")
                return False
                
        except requests.exceptions.RequestException as e:
            logger.error(f"❌ Download failed: {e}")
            return False


class PoseViewAnalyzer:
    """
    PoseView analyzer for generating 2D interaction diagrams.
    
    This class wraps the PoseView API client and provides methods
    for analyzing protein-ligand complexes from local PDB files.
    """
    
    def __init__(self, config: Dict = None):
        """
        Initialize PoseView analyzer.
        
        Parameters
        ----------
        config : Dict, optional
            Configuration dictionary
        """
        self.config = config or {}
        self.client = PoseViewClient(config)
        self.output_formats = self.config.get('output_formats', ['png', 'svg', 'pdf'])

    @staticmethod
    def _normalize_resname(raw_name: str) -> str:
        token = "".join(ch for ch in str(raw_name).upper() if ch.isalnum())
        if not token:
            return "UNK"
        candidates: List[str] = []
        if len(token) >= 3:
            candidates.extend(
                [
                    token[:3],
                    f"{token[0]}{token[-2:]}",
                    f"{token[:2]}{token[-1]}",
                    token[-3:],
                ]
            )
        elif len(token) == 2:
            candidates.extend([token + "X", f"{token[0]}X{token[1]}", f"X{token}"])
        else:
            candidates.extend([token + "XX", f"X{token}X", f"XX{token}"])

        normalized: List[str] = []
        seen = set()
        for candidate in candidates:
            code = "".join(ch for ch in candidate.upper() if ch.isalnum())[:3]
            if len(code) < 3:
                code = code.ljust(3, "X")
            if code and code not in seen:
                seen.add(code)
                normalized.append(code)

        aa_codes = {
            "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX",
            "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
        }
        for code in normalized:
            if code not in aa_codes:
                return code
        return normalized[0] if normalized else "UNK"

    def _infer_ligand_from_complex_name(self, complex_name: str) -> Optional[tuple]:
        """
        Infer ligand ID from naming conventions when PDB residue name is generic.
        """
        text = str(complex_name or "")
        match = re.search(r"_ligand_([A-Za-z0-9]{1,6})_([A-Za-z])_(\d+)", text)
        if match:
            return (
                self._normalize_resname(match.group(1)),
                match.group(2).upper(),
                match.group(3),
            )

        stem = Path(text).stem
        if not stem:
            return None
        token = stem.split("_")[-1]
        resname = self._normalize_resname(token)
        if resname == "UNK":
            return None
        return (resname, "B", "1")

    @staticmethod
    def _is_invalid_ligand(resname: str) -> bool:
        invalid = {"", "UNK", "UNX", "LIG"}
        aa_codes = {
            "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX",
            "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
        }
        token = str(resname or "").strip().upper()
        return token in invalid or token in aa_codes
    

    
    def _clean_pdb_file(self, input_file: Path) -> Path:
        """
        Clean PDB file by removing non-standard lines (ROOT, BRANCH, etc).
        Returns path to cleaned file (temp file if changes needed).
        """
        try:
            with open(input_file, 'r') as f:
                lines = f.readlines()
            
            cleaned_lines = []
            needs_cleaning = False
            
            for line in lines:
                if line.startswith(('ROOT', 'BRANCH', 'TORSDOF', 'ENDROOT', 'ENDBRANCH')):
                    needs_cleaning = True
                    continue
                cleaned_lines.append(line)
            
            if needs_cleaning:
                import tempfile
                cleaned_file = Path(tempfile.mktemp(suffix='.pdb', prefix=f"{input_file.stem}_cleaned_"))
                with open(cleaned_file, 'w') as f:
                    f.writelines(cleaned_lines)
                return cleaned_file
            
            return input_file
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to clean PDB file: {e}")
            return input_file

    def _detect_ligand(self, pdb_file: Path) -> tuple:
        """
        Detect the ligand in the PDB file.
        Returns (resname, chain, resnum).
        """
        try:
            candidates = self._detect_ligand_candidates(pdb_file)
            if not candidates:
                return ('UNK', 'B', '1')
            
            best_ligand = candidates[0]
            logger.info(f"🔍 Detected ligand: {best_ligand}")
            return best_ligand
            
        except Exception as e:
            logger.warning(f"⚠️ Failed to detect ligand: {e}")
            return ('UNK', 'B', '1')

    def _detect_ligand_candidates(self, pdb_file: Path) -> List[tuple]:
        """
        Detect plausible ligand candidates ordered by atom-count support.
        """
        ligand_counts: Dict[tuple, int] = {}
        exclude = {'HOH', 'WAT', 'NA', 'CL', 'MG', 'ZN', 'CA', 'MN', 'K'}

        with open(pdb_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.startswith('HETATM'):
                    continue
                raw_resname = line[17:20].strip()
                resname = self._normalize_resname(raw_resname)
                chain = (line[21:22].strip() or 'A').upper()
                resnum = (line[22:26].strip() or '1')

                if raw_resname.upper() in exclude:
                    continue
                if self._is_invalid_ligand(resname):
                    continue

                key = (resname, chain, resnum)
                ligand_counts[key] = ligand_counts.get(key, 0) + 1

        ranked = sorted(ligand_counts.items(), key=lambda kv: kv[1], reverse=True)
        return [item[0] for item in ranked]

    def _build_ligand_id_candidates(
        self,
        complex_name: str,
        pdb_file: Path,
        ligand_name: Optional[str],
        chain_id: Optional[str],
        residue_num: Optional[int],
    ) -> List[str]:
        """
        Build ordered ligand identifier attempts for PoseView job creation.
        """
        triples: List[tuple] = []
        seen_triples = set()

        def add_triple(name: Optional[str], chain: Optional[str], num: Optional[str]) -> None:
            norm_name = self._normalize_resname(name or "")
            if self._is_invalid_ligand(norm_name):
                return
            norm_chain = str(chain or "A").strip().upper() or "A"
            norm_num = str(num or "1").strip() or "1"
            key = (norm_name, norm_chain, norm_num)
            if key in seen_triples:
                return
            seen_triples.add(key)
            triples.append(key)

        add_triple(ligand_name, chain_id, residue_num)

        detected = self._detect_ligand(pdb_file)
        add_triple(*detected)

        inferred = self._infer_ligand_from_complex_name(complex_name)
        if inferred:
            add_triple(*inferred)

        for candidate in self._detect_ligand_candidates(pdb_file):
            add_triple(*candidate)

        max_candidates = int(self.config.get('max_ligand_candidates', 8))
        seen_ids = set()
        ids: List[str] = []

        for name, chain, num in triples[:max_candidates]:
            candidate_forms = [
                f"{name}_{chain}_{num}",  # canonical in docs
                f"{name}_{chain}",        # fallback: no residue number
                f"{name}",                # fallback: residue name only
            ]
            for candidate in candidate_forms:
                candidate = candidate.strip("_")
                if not candidate or candidate in seen_ids:
                    continue
                seen_ids.add(candidate)
                ids.append(candidate)

        # Final fallback: let server auto-select ligand.
        ids.append("")
        return ids

    @staticmethod
    def _get_result_url(payload: Dict[str, Any], keys: List[str]) -> str:
        """
        Resolve a result URL from multiple possible API field names.
        """
        for key in keys:
            value = payload.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()
        return ""

    def analyze_complex(
        self,
        pdb_file: Path,
        output_dir: Path,
        ligand_name: str = None,
        chain_id: str = None,
        residue_num: int = None
    ) -> PoseViewResult:
        """
        Analyze a protein-ligand complex and generate 2D diagram.
        """
        pdb_file = Path(pdb_file)
        output_dir = Path(output_dir)
        complex_name = pdb_file.stem
        
        # Create output directory
        complex_output_dir = output_dir / complex_name
        complex_output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"🔬 Analyzing: {complex_name}")
        
        result = PoseViewResult(
            complex_name=complex_name,
            pdb_file=pdb_file
        )
        
        try:
            # Step 0: Clean PDB file
            cleaned_pdb = self._clean_pdb_file(pdb_file)

            # Step 1: Detect ligand if not provided
            if not all([ligand_name, chain_id, residue_num]):
                d_name, d_chain, d_num = self._detect_ligand(cleaned_pdb)
                inferred = self._infer_ligand_from_complex_name(complex_name)
                if inferred is not None and self._is_invalid_ligand(d_name):
                    d_name, d_chain, d_num = inferred
                ligand_name = ligand_name or d_name
                chain_id = chain_id or d_chain
                residue_num = residue_num or d_num

            ligand_candidates = self._build_ligand_id_candidates(
                complex_name=complex_name,
                pdb_file=cleaned_pdb,
                ligand_name=ligand_name,
                chain_id=chain_id,
                residue_num=residue_num,
            )
            if not ligand_candidates:
                result.error_message = (
                    f"Could not resolve a valid ligand identifier for {complex_name}"
                )
                logger.warning(f"⚠️  {result.error_message}")
                return result

            # Step 2: Upload the PDB file
            file_id = self.client.upload_pdb_file(cleaned_pdb)

            # Clean up temp file if created
            if cleaned_pdb != pdb_file and cleaned_pdb.exists():
                try:
                    cleaned_pdb.unlink()
                except:
                    pass

            if not file_id:
                result.error_message = "Failed to upload PDB file"
                return result

            result.pdb_code = file_id

            # Step 3: Create PoseView job (try multiple ligand identifiers).
            job_id = None
            selected_ligand = None
            for idx, ligand_id in enumerate(ligand_candidates, start=1):
                shown = ligand_id or "<auto>"
                if idx > 1:
                    logger.info(f"↪️ PoseView retry {idx}/{len(ligand_candidates)} with ligand {shown}")
                job_id = self.client.create_poseview_job(file_id, ligand_id)
                if job_id:
                    selected_ligand = ligand_id or "AUTO"
                    break

            if not job_id:
                attempted = ", ".join(c if c else "<auto>" for c in ligand_candidates)
                detail = self.client.last_error_message or "unknown error"
                result.error_message = (
                    f"Failed to create PoseView job. "
                    f"Attempted ligands: {attempted}. Last error: {detail}"
                )
                logger.warning(f"⚠️  {result.error_message}")
                return result

            result.ligand = selected_ligand
            result.job_id = job_id
            result.endpoint = self.client._job_endpoint_labels.get(job_id)
            
            # Step 4: Wait for completion
            job_result = self.client.wait_for_completion(job_id)

            terminal_state = str(job_result.get("status_state") or "")
            if terminal_state != "completed":
                failure_state = terminal_state or "failed"
                failure_msg = str(job_result.get('message') or job_result.get("error") or 'Job failed')
                result.error_message = (
                    f"PoseView job failed "
                    f"(state={failure_state}, status={job_result.get('status_code')}) "
                    f"{failure_msg}"
                )
                return result
            
            # Step 5: Download result images
            png_url = self._get_result_url(
                job_result,
                ["result_png_picture", "result_png", "result_png_url"],
            )
            svg_url = self._get_result_url(
                job_result,
                ["result_svg_picture", "result_svg", "result_svg_url"],
            )
            pdf_url = self._get_result_url(
                job_result,
                ["result_pdf_picture", "result_pdf", "result_pdf_url"],
            )

            if 'png' in self.output_formats and png_url:
                png_path = complex_output_dir / f"{complex_name}_2d_diagram.png"
                if self.client.download_result_image(
                    png_url, png_path, 'png'
                ):
                    result.png_file = png_path
            
            if 'svg' in self.output_formats and svg_url:
                svg_path = complex_output_dir / f"{complex_name}_2d_diagram.svg"
                if svg_url.startswith("http"):
                    if self.client.download_result_image(svg_url, svg_path, 'svg'):
                        result.svg_file = svg_path
                else:
                    # Some APIs may return inline SVG content.
                    try:
                        svg_path.write_text(svg_url, encoding="utf-8")
                        result.svg_file = svg_path
                    except Exception:
                        pass
            
            if 'pdf' in self.output_formats and pdf_url:
                pdf_path = complex_output_dir / f"{complex_name}_2d_diagram.pdf"
                if self.client.download_result_image(
                    pdf_url, pdf_path, 'pdf'
                ):
                    result.pdf_file = pdf_path
            
            result.success = any([result.png_file, result.svg_file, result.pdf_file])
            
            if result.success:
                logger.info(f"  ✅ Generated 2D diagram for {complex_name}")
            else:
                result.error_message = "No images downloaded"
                logger.warning(f"  ⚠️ No images generated for {complex_name}")
            
            return result
            
        except Exception as e:
            logger.error(f"❌ Error analyzing {complex_name}: {e}")
            result.error_message = str(e)
            return result
    
    def analyze_directory(
        self,
        poses_dir: Path,
        output_dir: Path,
        binding_category: Optional[str] = None,
        ligand_name: str = None,
        chain_id: str = None
    ) -> Dict[str, PoseViewResult]:
        """
        Analyze all PDB files in a directory.
        """
        poses_dir = Path(poses_dir)
        output_dir = Path(output_dir)
        
        if binding_category:
            output_dir = output_dir / binding_category
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find PDB files
        pdb_files = sorted(poses_dir.glob("*.pdb"))
        
        if not pdb_files:
            logger.warning(f"⚠️ No PDB files found in {poses_dir}")
            return {}
        
        max_complexes = int(self.config.get("max_complexes", 0) or 0)
        if max_complexes > 0:
            pdb_files = pdb_files[:max_complexes]

        logger.info(f"📁 Processing {len(pdb_files)} PDB files from {poses_dir.name}")
        
        results = {}
        
        for pdb_file in pdb_files:
            result = self.analyze_complex(
                pdb_file, output_dir, ligand_name, chain_id
            )
            results[result.complex_name] = result
            
            # Rate limiting - wait between requests
            time.sleep(2)
        
        successful = sum(1 for r in results.values() if r.success)
        logger.info(f"✅ Completed: {successful}/{len(results)} successful")
        
        return results
    
    def analyze_best_poses_directory(
        self,
        best_poses_dir: Path,
        output_dir: Path,
        ligand_name: str = None,
        chain_id: str = None
    ) -> Dict[str, Dict[str, PoseViewResult]]:
        """
        Analyze the complete best_poses directory structure.
        """
        best_poses_dir = Path(best_poses_dir)
        output_dir = Path(output_dir)
        
        categories = ['strong_binders', 'moderate_binders', 'weak_binders']
        all_results = {}
        
        for category in categories:
            category_dir = best_poses_dir / category
            
            if category_dir.exists():
                logger.info(f"\n📂 Processing {category}...")
                results = self.analyze_directory(
                    category_dir, output_dir, category, ligand_name, chain_id
                )
                all_results[category] = results
            else:
                logger.info(f"ℹ️ Skipping {category} (not found)")
        
        # Generate summary
        self._generate_summary(all_results, output_dir)
        
        return all_results
    
    def _generate_summary(
        self,
        all_results: Dict[str, Dict[str, PoseViewResult]],
        output_dir: Path
    ):
        """Generate JSON summary of results."""
        summary = {
            'total_analyzed': 0,
            'total_successful': 0,
            'by_category': {},
            'complexes': []
        }
        
        for category, results in all_results.items():
            category_info = {
                'count': len(results),
                'successful': sum(1 for r in results.values() if r.success)
            }
            summary['by_category'][category] = category_info
            summary['total_analyzed'] += category_info['count']
            summary['total_successful'] += category_info['successful']
            
            for name, result in results.items():
                summary['complexes'].append({
                    'name': name,
                    'category': category,
                    'success': result.success,
                    'endpoint': result.endpoint,
                    'png_file': str(result.png_file) if result.png_file else None,
                    'error': result.error_message
                })
        
        summary_file = output_dir / 'poseview_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"📊 Summary saved to {summary_file}")


def run_poseview_analysis(
    poses_dir: Path,
    output_dir: Path,
    config: Dict = None
) -> Dict:
    """
    Main entry point for running PoseView analysis.
    
    Parameters
    ----------
    poses_dir : Path
        Path to poses directory
    output_dir : Path
        Output directory for results
    config : Dict, optional
        Configuration dictionary
        
    Returns
    -------
    Dict
        Analysis summary
    """
    print("🔬 Running PoseView 2D Diagram Analysis...")
    
    poses_dir = Path(poses_dir)
    output_dir = Path(output_dir)
    
    analyzer = PoseViewAnalyzer(config)
    
    # Check for categorized structure
    subdirs = ['strong_binders', 'moderate_binders', 'weak_binders']
    is_categorized = any((poses_dir / d).exists() for d in subdirs)
    
    if is_categorized:
        results = analyzer.analyze_best_poses_directory(poses_dir, output_dir)

        flat_results = [res for category in results.values() for res in category.values()]
        total = len(flat_results)
        successful = sum(1 for res in flat_results if res.success)
        reason_counts = Counter(
            (str(res.error_message).strip() or "unknown failure")
            for res in flat_results
            if not res.success
        )

        return {
            'success': True,
            'total_analyzed': total,
            'successful': successful,
            'by_category': {cat: len(res) for cat, res in results.items()},
            'failed': total - successful,
            'failure_reasons': dict(reason_counts),
            'output_dir': str(output_dir)
        }
    else:
        results = analyzer.analyze_directory(poses_dir, output_dir)
        flat_results = list(results.values())
        reason_counts = Counter(
            (str(res.error_message).strip() or "unknown failure")
            for res in flat_results
            if not res.success
        )

        return {
            'success': True,
            'total_analyzed': len(results),
            'successful': sum(1 for r in results.values() if r.success),
            'failed': sum(1 for r in results.values() if not r.success),
            'failure_reasons': dict(reason_counts),
            'output_dir': str(output_dir)
        }


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python poseview_integration.py <poses_dir> [output_dir]")
        print("\nExample:")
        print("  python poseview_integration.py ./best_poses ./poseview_output")
        sys.exit(1)
    
    poses_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("./poseview_output")
    
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    results = run_poseview_analysis(poses_dir, output_dir)
    print(f"\nResults: {json.dumps(results, indent=2)}")
