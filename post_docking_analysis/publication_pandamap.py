"""
Publication-Quality PandaMap Integration for Post-Docking Analysis.

Enhanced PandaMap integration with high-quality figure settings for publications:
- High DPI (300-600)
- Publication-ready formats (PDF, SVG, PNG)
- Customizable styling
- Consistent color schemes
- Proper figure sizing
"""
import pandas as pd
from pathlib import Path
from typing import Any, List, Dict, Optional, Tuple
import subprocess
import json
import logging
import re
from collections import Counter

try:
    from .protein_naming import resolve_protein_display_name, format_protein_label, extract_pdb_code
    from .pandamap_runner import PandaMapRunner
except ImportError:
    from protein_naming import resolve_protein_display_name, format_protein_label, extract_pdb_code
    from pandamap_runner import PandaMapRunner


class PublicationPandaMapAnalyzer:
    """
    Publication-quality PandaMap analyzer with enhanced figure settings.
    """
    
    # Publication-quality default settings
    PUBLICATION_CONFIG = {
        'dpi': 300,  # High resolution for print
        'formats': ['pdf', 'svg', 'png'],  # Multiple formats
        'figure_width': 10,  # inches
        'figure_height': 8,  # inches
        'font_size': 12,
        'font_family': 'Arial',
        'line_width': 2.0,
        'marker_size': 8,
        'color_scheme': 'publication',  # Consistent color scheme
        'background': 'white',
        'transparent': False,
        'bbox_inches': 'tight',
        'pad_inches': 0.1,
        # 3D quality/cue controls
        'show_surface': True,
        'show_3d_cues': True,
        # 2D map readability controls (from PandaMap CLI)
        'color_by': 'interaction',
        'size_scale': 1.35,
        # Optional per-project naming map for human-readable filenames/titles
        'protein_name_map': {},
        # Optional metadata from pipeline: {complex_name: {protein_label, ligand_display, affinity_category_label, best_affinity}}
        'complex_metadata': {},
        # When False, reuse already-generated PandaMap files in the output folder.
        'overwrite': False,
        # Enrich publication bundle with sidecar text reports when supported by the installed CLI.
        'generate_text_reports': True,
        # Include empirical ΔG estimate in report generation when supported.
        'estimate_delta_g': True,
        # Persist runtime option/capability snapshot alongside generated figures.
        'write_capability_manifest': True,
    }
    AMINO_ACID_CODES = {
        "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX",
        "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "SEC", "PYL",
    }
    
    def __init__(
        self,
        conda_env: str = "pandamap",
        config: Dict = None,
        publication_mode: bool = True
    ):
        """
        Initialize publication-quality PandaMap analyzer.
        
        Parameters
        ----------
        conda_env : str
            Conda environment name where PandaMap is installed
        config : Dict, optional
            Configuration dictionary (will merge with PUBLICATION_CONFIG)
        publication_mode : bool
            Enable publication-quality settings
        """
        self.conda_env = conda_env
        self.publication_mode = publication_mode
        
        # Merge user config with publication defaults.
        incoming = config or {}
        nested_cfg = incoming.get("visualization", {}).get("pandamap", {}) if isinstance(incoming, dict) else {}
        merged_cfg = {}
        if isinstance(incoming, dict):
            merged_cfg.update(incoming)
        merged_cfg.update(nested_cfg)
        self.config = {**self.PUBLICATION_CONFIG, **merged_cfg}
        
        self.logger = logging.getLogger(__name__)
        self.runner = PandaMapRunner(conda_env=conda_env, logger=self.logger)

    @staticmethod
    def _slugify(value: str) -> str:
        """Convert label text into a filesystem-safe slug."""
        slug = re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")
        return slug or "unknown"

    def _alias_map_name(self, map_name: str) -> str:
        """
        Replace PDB-like fragments in a filename with mapped protein display names.
        """
        protein_name_map = self.config.get('protein_name_map') or {}
        pdb_code = extract_pdb_code(str(map_name))
        if not pdb_code:
            return map_name

        display = resolve_protein_display_name(pdb_code, protein_name_map)
        if not display:
            return map_name

        display_slug = self._slugify(display)
        map_name_text = str(map_name)
        leading_token = self._slugify(map_name_text.split("__", 1)[0])
        if "__" in map_name_text and leading_token and (
            leading_token.lower() == display_slug.lower()
            or leading_token.lower() in display_slug.lower()
            or display_slug.lower() in leading_token.lower()
        ):
            return map_name
        if str(map_name).lower().startswith(f"{display_slug.lower()}__"):
            return map_name
        # Keep original stem for traceability while making the protein explicit.
        return f"{display_slug}__{map_name}"

    def _get_complex_metadata(self, complex_key: str) -> Dict[str, object]:
        """
        Fetch metadata for a complex stem if provided by the pipeline.
        """
        metadata = self.config.get("complex_metadata") or {}
        if not isinstance(metadata, dict):
            return {}
        value = metadata.get(str(complex_key), {})
        return value if isinstance(value, dict) else {}

    def _build_map_name(self, complex_key: str) -> str:
        """
        Build readable output stem that includes protein, ligand, and category.
        """
        meta = self._get_complex_metadata(complex_key)
        if not meta:
            return str(complex_key)

        protein = self._slugify(str(meta.get("protein_label") or "Protein"))
        ligand = self._slugify(str(meta.get("ligand_display") or "Ligand"))
        category = self._slugify(str(meta.get("affinity_category_label") or "Category"))
        label = f"{protein}__{ligand}__{category}__{complex_key}"
        if len(label) > 180:
            label = label[:180].rstrip("_")
        return label

    def _build_title(self, complex_key: str, ligand_name: str) -> str:
        """
        Build an informative title for PandaMap figures.
        """
        meta = self._get_complex_metadata(complex_key)
        if not meta:
            protein_text = format_protein_label(
                complex_key,
                self.config.get("protein_name_map")
            )
            return f"{protein_text} | ligand {ligand_name}"

        protein_label = str(meta.get("protein_label") or complex_key)
        ligand_label = str(meta.get("ligand_display") or ligand_name)
        category = str(meta.get("affinity_category_label") or "Uncategorized")
        affinity = meta.get("best_affinity")
        if affinity is None:
            return f"{protein_label} | Ligand {ligand_label} | {category}"

        try:
            affinity_value = float(affinity)
            return (
                f"{protein_label} | Ligand {ligand_label} | "
                f"{category} ({affinity_value:.2f} kcal/mol)"
            )
        except Exception:
            return f"{protein_label} | Ligand {ligand_label} | {category}"

    def _balanced_subset(self, pdb_files: List[Path], max_complexes: Optional[int]) -> List[Path]:
        """
        Select complexes in a protein-balanced order instead of alphabetical-only.
        """
        if not max_complexes or max_complexes <= 0 or len(pdb_files) <= max_complexes:
            return pdb_files

        buckets: Dict[str, List[Path]] = {}
        for pdb_file in pdb_files:
            meta = self._get_complex_metadata(pdb_file.stem)
            protein_key = str(meta.get("protein_label") or pdb_file.stem)
            buckets.setdefault(protein_key, []).append(pdb_file)

        # Keep deterministic output.
        for key in buckets:
            buckets[key] = sorted(buckets[key], key=lambda p: p.name)

        ordered_keys = sorted(buckets.keys())
        selected: List[Path] = []
        while len(selected) < max_complexes:
            added = 0
            for key in ordered_keys:
                if not buckets[key]:
                    continue
                selected.append(buckets[key].pop(0))
                added += 1
                if len(selected) >= max_complexes:
                    break
            if added == 0:
                break
        return selected

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

        for code in normalized:
            if code not in PublicationPandaMapAnalyzer.AMINO_ACID_CODES:
                return code
        return normalized[0] if normalized else "UNK"

    @staticmethod
    def _is_generic_ligand_name(value: str) -> bool:
        return str(value or "").strip().upper() in {"", "UNK", "UNX", "LIG", "UNL"}

    @classmethod
    def _should_pass_explicit_ligand(cls, value: str) -> bool:
        """
        Passing amino-acid-like residue names can misidentify protein residues as ligands.
        """
        token = str(value or "").strip().upper()
        if cls._is_generic_ligand_name(token):
            return False
        if token in cls.AMINO_ACID_CODES:
            return False
        return True

    def _infer_ligand_from_name(self, pdb_file: Path) -> Optional[str]:
        """
        Infer ligand residue name from complex naming conventions.
        """
        text = pdb_file.stem
        match = re.search(r"_ligand_([A-Za-z0-9]{1,6})_([A-Za-z])_(\d+)", text)
        if match:
            return self._normalize_resname(match.group(1))

        token = text.split("_")[-1] if text else ""
        if token.lower() in {"pdbqt", "poses", "top", "out"} and "_" in text:
            for part in reversed(text.split("_")):
                if part.lower() not in {"pdbqt", "poses", "top", "out"}:
                    token = part
                    break

        resname = self._normalize_resname(token)
        if self._is_generic_ligand_name(resname):
            return None
        return resname

    def _detect_ligand_from_pdb(self, pdb_file: Path) -> Optional[str]:
        """
        Detect most likely ligand residue name from HETATM records.
        """
        exclude = {"HOH", "WAT", "NA", "CL", "K", "CA", "MG", "MN", "ZN", "SO4", "PO4"}
        counts: Counter[str] = Counter()
        try:
            with pdb_file.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if not line.startswith("HETATM"):
                        continue
                    resname = self._normalize_resname(line[17:20].strip())
                    if resname in exclude or self._is_generic_ligand_name(resname):
                        continue
                    counts[resname] += 1
        except Exception as exc:
            self.logger.debug(f"PandaMap ligand detection failed for {pdb_file.name}: {exc}")
            return None

        if not counts:
            return None
        return counts.most_common(1)[0][0]

    def _resolve_ligand_name(self, pdb_file: Path, default_ligand: str) -> str:
        """
        Resolve ligand residue name per complex for PandaMap CLI.
        """
        ligand, _ = self._resolve_ligand_name_with_source(pdb_file, default_ligand)
        return ligand

    def _resolve_ligand_name_with_source(self, pdb_file: Path, default_ligand: str) -> Tuple[str, str]:
        """
        Resolve ligand residue name and provenance for diagnostics.
        """
        normalized_default = self._normalize_resname(default_ligand)
        if (
            not self._is_generic_ligand_name(normalized_default)
            and self._should_pass_explicit_ligand(normalized_default)
        ):
            return normalized_default, "pipeline_default"

        # Prefer residue identity found in the actual structure over filename heuristics.
        detected = self._detect_ligand_from_pdb(pdb_file)
        if detected and not self._is_generic_ligand_name(detected):
            return detected, "detected_from_pdb_hetatm"

        inferred = self._infer_ligand_from_name(pdb_file)
        if inferred and not self._is_generic_ligand_name(inferred):
            return inferred, "inferred_from_filename"

        return normalized_default, "generic_fallback"

    @staticmethod
    def _documented_feature_reference() -> Dict[str, Any]:
        """
        Static, version-agnostic PandaMap feature baseline from upstream docs.
        """
        return {
            "interaction_coverage": [
                "hydrogen_bond",
                "pi_pi_stacking",
                "cation_pi",
                "pi_cation",
                "carbon_pi",
                "donor_pi",
                "amide_pi",
                "alkyl_pi",
                "hydrophobic",
                "ionic",
                "salt_bridge",
                "halogen_bond",
                "metal_coordination",
                "covalent",
                "attractive_charge",
                "repulsive_charge",
            ],
            "outputs": ["2d_png", "3d_html", "text_report"],
            "analysis_modes": ["single_structure", "trajectory_occupancy", "empirical_delta_g"],
            "input_formats": ["pdb", "mmcif/cif", "pdbqt"],
            "optional_dependencies": {
                "rdkit": "chemically accurate 2D coordinates",
                "dssp": "preferred solvent accessibility calculation",
                "py3dmol_extra": "programmatic 3D viewer integration",
            },
            "notes": [
                "Exact option support is runtime-version dependent; check pandamap_capabilities.json for this run.",
                "Interaction cutoffs follow the PandaMap documentation baseline (v4.2 series).",
            ],
        }

    def generate_publication_2d_map(
        self,
        pdb_file: Path,
        ligand_name: str = "UNK",
        output_dir: Path = None,
        map_name: str = None,
        dpi: int = None,
        formats: List[str] = None
    ) -> Dict[str, object]:
        """
        Generate publication-quality 2D interaction map.
        
        Parameters
        ----------
        pdb_file : Path
            Path to the PDB file containing the complex
        ligand_name : str
            Name of the ligand residue
        output_dir : Path, optional
            Output directory for the map
        map_name : str, optional
            Name for the map file
        dpi : int, optional
            Resolution (default: from config, typically 300)
        formats : List[str], optional
            Output formats (default: ['pdf', 'svg', 'png'])
            
        Returns
        -------
        Dict[str, object]
            Payload containing generated files, optional text report path, and
            the command arguments used per output format.
        """
        self.logger.info(f"🐼 Generating publication-quality 2D interaction map for {pdb_file.name}...")
        
        if output_dir is None:
            output_dir = pdb_file.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if map_name is None or map_name == pdb_file.stem:
            map_name = self._build_map_name(pdb_file.stem)
            if map_name == pdb_file.stem:
                map_name = self._alias_map_name(map_name)
        else:
            map_name = self._alias_map_name(map_name)
        
        dpi = dpi or self.config['dpi']
        formats = formats or self.config['formats']
        
        generated_files: Dict[str, Path] = {}
        run_commands: Dict[str, List[str]] = {}
        report_file: Optional[Path] = None
        cli_mode = self.runner.detect_cli_mode()
        supported_options = self.runner.get_supported_options()
        overwrite_outputs = bool(self.config.get("overwrite", False) or self.config.get("force_regenerate", False))
        
        for fmt in formats:
            output_file = output_dir / f"{map_name}.{fmt}"
            if output_file.exists() and not overwrite_outputs:
                generated_files[fmt] = output_file
                self.logger.info(f"  ↪️ Reusing existing {fmt.upper()}: {output_file.name}")
                continue
            
            if cli_mode == "single":
                # New PandaMap CLI: positional structure_file, output format inferred from extension.
                title = self.config.get("title") or self._build_title(pdb_file.stem, ligand_name)
                args = [str(pdb_file)]
                if "--ligand" in supported_options and self._should_pass_explicit_ligand(ligand_name):
                    args.extend(["--ligand", ligand_name])
                if "--output" in supported_options:
                    args.extend(["--output", str(output_file)])
                if "--dpi" in supported_options:
                    args.extend(["--dpi", str(dpi)])
                if "--title" in supported_options and str(title).strip():
                    args.extend(["--title", str(title)])
                if bool(self.config.get("generate_text_reports", True)) and "--report" in supported_options:
                    args.append("--report")
                    if "--report-file" in supported_options:
                        report_file = output_dir / f"{map_name}.txt"
                        args.extend(["--report-file", str(report_file)])
                if bool(self.config.get("estimate_delta_g", True)) and "--deltaG" in supported_options:
                    args.append("--deltaG")
                if not self.config.get("show_3d_cues", True) and "--no-3d-cues" in supported_options:
                    args.append("--no-3d-cues")
            else:
                # Legacy PandaMap CLI (generate subcommand).
                args = [
                    "generate",
                    "--input", str(pdb_file),
                    "--ligand", ligand_name,
                    "--output", str(output_file),
                    "--format", fmt,
                    "--dpi", str(dpi),
                    "--width", str(self.config['figure_width']),
                    "--height", str(self.config['figure_height']),
                    "--font-size", str(self.config['font_size']),
                    "--font-family", self.config['font_family'],
                    "--line-width", str(self.config['line_width']),
                    "--background", self.config['background']
                ]
                if self.config.get('transparent', False) and fmt in ['png', 'svg']:
                    args.extend(["--transparent"])
                if bool(self.config.get("generate_text_reports", True)) and "--report" in supported_options:
                    args.append("--report")
                if bool(self.config.get("estimate_delta_g", True)) and "--deltaG" in supported_options:
                    args.append("--deltaG")
            run_commands[fmt] = list(args)
            
            try:
                result = self.runner.run(args, timeout=300, cwd=output_dir)
                
                if result.returncode == 0 and output_file.exists():
                    generated_files[fmt] = output_file
                    self.logger.info(f"  ✅ Generated {fmt.upper()}: {output_file.name}")
                else:
                    self.logger.warning(f"  ⚠️  Failed to generate {fmt}: {result.stderr}")
                    
            except subprocess.TimeoutExpired:
                self.logger.warning(f"  ⚠️  Generation timed out for {fmt}")
            except FileNotFoundError:
                self.logger.error("  ❌ PandaMap not found. Please install PandaMap.")
                break
            except Exception as e:
                self.logger.warning(f"  ⚠️  Error generating {fmt}: {e}")
        
        if generated_files:
            self.logger.info(f"✅ Generated {len(generated_files)} format(s) for {map_name}")
        else:
            self.logger.error(f"❌ Failed to generate any formats for {map_name}")

        return {
            "files": generated_files,
            "report_file": report_file if report_file and report_file.exists() else None,
            "commands": run_commands,
        }
    
    def generate_publication_3d_visualization(
        self,
        pdb_file: Path,
        ligand_name: str = "UNK",
        output_dir: Path = None,
        vis_name: str = None,
        interactive: bool = True
    ) -> Optional[Path]:
        """
        Generate publication-quality 3D visualization.
        
        Parameters
        ----------
        pdb_file : Path
            Path to the PDB file containing the complex
        ligand_name : str
            Name of the ligand residue
        output_dir : Path, optional
            Output directory for the visualization
        vis_name : str, optional
            Name for the visualization file
        interactive : bool
            Generate interactive HTML visualization
            
        Returns
        -------
        Path or None
            Path to generated visualization file
        """
        self.logger.info(f"🌐 Generating publication-quality 3D visualization for {pdb_file.name}...")
        
        if output_dir is None:
            output_dir = pdb_file.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if vis_name is None or vis_name == pdb_file.stem:
            vis_name = self._build_map_name(pdb_file.stem)
            if vis_name == pdb_file.stem:
                vis_name = self._alias_map_name(vis_name)
        else:
            vis_name = self._alias_map_name(vis_name)
        cli_mode = self.runner.detect_cli_mode()
        supported_options = self.runner.get_supported_options()
        overwrite_outputs = bool(self.config.get("overwrite", False) or self.config.get("force_regenerate", False))

        if interactive:
            output_file = output_dir / f"{vis_name}.html"
        else:
            output_file = output_dir / f"{vis_name}.png"
        if output_file.exists() and not overwrite_outputs:
            self.logger.info(f"↪️ Reusing existing 3D visualization: {output_file.name}")
            return output_file

        if cli_mode == "single":
            if "--3d" not in supported_options:
                self.logger.info("ℹ️  Installed PandaMap does not support --3d; skipping 3D visualization")
                return None
            if "--3d-output" not in supported_options and "--output" not in supported_options:
                self.logger.info("ℹ️  Installed PandaMap does not expose 3D output arguments; skipping 3D visualization")
                return None
            # New PandaMap CLI: enable 3D generation via --3d / --3d-output.
            title = self.config.get("title") or self._build_title(pdb_file.stem, ligand_name)
            args = [str(pdb_file), "--3d"]
            if "--ligand" in supported_options and self._should_pass_explicit_ligand(ligand_name):
                args.extend(["--ligand", ligand_name])
            if "--3d-output" in supported_options:
                args.extend(["--3d-output", str(output_file)])
            elif "--output" in supported_options:
                args.extend(["--output", str(output_file)])
            if "--title" in supported_options and str(title).strip():
                args.extend(["--title", str(title)])
            if "--width" in supported_options:
                args.extend(["--width", str(int(self.config['figure_width'] * 100))])
            if "--height" in supported_options:
                args.extend(["--height", str(int(self.config['figure_height'] * 100))])
            if not self.config.get("show_surface", True) and "--no-surface" in supported_options:
                args.append("--no-surface")
            if not self.config.get("show_3d_cues", True) and "--no-3d-cues" in supported_options:
                args.append("--no-3d-cues")
            if not interactive and "--dpi" in supported_options:
                args.extend(["--dpi", str(self.config['dpi'])])
        else:
            if interactive:
                args = [
                    "visualize",
                    "--input", str(pdb_file),
                    "--ligand", ligand_name,
                    "--output", str(output_file),
                    "--format", "html",
                    "--quality", "high"
                ]
            else:
                args = [
                    "visualize",
                    "--input", str(pdb_file),
                    "--ligand", ligand_name,
                    "--output", str(output_file),
                    "--format", "png",
                    "--dpi", str(self.config['dpi']),
                    "--width", str(self.config['figure_width']),
                    "--height", str(self.config['figure_height'])
                ]
        
        try:
            result = self.runner.run(args, timeout=300, cwd=output_dir)
            
            if result.returncode == 0 and output_file.exists():
                self.logger.info(f"✅ 3D visualization generated: {output_file.name}")
                return output_file
            else:
                self.logger.warning(f"⚠️  PandaMap execution failed: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            self.logger.warning("⚠️  PandaMap execution timed out")
            return None
        except FileNotFoundError:
            self.logger.error("⚠️  PandaMap not found. Please install PandaMap.")
            return None
        except Exception as e:
            self.logger.warning(f"⚠️  Error generating 3D visualization: {e}")
            return None
    
    def generate_comprehensive_publication_analysis(
        self,
        complexes_dir: Path,
        output_dir: Path,
        ligand_name: str = "UNK",
        max_complexes: int = None,
        generate_2d: bool = True,
        generate_3d: bool = True
    ) -> Dict:
        """
        Generate comprehensive publication-quality interaction analysis.
        
        Parameters
        ----------
        complexes_dir : Path
            Directory containing PDB complex files
        output_dir : Path
            Output directory for analysis results
        ligand_name : str
            Name of the ligand residue
        max_complexes : int, optional
            Maximum number of complexes to analyze (for performance)
        generate_2d : bool
            Generate 2D interaction maps
        generate_3d : bool
            Generate 3D visualizations
            
        Returns
        -------
        Dict
            Dictionary containing analysis summary
        """
        self.logger.info("🐼 Generating comprehensive publication-quality PandaMap analysis...")
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find all PDB files
        pdb_files = sorted(list(complexes_dir.glob("*.pdb")))
        
        if not pdb_files:
            self.logger.warning("⚠️  No PDB files found for analysis")
            return {}
        
        # Limit to max_complexes if specified
        pdb_files = self._balanced_subset(pdb_files, max_complexes)
        
        self.logger.info(f"📊 Analyzing {len(pdb_files)} complexes...")
        
        # Create output subdirectories
        maps_2d_dir = output_dir / "maps_2d" if generate_2d else None
        vis_3d_dir = output_dir / "maps_3d" if generate_3d else None
        
        if maps_2d_dir:
            maps_2d_dir.mkdir(exist_ok=True)
        if vis_3d_dir:
            vis_3d_dir.mkdir(exist_ok=True)
        
        generated_2d_maps = 0
        generated_3d_visualizations = 0
        analysis_results = []
        capabilities = self.runner.describe_capabilities()
        feature_reference = self._documented_feature_reference()
        
        # Generate analysis for each complex
        for i, pdb_file in enumerate(pdb_files, 1):
            self.logger.info(f"  [{i}/{len(pdb_files)}] Processing {pdb_file.name}...")
            
            complex_name = pdb_file.stem
            complex_meta = self._get_complex_metadata(complex_name)
            complex_ligand, ligand_source = self._resolve_ligand_name_with_source(pdb_file, ligand_name)
            if self._is_generic_ligand_name(complex_ligand):
                self.logger.warning(
                    f"⚠️  Could not resolve ligand for {pdb_file.name}; keeping fallback {complex_ligand}"
                )
            title = self._build_title(complex_name, complex_ligand)
            map_name = self._build_map_name(complex_name)
            
            # Generate 2D interaction map
            if generate_2d:
                maps_payload = self.generate_publication_2d_map(
                    pdb_file,
                    complex_ligand,
                    maps_2d_dir,
                    map_name=complex_name
                )
                maps = maps_payload.get("files", {}) if isinstance(maps_payload, dict) else {}
                report_file = maps_payload.get("report_file") if isinstance(maps_payload, dict) else None
                commands = maps_payload.get("commands", {}) if isinstance(maps_payload, dict) else {}
                if maps:
                    generated_2d_maps += 1
                    analysis_results.append({
                        'complex': complex_name,
                        'protein_label': complex_meta.get("protein_label"),
                        'protein_name': complex_meta.get("protein_name"),
                        'pdb_code': complex_meta.get("pdb_code"),
                        'ligand_display': complex_meta.get("ligand_display"),
                        'affinity_category': complex_meta.get("affinity_category"),
                        'affinity_category_label': complex_meta.get("affinity_category_label"),
                        'best_affinity': complex_meta.get("best_affinity"),
                        'ligand': complex_ligand,
                        'ligand_resolution_source': ligand_source,
                        'figure_title': title,
                        'map_stem': map_name,
                        '2d_maps': list(maps.keys()),
                        'report_file': str(report_file) if report_file else "",
                        'commands': commands,
                        '3d_vis': False
                    })
            
            # Generate 3D visualization
            if generate_3d:
                vis_file = self.generate_publication_3d_visualization(
                    pdb_file,
                    complex_ligand,
                    vis_3d_dir,
                    vis_name=complex_name,
                    interactive=True
                )
                if vis_file:
                    generated_3d_visualizations += 1
                    if analysis_results:
                        analysis_results[-1]['3d_vis'] = True
                    else:
                        analysis_results.append({
                            'complex': complex_name,
                            'protein_label': complex_meta.get("protein_label"),
                            'protein_name': complex_meta.get("protein_name"),
                            'pdb_code': complex_meta.get("pdb_code"),
                            'ligand_display': complex_meta.get("ligand_display"),
                            'affinity_category': complex_meta.get("affinity_category"),
                            'affinity_category_label': complex_meta.get("affinity_category_label"),
                            'best_affinity': complex_meta.get("best_affinity"),
                            'ligand': complex_ligand,
                            'ligand_resolution_source': ligand_source,
                            'figure_title': title,
                            'map_stem': map_name,
                            '2d_maps': [],
                            'report_file': "",
                            'commands': {},
                            '3d_vis': True
                        })
        
        # Generate summary report
        summary = {
            'total_poses_analyzed': len(pdb_files),
            'generated_2d_maps': generated_2d_maps,
            'generated_3d_visualizations': generated_3d_visualizations,
            'analysis_results': analysis_results,
            'publication_settings': {
                'dpi': self.config['dpi'],
                'formats': self.config['formats'],
                'figure_size': f"{self.config['figure_width']}x{self.config['figure_height']} inches"
            },
            'pandamap_runtime_capabilities': capabilities,
            'pandamap_documented_feature_reference': feature_reference,
            'analysis_timestamp': pd.Timestamp.now().isoformat()
        }
        
        # Save summary to JSON
        summary_file = output_dir / "pandamap_publication_analysis_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Save detailed results to CSV
        if analysis_results:
            results_df = pd.DataFrame(analysis_results)
            if "commands" in results_df.columns:
                results_df["commands"] = results_df["commands"].apply(
                    lambda value: json.dumps(value, ensure_ascii=False) if isinstance(value, dict) else value
                )
            results_csv = output_dir / "pandamap_analysis_results.csv"
            results_df.to_csv(results_csv, index=False)
        if bool(self.config.get("write_capability_manifest", True)):
            capability_file = output_dir / "pandamap_capabilities.json"
            with capability_file.open("w", encoding="utf-8") as handle:
                json.dump(
                    {
                        "runtime": capabilities,
                        "documented_feature_reference": feature_reference,
                    },
                    handle,
                    indent=2,
                )
            summary["capability_manifest"] = str(capability_file)
        
        self.logger.info("✅ Comprehensive publication-quality PandaMap analysis completed")
        self.logger.info(f"   📊 Generated {generated_2d_maps} 2D interaction maps")
        self.logger.info(f"   🌐 Generated {generated_3d_visualizations} 3D visualizations")
        self.logger.info(f"   📄 Summary saved to: {summary_file}")
        
        return summary


def run_publication_pandamap_analysis(
    complexes_dir: Path,
    output_dir: Path,
    ligand_name: str = "UNK",
    conda_env: str = "pandamap",
    config: Dict = None,
    max_complexes: int = None
) -> Dict:
    """
    Run comprehensive publication-quality PandaMap analysis on complexes.
    
    Parameters
    ----------
    complexes_dir : Path
        Directory containing PDB complex files
    output_dir : Path
        Output directory for analysis results
    ligand_name : str
        Name of the ligand residue
    conda_env : str
        Conda environment name where PandaMap is installed
    config : Dict, optional
        Configuration dictionary for publication settings
    max_complexes : int, optional
        Maximum number of complexes to analyze
        
    Returns
    -------
    Dict
        Dictionary containing analysis results
    """
    logger = logging.getLogger(__name__)
    logger.info("🐼 Running publication-quality PandaMap interaction analysis...")
    
    # Initialize analyzer
    analyzer = PublicationPandaMapAnalyzer(conda_env, config, publication_mode=True)
    
    # Generate comprehensive analysis
    summary = analyzer.generate_comprehensive_publication_analysis(
        complexes_dir=complexes_dir,
        output_dir=output_dir,
        ligand_name=ligand_name,
        max_complexes=max_complexes
    )
    
    return summary
