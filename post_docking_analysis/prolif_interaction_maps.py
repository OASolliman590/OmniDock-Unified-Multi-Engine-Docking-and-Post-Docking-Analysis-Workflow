"""
ProLIF 2.x interaction map generation for protein-ligand complexes.

Creates:
- static barcode PNG (matplotlib)
- optional LigNetwork HTML (interactive)
- flattened interaction-table CSV
"""
from pathlib import Path
from typing import Dict, Optional, Tuple, List, Any
import logging
import os
import tempfile
import re

try:
    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    import prolif as plf
    from MDAnalysis.topology import tables as mda_tables
    PROLIF_AVAILABLE = True
except ImportError:
    PROLIF_AVAILABLE = False

logger = logging.getLogger(__name__)

_VDW_PATCH = {
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
    "Zn": 1.39,
    "Mg": 1.73,
    "Na": 2.27,
    "K": 2.75,
    "Ca": 2.31,
    "Fe": 1.56,
    "Cu": 1.40,
    "Mn": 1.61,
}

_VALID_ELEMENTS = {
    "H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
    "Na", "K", "Ca", "Mg", "Zn", "Mn", "Fe", "Cu", "Co", "Ni",
}


def _slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(value or "")).strip("_")
    return text or "unknown"


def _metadata_output_stem(
    complex_name: str,
    complex_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
) -> str:
    if not complex_metadata:
        return f"{complex_name}_interaction_map"
    meta = complex_metadata.get(str(complex_name), {})
    if not isinstance(meta, dict) or not meta:
        return f"{complex_name}_interaction_map"

    protein_slug = _slugify(str(meta.get("protein_label") or "Protein"))
    ligand_slug = _slugify(str(meta.get("ligand_display") or "Ligand"))
    category_slug = _slugify(str(meta.get("affinity_category_label") or "Uncategorized"))
    stem = f"{protein_slug}__{ligand_slug}__{category_slug}__{complex_name}_interaction_map"
    if len(stem) > 190:
        stem = stem[:190].rstrip("_")
    return stem


def _metadata_display_title(
    complex_name: str,
    ligand_resname: str,
    complex_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
) -> str:
    if not complex_metadata:
        return f"ProLIF Docking IFP ({complex_name}) | ligand={ligand_resname}"
    meta = complex_metadata.get(str(complex_name), {})
    if not isinstance(meta, dict) or not meta:
        return f"ProLIF Docking IFP ({complex_name}) | ligand={ligand_resname}"

    protein_label = str(meta.get("protein_label") or complex_name)
    ligand_label = str(meta.get("ligand_display") or ligand_resname)
    category_label = str(meta.get("affinity_category_label") or "Uncategorized")
    affinity = meta.get("best_affinity")
    if affinity is None:
        return f"ProLIF Docking IFP: {protein_label} | Ligand {ligand_label} | {category_label}"
    try:
        return (
            f"ProLIF Docking IFP: {protein_label} | Ligand {ligand_label} | "
            f"{category_label} ({float(affinity):.2f} kcal/mol)"
        )
    except Exception:
        return f"ProLIF Docking IFP: {protein_label} | Ligand {ligand_label} | {category_label}"


def _patch_mda_vdwradii() -> None:
    """Extend MDAnalysis vdw table for common halogens/metals in docking structures."""
    for key, value in _VDW_PATCH.items():
        mda_tables.vdwradii[key] = value


def _infer_element(atom_name: str, atom_type: str) -> str:
    """Infer a valid element symbol from atom name/type strings."""
    atype = str(atom_type or "").strip()
    up_type = atype.upper()
    low_type = atype.lower()

    if low_type in {"cl", "br", "na", "mg", "zn", "mn", "fe", "cu", "co", "ni", "cd", "hg"}:
        return low_type.capitalize()
    if low_type == "ca":
        # In docking files CA may represent aromatic carbon.
        if atype == "Ca":
            return "Ca"
        return "C"
    if up_type in {"A", "C", "CA", "CG0", "CG1", "CG2"}:
        return "C"
    if up_type.startswith("N"):
        return "N"
    if up_type.startswith("O"):
        return "O"
    if up_type.startswith("S"):
        return "S"
    if up_type.startswith("P"):
        return "P"
    if up_type.startswith("H"):
        return "H"
    if up_type in {"F", "I"}:
        return up_type

    name = "".join(ch for ch in str(atom_name).strip() if ch.isalpha()).upper()
    if name.startswith("CL"):
        return "Cl"
    if name.startswith("BR"):
        return "Br"
    if name.startswith("ZN"):
        return "Zn"
    if name.startswith("MG"):
        return "Mg"
    if name.startswith("MN"):
        return "Mn"
    if name.startswith("FE"):
        return "Fe"
    if name.startswith("CU"):
        return "Cu"
    if name.startswith("CO"):
        return "Co"
    if name.startswith("NI"):
        return "Ni"
    if name.startswith("NA"):
        return "Na"
    if name.startswith("CD"):
        return "Cd"
    if name.startswith("HG"):
        return "Hg"
    if name:
        return name[0]
    return ""


def _sanitize_pdb_elements(input_pdb: Path) -> Path:
    """
    Ensure ATOM/HETATM rows contain valid element symbols (cols 77-78).
    Returns original file if no edits are needed.
    """
    fixed_lines: List[str] = []
    changed = False

    with input_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                fixed_lines.append(line)
                continue

            raw = line.rstrip("\n")
            newline = "\n" if line.endswith("\n") else ""
            if len(raw) < 80:
                raw = raw.ljust(80)

            current = raw[76:78].strip()
            current_norm = current.capitalize() if current else ""
            if current_norm not in _VALID_ELEMENTS:
                atom_name = raw[12:16]
                atom_type = raw.split()[-1] if raw.split() else ""
                inferred = _infer_element(atom_name, atom_type)
                if inferred:
                    raw = raw[:76] + f"{inferred:>2}" + raw[78:]
                    changed = True

            fixed_lines.append(raw + newline)

    if not changed:
        return input_pdb

    fd, tmp_name = tempfile.mkstemp(prefix=f"{input_pdb.stem}_elements_", suffix=".pdb")
    tmp_path = Path(tmp_name)
    os.close(fd)
    with tmp_path.open("w", encoding="utf-8") as out:
        out.writelines(fixed_lines)
    return tmp_path


class ProLifInteractionMapper:
    """Generate ProLIF barcode/network outputs for a single complex file."""

    def __init__(self):
        if not PROLIF_AVAILABLE:
            logger.warning(
                "ProLIF stack not available. Install with:\n"
                "conda install -c conda-forge prolif mdanalysis rdkit"
            )

    @staticmethod
    def _detect_ligand_resname(universe: "mda.Universe", preferred: str = "UNK") -> Optional[str]:
        # Prefer explicit ligand name when present.
        preferred_atoms = universe.select_atoms(f"resname {preferred}")
        if len(preferred_atoms) > 0:
            return preferred

        # Fallback: use largest non-protein hetero residue.
        residue_sizes: Dict[str, int] = {}
        hetero = universe.select_atoms("not protein and not resname HOH WAT SOL")
        for residue in hetero.residues:
            residue_sizes[residue.resname] = residue_sizes.get(residue.resname, 0) + len(residue.atoms)
        if not residue_sizes:
            return None
        return sorted(residue_sizes.items(), key=lambda item: item[1], reverse=True)[0][0]

    @staticmethod
    def _save_lignetwork_html(html_obj: Any, output_html: Path) -> bool:
        try:
            output_html.parent.mkdir(parents=True, exist_ok=True)
            payload = getattr(html_obj, "data", None)
            if not payload and hasattr(html_obj, "_repr_html_"):
                payload = html_obj._repr_html_()
            if not payload:
                return False
            output_html.write_text(str(payload), encoding="utf-8")
            return True
        except Exception as exc:
            logger.warning(f"⚠️  Could not write ProLIF LigNetwork HTML: {exc}")
            return False

    @staticmethod
    def _save_interaction_table(df, output_csv: Path) -> bool:
        try:
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            flat_df = df.copy()
            flat_df.index = [str(v) for v in flat_df.index]
            flat_df.columns = [
                "|".join([str(x) for x in col]) if isinstance(col, tuple) else str(col)
                for col in flat_df.columns
            ]
            flat_df.to_csv(output_csv)
            return True
        except Exception as exc:
            logger.warning(f"⚠️  Could not write ProLIF interaction table: {exc}")
            return False

    @staticmethod
    def _save_placeholder_png(output_png: Path, title: str, message: str, dpi: int, figsize: Tuple[int, int]) -> None:
        plt.close("all")  # ensure no stale renderer from a previous failed plot
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.axis("off")
        output_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def _normalize_resname(raw_name: str) -> str:
        token = "".join(ch for ch in str(raw_name).upper() if ch.isalnum())
        if not token:
            return "UNK"
        if len(token) >= 3:
            return token[:3]
        if len(token) == 2:
            return token + "X"
        return token + "XX"

    def _infer_ligand_resname_from_name(self, name: str) -> Optional[str]:
        """
        Infer ligand token from common docking naming patterns.
        """
        text = str(name or "")
        match = re.search(r"_ligand_([A-Za-z0-9]{1,6})_([A-Za-z])_(\d+)", text)
        if match:
            return self._normalize_resname(match.group(1))

        stem = Path(text).stem
        token = stem.split("_")[-1] if stem else ""
        if token.lower() in {"pdbqt", "poses", "top", "out"} and "_" in stem:
            for part in reversed(stem.split("_")):
                if part.lower() not in {"pdbqt", "poses", "top", "out"}:
                    token = part
                    break
        resname = self._normalize_resname(token)
        if resname in {"UNK", "UNX", "LIG"}:
            return None
        return resname

    def _load_pose_iterable_from_sdf(
        self,
        poses_sdf: Path,
        max_poses: Optional[int] = None,
    ) -> List[Any]:
        """
        Load docking poses from SDF as ProLIF-compatible molecules.
        """
        pose_iterable: List[Any] = []
        if hasattr(plf, "sdf_supplier"):
            supplier = plf.sdf_supplier(str(poses_sdf))
            for pose in supplier:
                if pose is None:
                    continue
                pose_iterable.append(pose)
                if max_poses and len(pose_iterable) >= max_poses:
                    break
            return pose_iterable

        # Fallback for environments where sdf_supplier is unavailable.
        try:
            from rdkit import Chem
        except Exception as exc:
            logger.warning(f"⚠️  RDKit not available for SDF fallback loading: {exc}")
            return []

        supplier = Chem.SDMolSupplier(str(poses_sdf), removeHs=False)
        for mol in supplier:
            if mol is None:
                continue
            pose_iterable.append(plf.Molecule.from_rdkit(mol))
            if max_poses and len(pose_iterable) >= max_poses:
                break
        return pose_iterable

    def _save_interaction_frequency(self, df, output_csv: Path) -> bool:
        """
        Save per-interaction frequencies across docking poses.
        """
        try:
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            freq = df.astype(float).mean(axis=0)
            freq_df = freq.rename("frequency").to_frame().reset_index()
            freq_df.columns = ["ligand", "protein", "interaction", "frequency"]
            freq_df.to_csv(output_csv, index=False)
            return True
        except Exception as exc:
            logger.warning(f"⚠️  Could not write ProLIF interaction frequency table: {exc}")
            return False

    def create_interaction_map(
        self,
        complex_pdb: Path,
        output_png: Path,
        ligand_resname: str = "UNK",
        dpi: int = 300,
        figsize: Tuple[int, int] = (12, 8),
        write_html: bool = True,
        display_title: Optional[str] = None,
    ) -> bool:
        if not PROLIF_AVAILABLE:
            logger.error("ProLIF stack not available")
            return False
        if not complex_pdb.exists():
            logger.error(f"Complex PDB file not found: {complex_pdb}")
            return False

        output_png.parent.mkdir(parents=True, exist_ok=True)
        _patch_mda_vdwradii()
        working_pdb = complex_pdb

        try:
            working_pdb = _sanitize_pdb_elements(complex_pdb)
            universe = mda.Universe(str(working_pdb))
            selected_resname = self._detect_ligand_resname(universe, ligand_resname)
            if not selected_resname:
                self._save_placeholder_png(
                    output_png,
                    title=display_title or f"ProLIF Interactions: {complex_pdb.stem}",
                    message="No ligand residue detected in complex",
                    dpi=dpi,
                    figsize=figsize,
                )
                return True

            ligand_atoms = universe.select_atoms(f"resname {selected_resname}")
            protein_atoms = universe.select_atoms("protein")
            if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
                self._save_placeholder_png(
                    output_png,
                    title=display_title or f"ProLIF Interactions: {complex_pdb.stem}",
                    message="Missing protein or ligand atoms for ProLIF",
                    dpi=dpi,
                    figsize=figsize,
                )
                return True

            ligand_mol = plf.Molecule.from_mda(ligand_atoms, NoImplicit=False)
            protein_mol = plf.Molecule.from_mda(protein_atoms, NoImplicit=False)

            fp = plf.Fingerprint()
            fp.run_from_iterable([ligand_mol], protein_mol)
            interaction_df = fp.to_dataframe()

            # Static barcode PNG
            # plot_barcode() crashes with KeyError:'ligand' on an empty fingerprint.
            # Guard: only call when interactions exist; always close all figures first.
            try:
                has_interactions = (
                    not interaction_df.empty and len(interaction_df.columns) > 0
                )
                if not has_interactions:
                    raise ValueError("no interactions")
                ax = fp.plot_barcode()
                fig = ax.figure
                fig.set_size_inches(figsize[0], figsize[1])
                if display_title:
                    fig.suptitle(display_title, fontsize=11)
                fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
            except Exception:
                plt.close("all")  # clear any dirty matplotlib state before placeholder
                self._save_placeholder_png(
                    output_png,
                    title=display_title or f"ProLIF Interactions: {complex_pdb.stem}",
                    message="No interactions detected for barcode plot",
                    dpi=dpi,
                    figsize=figsize,
                )

            # Interaction network HTML
            # plot_lignetwork(kind="aggregate") raises KeyError:'ligand' on an
            # empty fingerprint (no interactions detected). Guard the call.
            if write_html:
                try:
                    has_interactions = (
                        not interaction_df.empty and len(interaction_df.columns) > 0
                    )
                    if has_interactions:
                        network_html = fp.plot_lignetwork(ligand_mol, kind="aggregate")
                        html_file = output_png.with_suffix(".html")
                        self._save_lignetwork_html(network_html, html_file)
                    else:
                        logger.info(
                            "ℹ️  ProLIF aggregate LigNetwork skipped for %s (no interactions detected)",
                            complex_pdb.name,
                        )
                except Exception as exc:
                    logger.warning(f"⚠️  ProLIF LigNetwork generation failed for {complex_pdb.name}: {exc}")

            # Flattened interactions table
            table_file = output_png.with_suffix(".csv")
            self._save_interaction_table(interaction_df, table_file)

            logger.info(f"✅ Created ProLIF map: {output_png.name}")
            return True

        except Exception as exc:
            logger.error(f"❌ Error creating ProLIF map for {complex_pdb.name}: {exc}", exc_info=True)
            return False
        finally:
            if working_pdb != complex_pdb:
                try:
                    working_pdb.unlink(missing_ok=True)
                except Exception:
                    pass

    def create_docking_interaction_map(
        self,
        complex_pdb: Path,
        poses_sdf: Path,
        output_png: Path,
        ligand_resname: str = "UNK",
        dpi: int = 300,
        figsize: Tuple[int, int] = (12, 8),
        write_html: bool = True,
        max_poses: Optional[int] = None,
        lignetwork_threshold: float = 0.3,
        count_occurrences: bool = False,
        display_title: Optional[str] = None,
    ) -> bool:
        """
        ProLIF docking workflow:
        - load all poses from SDF
        - run IFP on the pose iterable against one protein
        - export pose-wise barcode and aggregate/frame networks
        """
        if not PROLIF_AVAILABLE:
            logger.error("ProLIF stack not available")
            return False
        if not complex_pdb.exists():
            logger.error(f"Complex PDB file not found: {complex_pdb}")
            return False
        if not poses_sdf.exists():
            logger.error(f"Docking poses SDF not found: {poses_sdf}")
            return False

        output_png.parent.mkdir(parents=True, exist_ok=True)
        _patch_mda_vdwradii()
        working_pdb = complex_pdb

        try:
            working_pdb = _sanitize_pdb_elements(complex_pdb)
            universe = mda.Universe(str(working_pdb))
            protein_atoms = universe.select_atoms("protein")
            if len(protein_atoms) == 0:
                self._save_placeholder_png(
                    output_png,
                    title=display_title or f"ProLIF Docking Interactions: {complex_pdb.stem}",
                    message="Missing protein atoms for ProLIF docking workflow",
                    dpi=dpi,
                    figsize=figsize,
                )
                return True

            protein_mol = plf.Molecule.from_mda(protein_atoms, NoImplicit=False)
            pose_iterable = self._load_pose_iterable_from_sdf(poses_sdf, max_poses=max_poses)
            if not pose_iterable:
                self._save_placeholder_png(
                    output_png,
                    title=display_title or f"ProLIF Docking Interactions: {complex_pdb.stem}",
                    message="No valid docking poses found in SDF",
                    dpi=dpi,
                    figsize=figsize,
                )
                return True

            ligand_hint = self._normalize_resname(ligand_resname)
            if ligand_hint in {"UNK", "UNX", "LIG"}:
                inferred = self._infer_ligand_resname_from_name(poses_sdf.stem)
                if inferred:
                    ligand_hint = inferred

            fp = plf.Fingerprint(count=bool(count_occurrences))
            fp.run_from_iterable(pose_iterable, protein_mol)
            interaction_df = fp.to_dataframe(index_col="Pose")

            # Static barcode over docking poses.
            try:
                ax = fp.plot_barcode(xlabel="Pose")
                fig = ax.figure
                fig.set_size_inches(figsize[0], figsize[1])
                title_text = display_title or f"ProLIF Docking IFP ({complex_pdb.stem}) | ligand={ligand_hint}"
                fig.suptitle(f"{title_text} | poses={len(pose_iterable)}", fontsize=11)
                fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
            except Exception:
                self._save_placeholder_png(
                    output_png,
                    title=display_title or f"ProLIF Docking Interactions: {complex_pdb.stem}",
                    message="No interactions detected for docking-pose barcode plot",
                    dpi=dpi,
                    figsize=figsize,
                )

            # Aggregated + first-pose networks.
            # Guard: plot_lignetwork(kind="aggregate") raises KeyError:'ligand'
            # when there are no interactions across all poses.
            if write_html:
                try:
                    has_interactions = (
                        not interaction_df.empty and len(interaction_df.columns) > 0
                    )
                    if has_interactions:
                        aggregate_network = fp.plot_lignetwork(
                            pose_iterable[0],
                            kind="aggregate",
                            threshold=float(lignetwork_threshold),
                        )
                        aggregate_file = output_png.with_name(f"{output_png.stem}_aggregate.html")
                        self._save_lignetwork_html(aggregate_network, aggregate_file)
                    else:
                        logger.info(
                            "ℹ️  ProLIF aggregate LigNetwork skipped for %s (no interactions detected)",
                            complex_pdb.name,
                        )
                except Exception as exc:
                    logger.warning(f"⚠️  ProLIF aggregate LigNetwork failed for {complex_pdb.name}: {exc}")

                try:
                    # ProLIF 2.x _make_frame_df_from_fp crashes on an empty frame
                    # (KeyError: "None of ['ligand','protein','interaction','atoms']").
                    # Guard: only attempt the frame network when frame 0 has interactions.
                    frame_ifp = (fp.ifp[0] if (hasattr(fp, "ifp") and fp.ifp) else {}) or {}
                    if frame_ifp:
                        frame_network = fp.plot_lignetwork(
                            pose_iterable[0],
                            kind="frame",
                            frame=0,
                        )
                        frame_file = output_png.with_name(f"{output_png.stem}_pose0.html")
                        self._save_lignetwork_html(frame_network, frame_file)
                    else:
                        logger.info(
                            "ℹ️  ProLIF frame LigNetwork skipped for %s (no interactions in pose 0)",
                            complex_pdb.name,
                        )
                except Exception as exc:
                    message = str(exc)
                    if "None of ['ligand', 'protein', 'interaction', 'atoms']" in message:
                        logger.info(
                            "ℹ️  ProLIF frame LigNetwork skipped for %s (pose-level schema unavailable)",
                            complex_pdb.name,
                        )
                    else:
                        logger.warning(f"⚠️  ProLIF frame LigNetwork failed for {complex_pdb.name}: {exc}")

            table_file = output_png.with_suffix(".csv")
            self._save_interaction_table(interaction_df, table_file)
            freq_file = output_png.with_name(f"{output_png.stem}_frequency.csv")
            self._save_interaction_frequency(interaction_df, freq_file)

            logger.info(
                f"✅ Created ProLIF docking map: {output_png.name} "
                f"(poses={len(pose_iterable)}, ligand={ligand_hint})"
            )
            return True

        except Exception as exc:
            logger.error(f"❌ Error creating ProLIF docking map for {complex_pdb.name}: {exc}", exc_info=True)
            return False
        finally:
            if working_pdb != complex_pdb:
                try:
                    working_pdb.unlink(missing_ok=True)
                except Exception:
                    pass


def create_interaction_maps_for_all_complexes(
    complexes_dir: Path,
    output_dir: Path,
    ligand_resname: str = "UNK",
    dpi: int = 300,
    figsize: Tuple[int, int] = (14, 10),
    max_complexes: Optional[int] = None,
    overwrite: bool = False,
    write_html: bool = True,
    complex_manifest: Optional[List[Dict[str, Any]]] = None,
    complex_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
    poses_per_complex: Optional[int] = None,
    lignetwork_threshold: float = 0.3,
    count_occurrences: bool = False,
) -> Dict[str, Path]:
    """
    Create ProLIF outputs for all complex PDB files in a directory.
    """
    if not PROLIF_AVAILABLE:
        logger.warning("ProLIF stack not available - skipping interaction map generation")
        return {}

    mapper = ProLifInteractionMapper()
    output_dir.mkdir(parents=True, exist_ok=True)

    created_files: Dict[str, Path] = {}
    if complex_manifest:
        manifest_entries = list(complex_manifest)
        if max_complexes is not None:
            manifest_entries = manifest_entries[:max_complexes]

        for complex_info in manifest_entries:
            complex_name = str(complex_info.get("complex_name") or "").strip()
            if not complex_name:
                continue
            complex_file = complexes_dir / f"{complex_name}.pdb"
            pose_file = Path(str(complex_info.get("pose_file") or ""))
            if not complex_file.exists():
                logger.warning(f"⚠️  ProLIF skip (missing complex PDB): {complex_file}")
                continue
            if not pose_file.exists():
                logger.warning(f"⚠️  ProLIF skip (missing docking poses SDF): {pose_file}")
                continue

            output_stem = _metadata_output_stem(complex_name, complex_metadata)
            output_png = output_dir / f"{output_stem}.png"
            if output_png.exists() and not overwrite:
                created_files[complex_name] = output_png
                continue

            ligand_hint = str(complex_info.get("ligand_name") or ligand_resname)
            display_title = _metadata_display_title(complex_name, ligand_hint, complex_metadata)
            ok = mapper.create_docking_interaction_map(
                complex_file,
                pose_file,
                output_png,
                ligand_resname=ligand_hint,
                dpi=dpi,
                figsize=figsize,
                write_html=write_html,
                max_poses=poses_per_complex,
                lignetwork_threshold=lignetwork_threshold,
                count_occurrences=count_occurrences,
                display_title=display_title,
            )
            if ok:
                created_files[complex_name] = output_png
    else:
        complex_files = sorted(complexes_dir.glob("*.pdb"))
        if max_complexes is not None:
            complex_files = complex_files[:max_complexes]

        for complex_file in complex_files:
            output_png = output_dir / f"{complex_file.stem}_interaction_map.png"
            if output_png.exists() and not overwrite:
                created_files[complex_file.stem] = output_png
                continue
            ok = mapper.create_interaction_map(
                complex_file,
                output_png,
                ligand_resname=ligand_resname,
                dpi=dpi,
                figsize=figsize,
                write_html=write_html,
                display_title=f"ProLIF Interactions: {complex_file.stem}",
            )
            if ok:
                created_files[complex_file.stem] = output_png

    return created_files
