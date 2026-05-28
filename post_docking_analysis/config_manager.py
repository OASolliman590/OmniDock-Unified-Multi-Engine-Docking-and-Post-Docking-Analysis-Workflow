"""
Configuration management for post-docking analysis pipeline.
Supports both JSON and YAML configuration files.
"""
from pathlib import Path
import json
from typing import Dict, Any
import yaml

# Default configuration
DEFAULT_CONFIG = {
    # Analysis Parameters
    "analysis": {
        "docking_types": ["vina", "gnina"],
        "comparative_benchmark": "*",
        "binding_affinity_analysis": True,
        "rmsd_analysis": True,
        "generate_visualizations": True,
        "extract_poses": True
    },
    
    # Input/Output Directories
    "paths": {
        "input_dir": "",
        "output_dir": "./post_docking_results",
        "receptors_dir": "",
        "gnina_out_dir": ""
    },
    
    # Pose Extraction Parameters
    "pose_extraction": {
        "extract_all_poses": False,
        "best_pose_criteria": "affinity",
        "output_formats": ["pdb"]
    },
    
    # Binding Affinity Analysis Parameters
    "binding_affinity": {
        "strong_binder_threshold": -8.0,
        "top_performers_count": 10,
        "analyze_by_protein": True,
        "analyze_by_ligand": True
    },
    
    # RMSD Analysis Parameters
    "rmsd": {
        "clustering_method": "kmeans",
        "kmeans_clusters": 3,
        "dbscan_epsilon": 2.0,
        "dbscan_min_samples": 2
    },
    
    # Visualization Parameters
    "visualization": {
        "output_formats": ["png"],
        "dpi": 300,
        "generate_3d": True,
        "generate_2d_interactions": True,
        "pandamap": {
            "show_surface": True,
            "show_3d_cues": True,
            "width": 1200,
            "height": 900,
            "generate_text_reports": True,
            "estimate_delta_g": True,
            "write_capability_manifest": True,
            "qc_min_non_generic_ligand_ratio": 0.60,
            "qc_require_reports_when_supported": True
        },
        "prolif": {
            "dpi": 350,
            "figsize": [14, 10],
            "max_complexes": 60
        },
        "ligplot": {
            "contact_type": "2",
            "hydrogenate": "auto",
            "strip_metals": "auto",
            "no_abort": True
        },
        "plip": {
            "enabled": True,
            "output_formats": ["png", "xml", "txt"],
            "generate_pymol_session": True,
            "timeout_seconds": 120,
            "interaction_types": [
                "hydrophobic", "hbond", "waterbridge",
                "saltbridge", "pistacking", "pication",
                "halogen", "metal"
            ]
        }
    },
    
    # Advanced Options
    "advanced": {
        "fix_chains": False,
        "run_additional_docking": False,
        "directory_structure": "AUTO",
        "receptor_pattern": "*receptor*.pdb",
        "ligand_pattern": "*ligand*.sdf",
        "docking_result_pattern": "*out*.pdbqt"
    },
    "quality_control": {
        "ligand_qc": {
            "enabled": True,
            "admet_filters_enabled": True,
            "max_lipinski_violations": 1,
            "max_molecular_weight": 650.0,
            "max_logp": 6.0,
            "max_tpsa": 180.0,
            "max_rotatable_bonds": 15,
            "max_formal_charge_abs": 2,
            "min_heavy_atom_count": 6,
            "block_pains": True,
            "block_brenk": True,
            "block_reactive": True,
        },
        "receptor_qc": {
            "enabled": True,
            "min_atom_count": 100,
            "min_heavy_atom_count": 60,
            "min_chain_count": 1,
            "max_coordinate_span": 500.0,
        },
        "hit_classification": {
            "policy": "target_aware",
            "strong_percentile": 0.10,
            "moderate_percentile": 0.35,
        },
    }
}

class ConfigManager:
    """
    Configuration manager for the post-docking analysis pipeline.
    Supports both JSON and YAML configuration files.
    """
    
    def __init__(self, config_file: str = None):
        """
        Initialize the configuration manager.
        
        Parameters
        ----------
        config_file : str, optional
            Path to configuration file
        """
        self.config = self._deep_copy_dict(DEFAULT_CONFIG)
        self._normalize_quality_control()
        
        if config_file:
            self.load_config(config_file)
    
    def _deep_copy_dict(self, d):
        """Create a deep copy of a dictionary."""
        if isinstance(d, dict):
            return {k: self._deep_copy_dict(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [self._deep_copy_dict(item) for item in d]
        else:
            return d
    
    def load_config(self, config_file: str):
        """
        Load configuration from a file (JSON or YAML).
        
        Parameters
        ----------
        config_file : str
            Path to configuration file
        """
        try:
            file_config = read_config_file(config_file)
            
            # Update default config with file config
            self._update_nested_dict(self.config, file_config)
            self._normalize_quality_control()
            print(f"✅ Configuration loaded from: {config_file}")
        except Exception as e:
            print(f"❌ Error loading configuration file: {e}")
    
    def _update_nested_dict(self, base_dict, update_dict):
        """Update nested dictionary values."""
        if not isinstance(update_dict, dict):
            raise TypeError(
                f"Expected nested configuration mapping, got {type(update_dict).__name__}"
            )
        for key, value in update_dict.items():
            if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
                self._update_nested_dict(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def save_config(self, config_file: str):
        """
        Save current configuration to a file.
        
        Parameters
        ----------
        config_file : str
            Path to configuration file
        """
        config_path = Path(config_file)
        try:
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(config_path, 'w') as f:
                suffix = config_path.suffix.lower()
                if suffix in {'.yaml', '.yml'}:
                    yaml.dump(self.config, f, indent=2, default_flow_style=False)
                else:
                    json.dump(self.config, f, indent=2)
            print(f"✅ Configuration saved to: {config_file}")
        except Exception as e:
            print(f"❌ Error saving configuration file: {e}")
    
    def get(self, key_path: str, default=None):
        """
        Get a configuration value using dot notation (e.g., 'analysis.docking_types').
        
        Parameters
        ----------
        key_path : str
            Configuration key path (dot-separated)
        default : any, optional
            Default value if key not found
            
        Returns
        -------
        any
            Configuration value
        """
        keys = key_path.split('.')
        current = self.config
        
        try:
            for key in keys:
                current = current[key]
            return current
        except (KeyError, TypeError):
            return default
    
    def set(self, key_path: str, value):
        """
        Set a configuration value using dot notation.
        
        Parameters
        ----------
        key_path : str
            Configuration key path (dot-separated)
        value : any
            Configuration value
        """
        keys = key_path.split('.')
        current = self.config
        
        # Navigate to the parent dictionary
        for key in keys[:-1]:
            if key not in current or not isinstance(current[key], dict):
                current[key] = {}
            current = current[key]
        
        # Set the final value
        current[keys[-1]] = value
    
    def update(self, config_dict: Dict[str, Any]):
        """
        Update multiple configuration values.
        
        Parameters
        ----------
        config_dict : Dict[str, Any]
            Dictionary of configuration values to update
        """
        if not isinstance(config_dict, dict):
            raise TypeError(
                f"ConfigManager.update expects a dictionary, got {type(config_dict).__name__}"
            )
        self._update_nested_dict(self.config, config_dict)
        self._normalize_quality_control()

    def _normalize_quality_control(self):
        """Normalize QC and classification controls to safe numeric/boolean types."""
        qc = self.config.setdefault("quality_control", {})
        if not isinstance(qc, dict):
            qc = {}
            self.config["quality_control"] = qc

        ligand_qc = qc.setdefault("ligand_qc", {})
        receptor_qc = qc.setdefault("receptor_qc", {})
        hit_cls = qc.setdefault("hit_classification", {})
        if not isinstance(ligand_qc, dict):
            ligand_qc = {}
            qc["ligand_qc"] = ligand_qc
        if not isinstance(receptor_qc, dict):
            receptor_qc = {}
            qc["receptor_qc"] = receptor_qc
        if not isinstance(hit_cls, dict):
            hit_cls = {}
            qc["hit_classification"] = hit_cls

        def _to_bool(value, default):
            return bool(default if value is None else value)

        def _to_int(value, default, minimum=0):
            try:
                return max(minimum, int(value))
            except Exception:
                return int(default)

        def _to_float(value, default, minimum=None, maximum=None):
            try:
                number = float(value)
            except Exception:
                number = float(default)
            if minimum is not None:
                number = max(float(minimum), number)
            if maximum is not None:
                number = min(float(maximum), number)
            return number

        ligand_qc["enabled"] = _to_bool(ligand_qc.get("enabled"), True)
        ligand_qc["admet_filters_enabled"] = _to_bool(ligand_qc.get("admet_filters_enabled"), True)
        ligand_qc["max_lipinski_violations"] = _to_int(ligand_qc.get("max_lipinski_violations"), 1, minimum=0)
        ligand_qc["max_molecular_weight"] = _to_float(ligand_qc.get("max_molecular_weight"), 650.0, minimum=1.0)
        ligand_qc["max_logp"] = _to_float(ligand_qc.get("max_logp"), 6.0)
        ligand_qc["max_tpsa"] = _to_float(ligand_qc.get("max_tpsa"), 180.0, minimum=0.0)
        ligand_qc["max_rotatable_bonds"] = _to_int(ligand_qc.get("max_rotatable_bonds"), 15, minimum=0)
        ligand_qc["max_formal_charge_abs"] = _to_int(ligand_qc.get("max_formal_charge_abs"), 2, minimum=0)
        ligand_qc["min_heavy_atom_count"] = _to_int(ligand_qc.get("min_heavy_atom_count"), 6, minimum=0)
        ligand_qc["block_pains"] = _to_bool(ligand_qc.get("block_pains"), True)
        ligand_qc["block_brenk"] = _to_bool(ligand_qc.get("block_brenk"), True)
        ligand_qc["block_reactive"] = _to_bool(ligand_qc.get("block_reactive"), True)

        receptor_qc["enabled"] = _to_bool(receptor_qc.get("enabled"), True)
        receptor_qc["min_atom_count"] = _to_int(receptor_qc.get("min_atom_count"), 100, minimum=1)
        receptor_qc["min_heavy_atom_count"] = _to_int(receptor_qc.get("min_heavy_atom_count"), 60, minimum=1)
        receptor_qc["min_chain_count"] = _to_int(receptor_qc.get("min_chain_count"), 1, minimum=1)
        receptor_qc["max_coordinate_span"] = _to_float(
            receptor_qc.get("max_coordinate_span"),
            500.0,
            minimum=1.0,
        )

        hit_cls["policy"] = str(hit_cls.get("policy") or "target_aware").strip().lower() or "target_aware"
        hit_cls["strong_percentile"] = _to_float(hit_cls.get("strong_percentile"), 0.10, minimum=0.0, maximum=1.0)
        hit_cls["moderate_percentile"] = _to_float(
            hit_cls.get("moderate_percentile"),
            0.35,
            minimum=0.0,
            maximum=1.0,
        )

def load_config(config_file: str = None):
    """
    Load configuration from a file or return default configuration.
    
    Parameters
    ----------
    config_file : str, optional
        Path to configuration file
        
    Returns
    -------
    ConfigManager
        Configuration manager instance
    """
    return ConfigManager(config_file)


def read_config_file(config_file: str) -> Dict[str, Any]:
    """
    Read a JSON/YAML config file as a raw mapping.

    Raises FileNotFoundError / ValueError when unreadable or malformed.
    """
    config_path = Path(config_file)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_file}")

    with open(config_path, "r", encoding="utf-8") as handle:
        suffix = config_path.suffix.lower()
        if suffix in {".yaml", ".yml"}:
            parsed = yaml.safe_load(handle)
        else:
            parsed = json.load(handle)

    if parsed is None:
        return {}
    if not isinstance(parsed, dict):
        raise ValueError(
            f"Configuration root must be a mapping/object, got {type(parsed).__name__}"
        )
    return parsed


def _select_overlay_from_merged(raw: Dict[str, Any], merged: Dict[str, Any]) -> Dict[str, Any]:
    """
    Keep only keys explicitly provided in `raw`, using normalized values from `merged`.
    """
    overlay: Dict[str, Any] = {}
    for key, raw_value in raw.items():
        if isinstance(raw_value, dict) and isinstance(merged.get(key), dict):
            overlay[key] = _select_overlay_from_merged(raw_value, merged[key])
        else:
            overlay[key] = merged.get(key)
    return overlay


def load_config_overrides(config_file: str) -> Dict[str, Any]:
    """
    Load user config and return normalized explicit overrides only.

    This preserves default behavior in call sites unless the user provided the key.
    """
    raw_config = read_config_file(config_file)
    manager = ConfigManager()
    manager.update(raw_config)
    return _select_overlay_from_merged(raw_config, manager.config)
