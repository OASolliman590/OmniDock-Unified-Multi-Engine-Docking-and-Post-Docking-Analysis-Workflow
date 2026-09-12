#!/usr/bin/env python3
"""
Setup script for PDB Prepare Wizard
"""

from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from pathlib import Path


class BuildWithPreparationScript(build_py):
    """Keep the shell backend beside the installed preparation Python module."""
    def run(self):
        super().run()
        self.copy_file("prep_autodock_enhanced.sh", str(Path(self.build_lib) / "prep_autodock_enhanced.sh"))

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="pdb-prepare-wizard",
    version="3.0.1",
    author="Molecular Docking Pipeline",
    author_email="",
    description="A comprehensive tool for preparing PDB files for molecular docking studies with advanced PLIP integration and comprehensive interaction analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(exclude=("test", "test.*", "audit", "audit.*")),
    py_modules=["main", "core_pipeline", "cli_pipeline", "interactive_pipeline",
                "batch_pdb_preparation", "autodock_preparation"],
    package_data={"post_docking_analysis": ["config/*.yaml"]},
    cmdclass={"build_py": BuildWithPreparationScript},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    extras_require={
        "chemistry": ["rdkit>=2024.3", "meeko>=0.6", "gemmi>=0.6"],
        "interactions": ["rdkit>=2024.3", "plip>=2.2", "prolif>=2.0", "MDAnalysis>=2.6"],
        "dev": ["pytest>=8", "build>=1", "wheel"],
        "notebooks": ["jupyter>=1.0"],
    },
    entry_points={
        "console_scripts": [
            "pdb-prepare-wizard=main:main",
            "pdb-wizard-workflow=workflow.cli:main",
            "pdb-wizard-interactive=workflow.interactive:run_interactive_workflow",
            "pdb-wizard-legacy-interactive=interactive_pipeline:run_interactive_pipeline",
            "pdb-wizard-cli=cli_pipeline:main",
            "pdb-wizard-batch=batch_pdb_preparation:main",
            "pdb-wizard-prepare-docking=docking.cli:prepare_docking_main",
            "pdb-wizard-dock=docking.cli:dock_main",
            "post-docking-analysis=post_docking_analysis.cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
) 
