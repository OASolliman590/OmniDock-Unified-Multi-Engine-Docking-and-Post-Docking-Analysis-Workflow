#!/usr/bin/env python3
"""
Setup script for PDB Prepare Wizard
"""

from setuptools import setup, find_packages

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
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
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
