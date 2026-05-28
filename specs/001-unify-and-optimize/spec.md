# Feature Specification: Unified Visualization Overhaul

**Feature Branch**: `001-unify-and-optimize`  
**Created**: 2026-02-11  
**Status**: Draft  
**Input**: User description: "first we need to add a feature to visualziaiton if the file name is pdb code to reflect the protien name... adjust pandamap/prolif/ligplot parameters... correct RMSD figure... unify visualizations/raw data folders"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Protein Naming In Visuals (Priority: P1)

As a researcher, I need plots and interaction outputs to show human-readable protein names instead of raw PDB-code-like labels so I can interpret figures quickly.

**Why this priority**: It directly affects interpretability of every figure.

**Independent Test**: Run simplified pipeline on receptor files named with PDB codes and confirm generated plot axes/titles and manifests contain display names.

**Acceptance Scenarios**:

1. **Given** receptor identifiers with embedded PDB codes, **When** analysis runs, **Then** a protein display-name mapping is created and used in visual outputs.
2. **Given** no explicit name beyond PDB code, **When** mapping is generated, **Then** a stable fallback display name is produced.

---

### User Story 2 - Better Interaction Visual Quality (Priority: P1)

As a researcher, I need PandaMap, ProLIF, and LigPlot outputs to emphasize ligand shape/interactions with tuned defaults suitable for publication and batch runs.

**Why this priority**: Interaction figures are the main decision artifact.

**Independent Test**: Generate PandaMap/ProLIF/LigPlot outputs and verify parameterized behavior and non-empty output files.

**Acceptance Scenarios**:

1. **Given** PandaMap runs, **When** 3D output is generated, **Then** configured shape/surface-related flags are respected.
2. **Given** ProLIF/LigPlot runs, **When** defaults are used, **Then** output quality-oriented defaults and optimization behavior are applied.

---

### User Story 3 - Correct and Readable RMSD Figures (Priority: P1)

As a researcher, I need RMSD figures to be readable and correct for the analyzed set so clustering/diversity interpretation is reliable.

**Why this priority**: RMSD interpretation can mislead decisions if plots are unclear.

**Independent Test**: Run RMSD stage and verify heatmap/cluster/diversity plots render with proper labels and valid matrix handling.

**Acceptance Scenarios**:

1. **Given** a valid RMSD matrix, **When** figures are generated, **Then** heatmap labels/ticks are legible and consistent with poses.
2. **Given** sparse or degenerate data, **When** figures are generated, **Then** plotting code handles edge cases without crashing.

---

### User Story 4 - Unified Output Topology (Priority: P1)

As a researcher, I need all rendered visuals in one folder and all raw/tabular artifacts in another, with a linkage manifest between them.

**Why this priority**: Makes downstream report assembly and auditing much faster.

**Independent Test**: Run pipeline and verify `visualizations/` and `raw_data/` contain expected grouped outputs plus linkage metadata.

**Acceptance Scenarios**:

1. **Given** a completed run, **When** output organization executes, **Then** visualization files are consolidated under `visualizations/`.
2. **Given** a completed run, **When** output organization executes, **Then** raw CSV/JSON/TXT/log artifacts are consolidated under `raw_data/` and matched through a manifest.

### Edge Cases

- Missing or ambiguous receptor naming sources (pairlist absent or minimal receptor filenames).
- Missing optional tools (`PandaMap`, `py3Dmol`, `ProLIF`, `LigPlot`) should not hard-fail full pipeline.
- Existing output files must not be destructively overwritten without deterministic behavior.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST generate and persist a receptor/protein display-name mapping from detected receptors and available metadata.
- **FR-002**: System MUST apply display names in visualization labels/titles where protein identifiers are shown.
- **FR-003**: System MUST expose tuned PandaMap parameters for 2D/3D generation emphasizing ligand shape cues.
- **FR-004**: System MUST expose tuned ProLIF defaults and robust batch behavior for interaction-map generation.
- **FR-005**: System MUST expose tuned LigPlot defaults aimed at more complete interaction depiction and robust execution.
- **FR-006**: System MUST improve RMSD plot readability and robustness for valid and edge-case matrices.
- **FR-007**: System MUST consolidate rendered visuals into `visualizations/` and raw/tabular artifacts into `raw_data/`.
- **FR-008**: System MUST generate a manifest linking visualization files to source/raw analysis artifacts and relevant receptor/protein identifiers.

### Key Entities *(include if feature involves data)*

- **ProteinNameMapEntry**: Auto-detected mapping between receptor identifier, PDB code, and display name.
- **VisualizationManifestEntry**: Link row connecting visualization file path, visualization type, source artifact path(s), and protein/receptor metadata.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of generated figures that include protein labels use display names when a mapping is available.
- **SC-002**: PandaMap publication run produces non-empty outputs with configured 3D shape/surface settings for target complexes.
- **SC-003**: RMSD visualization stage completes without plotting exceptions on the target project dataset.
- **SC-004**: After run completion, `visualizations/` and `raw_data/` both exist and include a manifest that links generated visuals to raw artifacts.
