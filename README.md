# Ultrafast Laser TTM Toolbox

[![DOI](https://zenodo.org/badge/1249804402.svg)](https://doi.org/10.5281/zenodo.20389305)

This folder is a clean staging copy for a public repository built from the reusable parts of the original workspace. The goal is to preserve the core modeling functionality while removing project-specific clutter such as manuscript validation campaigns, archived backups, and generated result folders.

The codebase centers on MATLAB solvers for pulsed-laser heating in metals using two-temperature-model workflows. The staged version is still tungsten-oriented in its defaults, but the structure is intended to be reusable for other materials, pulse widths, spot sizes, repetition rates, and scanning conditions.

## Repository Purpose

This staging copy is meant to become a portable, public-facing repository candidate with these priorities:

- keep the general solver stack
- keep a small number of reusable batch entry points
- keep enough documentation for another user to understand what each script does
- make the repository legible to coding agents as well as human users
- avoid embedding one specific manuscript as the repository's organizing principle
- avoid shipping generated output data as source content

## Agent-Friendly Design

This repository is intentionally being shaped so coding agents can work in it reliably rather than treating agent support as an accident.

Agent-friendly features already present in the staged copy:

- relative-path output handling instead of machine-local hardcoded paths
- solver code concentrated in `src/` rather than scattered across output or archive folders
- lightweight editable examples in `examples/` that call into reusable solver functions
- reusable batch workflows in `scripts/batch/` with shared helper utilities
- repository-scoped guidance in `AGENTS.md`
- a dedicated agent quickstart in `docs/AGENT_QUICKSTART.md`

Design choice for the 0D solver:

- `Surface_Point_Solver.m` remains in `src/` because it now behaves like a supported solver entry point with a config-struct interface and returned results, not just as a teaching-only script

## What Is Included

### Core solvers in `src/`

- `Depth_Profile_Solver.m` — 1D depth-resolved multi-pulse TTM solver for coupled electron and lattice temperature evolution through depth
- `Radial_Profile_Solver.m` — radial surface-profile solver for steady-state or pulse-accumulation style radial temperature studies
- `Inversion_Quantifier.m` — post-processing and analysis of electron-lattice inversion behavior
- `Single_Pulse_Visualizer.m` — single-pulse visualization workflow for inspecting early-time temperature evolution
- `Surface_Point_Solver.m` — reduced 0D surface-point TTM model for simpler studies and quick parameter sweeps
- `Scanning_Beam_Solver.m` — moving-beam surface scan model with scanning kinematics and diffusion

### Example entry points in `examples/`

- `README.md` — guide to the example scripts and recommended starting order
- `Example_Surface_Point_Baseline.m` — simplest baseline case using the 0D surface-point solver
- `Example_Depth_Profile_Baseline.m` — baseline case for the main 1D depth solver
- `Example_Radial_Profile_Baseline.m` — baseline case for the radial-profile solver
- `Example_Scanning_Beam_Baseline.m` — canonical moving-beam baseline example
- `Scanning_Beam_Single_Run.m` — legacy compatibility alias for the scanning-beam example

### Scripts in `scripts/`

- `README.md` — guide to validation and batch-runner entry points
- `Verify_Public_Repo_Smoke.m` — lightweight runtime smoke test for the public-facing solver surfaces

### Batch runners in `scripts/batch/`

- `Batch_Depth_Profile.m` — generic depth-only parameter sweep entry point
- `Batch_Multi_Solver_Study.m` — combined depth, radial, single-pulse, and inversion workflow for multi-case studies
- `Batch_Inversion_Analysis.m` — inversion-analysis batch runner across multiple cases
- `Batch_Depth_and_Radial_Profile.m` — legacy compatibility alias retained for older workflow references

### Supporting documentation in `docs/`

- `AGENT_QUICKSTART.md` — fast map of the repo for coding agents and automation tools
- `RESULT_CONTRACT_SCHEMA.json` — machine-readable note describing the shared solver result contract
- `REPOSITORY_NOTES.md` — repository curation notes and publishing guidance
- `PROJECT_CONTEXT.md` — context for what this repository represents and what has been intentionally removed
- `VALIDATION_GUIDE.md` — lightweight manual validation workflow for the canonical examples and solver outputs

## What Was Left Out

The following content was intentionally excluded from the staged public-repo copy:

- manuscript-specific batch orchestration
- validation plans tied to one experimental study
- archived backups and duplicate "Finished Files"
- generated output folders and saved result snapshots
- one-off experiment directories whose main value is historical context rather than reusable functionality

This is a scope decision, not a claim that the removed content is unimportant. It simply does not belong in the first pass of a general-purpose public source repository.

## Folder Layout

```text
Public_Repo_Staging/
  src/              % reusable solver functions
  examples/         % editable single-run examples
  scripts/
    README.md       % script-level entry point guide
    Verify_Public_Repo_Smoke.m
    batch/          % reusable batch entry points
  docs/             % repository context and publishing notes
  outputs/          % default destination for generated results
```

All staged scripts now write default results into `outputs/` at the repository root rather than into the original local workspace layout.

## How To Start

### MATLAB setup

1. Open MATLAB with this folder as the working directory.
2. Add the solver folder to the path:

```matlab
addpath(genpath('src'))
```

3. Run one of the main entry points below.

### Suggested first runs

- `examples/Example_Surface_Point_Baseline` for the simplest starting case
- `examples/Example_Depth_Profile_Baseline` for the main 1D depth model
- `examples/Example_Radial_Profile_Baseline` when the main question is radial spread
- `Single_Pulse_Visualizer()` for early-time pulse behavior and figures
- `examples/Example_Scanning_Beam_Baseline` for a simple moving-beam example
- `Batch_Depth_Profile` for a small multi-case batch run

### Parameter editing model

The repository currently uses a mixed style:

- some files are function-based and accept a config struct
- some files are script-style entry points with editable input sections near the top

That means the fastest way to explore is often to edit the user input block in a script, while more structured automation is usually easiest through the function-based solvers.

## Solver Interface Matrix

| Solver | File | Primary interface | Returns results struct | Common contract fields | Best use |
| --- | --- | --- | --- | --- | --- |
| Surface point | `src/Surface_Point_Solver.m` | `cfg` struct | yes | `solver`, `solverId`, `contractVersion`, `material`, `outputFile`, `outputDir`, `inputConfig` | fastest 0D pulse-accumulation studies |
| Depth profile | `src/Depth_Profile_Solver.m` | `cfg` struct | yes | `solver`, `solverId`, `contractVersion`, `material`, `nPulses`, `wallTime_s`, `outputFile`, `outputDir`, `inputConfig` | main 1D depth-resolved multi-pulse workflow |
| Radial profile | `src/Radial_Profile_Solver.m` | `cfg` struct | yes | `solver`, `solverId`, `contractVersion`, `material`, `nPulses`, `wallTime_s`, `outputFile`, `outputDir`, `inputConfig` | radial spread and footprint studies |
| Single pulse | `src/Single_Pulse_Visualizer.m` | `cfg` struct | yes | `solver`, `solverId`, `contractVersion`, `material`, `wallTime_s`, `outputFile`, `outputDir`, `inputConfig` | early-time single-pulse inspection |
| Inversion analysis | `src/Inversion_Quantifier.m` | `cfg` struct | yes | `solver`, `solverId`, `contractVersion`, `material`, `nPulses`, `wallTime_s`, `outputFile`, `outputDir`, `depthOutputFile` | per-pulse inversion statistics |
| Scanning beam | `src/Scanning_Beam_Solver.m` | `params` struct plus optional `outputDir`, `savePlots` | yes | `solver`, `solverId`, `contractVersion`, `material`, `nPulses`, `wallTime_s`, `outputFile`, `outputDir`, `inputConfig` | moving-beam surface scan modeling |

## Common Result Contract

The solvers do not return identical full payloads, but the staged public repo now exposes a shared minimal result contract to make scripting and agent-driven tooling more reliable.

Fields you can expect on the main solver entry points:

- `solver` — human-readable solver label
- `solverId` — stable machine-friendly identifier
- `contractVersion` — current shared result-contract version
- `material` — active material preset or mode
- `outputFile` — primary text output path when one is written
- `outputDir` — output folder path
- `inputConfig` — input struct used to invoke the solver

Additional common fields appear where they are physically meaningful, especially `nPulses` and `wallTime_s`.

## Recommended Workflow

For a new user, the most sensible progression is:

1. Run `examples/Example_Surface_Point_Baseline` to confirm the environment and inspect basic pulse behavior.
2. Move to `examples/Example_Depth_Profile_Baseline` for depth-resolved multi-pulse accumulation.
3. Use `examples/Example_Radial_Profile_Baseline` if the main question is radial spread or footprint.
4. Use `Scanning_Beam_Solver` or `examples/Example_Scanning_Beam_Baseline` for translating stationary heating logic into a moving-laser process.
5. Use the batch runners only after the single-case behavior is understood.

## Design Context

The staged repository is broader than the manuscript that motivated some of the original work. The intent here is not to freeze the code around one experiment, one material, or one laser system. The current defaults reflect the original development history, but the solver concepts are applicable to broader pulsed-laser thermal studies provided the underlying assumptions remain appropriate.

Areas that are most likely to need user adaptation include:

- material properties and presets
- optical absorption assumptions
- pulse profile and repetition rate ranges
- diffusion-domain sizing
- output naming and metadata conventions

## Outputs And Reproducibility

Generated files are directed to `outputs/` by default. That folder is kept in the staged layout because the scripts expect a writable destination, but generated results themselves are excluded from version control.

If this staging copy becomes the actual public repository, a later cleanup pass should likely add:

- a small set of canonical example cases
- a concise result interpretation guide
- stronger input validation for public users

## License

This repository is released under the MIT License. See `LICENSE`.

## Citation

If this repository contributes to published work, please cite the software record in `CITATION.cff`.

- Concept DOI (latest release): `10.5281/zenodo.20389305`
- Version DOI (`v0.1.0`): `10.5281/zenodo.20389306`

## Acknowledgments and Funding

This work was supported by the National Science Foundation under Award No. CMMI-2412544, Collaborative Research: Additive Manufacturing of Crack-Free Tungsten Using Ultrashort Pulsed Lasers (PI: Dr. Anming Hu, Division of Civil, Mechanical, and Manufacturing Innovation, NSF Program: AM-Advanced Manufacturing).

The authors gratefully acknowledge Drs. Yanfei Gao, Wenda Tan, and Seungha Shin for their contributions and collaboration on this project.

Additional support was provided by the University of Tennessee, Knoxville, through a hiring package. D.F. gratefully acknowledges support from the UTK 100 Talented PhD Scholarship.

Support for the Center for Materials Processing from the State of Tennessee and the Tennessee Higher Education Commission is also gratefully acknowledged.

## Contributing

For small fixes, documentation improvements, and workflow polish, see `CONTRIBUTING.md`.

## Public-Release Caveats

- The defaults are tungsten-focused and reflect the development history of the project.
- Naming is still partly shaped by the original research workflow rather than by a polished public API.
- Some files under `examples/` are intentionally user-edited entry points rather than polished software interfaces.
- MATLAB editor diagnostics currently report a few existing unused-variable and formatting warnings in the solver files; those were not part of this repository-organization pass.

## Related Documentation

For more context on how this staging repository was curated, see `docs/REPOSITORY_NOTES.md` and `docs/PROJECT_CONTEXT.md`. For example selection and starter workflows, see `examples/README.md`. For validation workflow, see `docs/VALIDATION_GUIDE.md`. For agent-specific editing and maintenance guidance, see `AGENTS.md`, `docs/AGENT_QUICKSTART.md`, and `docs/RESULT_CONTRACT_SCHEMA.json`. For citation and reuse metadata, see `CITATION.cff`, `LICENSE`, and `CONTRIBUTING.md`.