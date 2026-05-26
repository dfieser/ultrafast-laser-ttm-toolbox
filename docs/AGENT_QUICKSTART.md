# Agent Quickstart

This document is a fast repo map for coding agents and automation tools.

## Primary Intent

Treat agent usability as a supported repository feature. The goal is not only to make the MATLAB code runnable, but also to make the project understandable and maintainable by an automated assistant with minimal repo archaeology.

## Start Here By Task Type

### User wants to run something quickly

- `examples/README.md`
- `examples/Example_Surface_Point_Baseline.m`
- `examples/Example_Depth_Profile_Baseline.m`
- `examples/Example_Radial_Profile_Baseline.m`
- `examples/Example_Scanning_Beam_Baseline.m`

### User wants solver behavior changed

- `src/Surface_Point_Solver.m`
- `src/Depth_Profile_Solver.m`
- `src/Radial_Profile_Solver.m`
- `src/Single_Pulse_Visualizer.m`
- `src/Inversion_Quantifier.m`
- `src/Scanning_Beam_Solver.m`

### User wants batch-study behavior changed

- `scripts/batch/Batch_Depth_Profile.m`
- `scripts/batch/Batch_Multi_Solver_Study.m`
- `scripts/batch/Batch_Inversion_Analysis.m`
- `src/RepoUtils.m`

### User wants repository positioning or onboarding changed

- `README.md`
- `AGENTS.md`
- `docs/PROJECT_CONTEXT.md`
- `docs/REPOSITORY_NOTES.md`

## Structural Rules

- Keep reusable solver logic in `src/`.
- Keep editable walkthroughs and starter runs in `examples/`.
- Keep repeated multi-case workflows in `scripts/batch/`.
- Keep generated artifacts in `outputs/`.
- Keep agent-facing maintenance rules in `AGENTS.md` and this file.

## Important Design Decisions

- `Surface_Point_Solver.m` stays in `src/` because it is now a config-driven solver API that returns results, not just a tutorial script.
- `Example_Scanning_Beam_Baseline.m` is the canonical moving-beam example in `examples/` and delegates the reusable logic to `Scanning_Beam_Solver.m`.
- `Scanning_Beam_Single_Run.m` is retained as a compatibility alias for older references.
- Batch scripts should reuse `RepoUtils.m` rather than carrying copy-pasted local helpers.

## Shared Result Contract

The main solver entry points now expose a shared minimal result contract for automation.

Expected cross-solver fields:

- `solver`
- `solverId`
- `contractVersion`
- `material`
- `outputFile`
- `outputDir`
- `inputConfig`

Additional fields such as `nPulses` and `wallTime_s` are standardized where they make sense.

For automation that benefits from a machine-readable reference, see `docs/RESULT_CONTRACT_SCHEMA.json`.

## What Good Agent Changes Look Like

- small, local edits tied to one solver or one documentation surface
- preserving relative-path portability
- keeping examples thin
- reducing duplication instead of adding near-copies
- updating docs when file roles or workflow entry points change

## What Bad Agent Changes Look Like

- copying solver internals into examples or batch scripts
- reintroducing manuscript-specific folders or output dumps
- mixing generated artifacts with source files
- adding machine-specific paths
- moving supported solver APIs into `examples/` without a strong reason

## Validation Guidance

- validate touched MATLAB files with editor diagnostics when possible
- prefer file-scoped validation after each meaningful edit
- for docs, verify that names and folder references still match the actual repo layout

## If Unsure

Choose the smallest change that preserves:

- solver portability
- clear separation between `src/`, `examples/`, and `scripts/batch/`
- the repo's public-facing and agent-friendly structure