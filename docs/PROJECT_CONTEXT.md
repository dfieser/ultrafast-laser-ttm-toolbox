# Project Context

This repository was curated from a larger research workspace for pulsed-laser thermal modeling in MATLAB. The original workspace mixed several kinds of content:

- core solver code
- batch automation for broad parameter sweeps
- manuscript-specific validation workflows
- archived backups and duplicate files
- generated outputs from many runs

For a public repository, that mixture is not ideal. This repository separates the reusable modeling core from the research-history residue.

## Conceptual Scope

The repository focuses on two-temperature-model simulations for pulsed-laser heating of metals, at several levels of fidelity:

- 0D surface-point modeling
- 1D depth-resolved modeling
- radial-profile analysis
- single-pulse visualization
- scanning-beam surface heating
- multi-case batch execution

The defaults and comments still reflect the project's tungsten-focused development history. The modeling structure, though, is broader than one material system or one laser setup.

## What This Repository Is

This repository is best understood as a reusable research codebase for exploring pulsed-laser thermal behavior under configurable conditions.

The repository is also organized so coding agents can maintain and extend it as a first-class environment rather than an opaque research dump.

The repository aims to support questions like:

- how electron and lattice temperatures evolve during and after a pulse
- how residual heating accumulates over many pulses
- how temperature varies with depth or radius
- how thermal behavior changes with repetition rate, average power, spot size, and pulse width
- how a moving beam changes the thermal footprint relative to a stationary beam

## What This Repository Is Not

In its current form, this repository is not meant to be:

- a paper-specific artifact bundle
- a complete archival record of every simulation ever run
- a packaged application with a GUI
- a fully generalized materials database
- a ready-to-use industrial process simulator

## Why Manuscript-Specific Content Was Removed

Some of the original work was organized around one manuscript and its validation structure. That material can still be valuable internally, but in a public repository it:

- makes the repository look narrower than it really is
- introduces folders that are meaningful only with manuscript context
- buries the reusable solver entry points under study-specific orchestration
- encourages users to treat one validation campaign as the only intended use case

The public repository therefore keeps the general solver capabilities and leaves the manuscript-specific automation behind.

## Solver Relationship Overview

Think of the files in layers:

1. simple pulse-level or point-level inspection
2. depth-resolved or radial thermal modeling
3. moving-beam process modeling
4. example entry scripts for editable runs
5. batch orchestration for repeated studies

The new example layer exists to give public users a stable starting surface without forcing them to edit the larger solver files immediately.

That hierarchy matters because most users should understand the lower-level behavior before relying on the higher-level batch workflows.

## Assumptions To Keep In Mind

The code still reflects modeling assumptions that may not fit a future use case, even as the repository is generalized. A public user should review at least:

- material property definitions
- absorption assumptions
- boundary conditions
- spatial and temporal resolution choices
- solver tolerances and stopping criteria
- how to interpret outputs physically

## Recommended Direction For Future Cleanup

As the repository matures, the best long-term direction is to shift from research-script organization toward clearer public structure:

- common configuration conventions across solvers
- fewer duplicated helper patterns
- consistent naming and units in outputs
- documented example cases with expected behavior
- a cleaner distinction between reusable solver functions in `src/` and editable walkthroughs in `examples/`
- explicit agent-facing guidance and stable edit surfaces for future automated maintenance