# Project Context

This repository staging copy comes from a larger research workspace built around pulsed-laser thermal modeling in MATLAB. The original workspace mixed together several different kinds of content:

- core solver code
- batch automation for broad parameter sweeps
- manuscript-specific validation workflows
- archived backups and duplicate files
- generated outputs from many runs

For a public repository, that mixture is not ideal. This staged copy exists to separate the reusable modeling core from the research-history residue.

## Conceptual Scope

The repository focuses on two-temperature-model style simulations for pulsed-laser heating of metals, including several levels of fidelity:

- 0D surface-point modeling
- 1D depth-resolved modeling
- radial-profile analysis
- single-pulse visualization
- scanning-beam surface heating
- multi-case batch execution

The defaults and comments still reflect the project's tungsten-focused development history, but the modeling structure is broader than one material system or one laser setup.

## What This Repository Is

This repository is best understood as a reusable research codebase for exploring pulsed-laser thermal behavior under configurable conditions.

It is also being intentionally organized so coding agents can work in it as a first-class maintenance and extension environment rather than as an opaque research dump.

It is intended to support questions like:

- how electron and lattice temperatures evolve during and after a pulse
- how residual heating accumulates over many pulses
- how temperature varies with depth or radius
- how thermal behavior changes with repetition rate, average power, spot size, and pulse width
- how a moving beam changes the thermal footprint relative to a stationary beam

## What This Repository Is Not

This repository is not intended, in its current form, to be:

- a paper-specific artifact bundle
- a complete archival record of every simulation ever run
- a packaged application with a GUI
- a fully generalized materials database
- a turnkey industrial process simulator

## Why Manuscript-Specific Content Was Removed

Some of the original work was organized around one manuscript and its validation structure. That material can still be valuable internally, but it creates several problems in a public repository:

- it makes the repository look narrower than it really is
- it introduces folders that are meaningful only with manuscript context
- it buries the reusable solver entry points under study-specific orchestration
- it encourages users to treat one validation campaign as the only intended use case

The staging copy therefore keeps the general solver capabilities and leaves the manuscript-specific automation behind.

## Solver Relationship Overview

The staged files can be thought of in layers:

1. simple pulse-level or point-level inspection
2. depth-resolved or radial thermal modeling
3. moving-beam process modeling
4. example entry scripts for editable runs
5. batch orchestration for repeated studies

The new example layer exists to give public users a stable starting surface without forcing them to edit the larger solver files immediately.

That hierarchy matters because most users should understand the lower-level behavior before relying on the higher-level batch workflows.

## Assumptions To Keep In Mind

Even though the repository is being generalized, the code still reflects modeling assumptions that may or may not fit a future use case. A public user should review at least:

- material property definitions
- absorption assumptions
- boundary conditions
- spatial and temporal resolution choices
- solver tolerances and stopping criteria
- how outputs are interpreted physically

## Recommended Direction For Future Cleanup

If this becomes an active public repository, the best long-term direction would be to shift from research-script organization toward clearer public structure:

- common configuration conventions across solvers
- fewer duplicated helper patterns
- consistent naming and units in outputs
- documented example cases with expected behavior
- a cleaner distinction between reusable solver functions in `src/` and editable walkthroughs in `examples/`
- explicit agent-facing guidance and stable edit surfaces for future automated maintenance