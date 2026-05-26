# Validation Guide

This document defines a lightweight validation workflow for the staged public repository.

It is intentionally simple. The goal is not exhaustive scientific verification of every physics assumption, but a repeatable way to confirm that the public-facing solver surfaces and example entry points still behave coherently after cleanup work.

## Validation Strategy

Use the canonical examples under `examples/` as the first-pass smoke tests.

They are suitable because they:

- exercise the supported public entry points
- write into repo-relative `outputs/`
- avoid manuscript-specific data dependencies
- provide a stable surface for both human users and coding agents

## Recommended Validation Order

1. `examples/Example_Surface_Point_Baseline.m`
2. `examples/Example_Depth_Profile_Baseline.m`
3. `examples/Example_Radial_Profile_Baseline.m`
4. `examples/Example_Scanning_Beam_Baseline.m`

## What To Check

### General checks for every run

- the script runs without path errors
- outputs are written into `outputs/` rather than a machine-local folder
- the returned `results` struct is created
- the shared result-contract fields are present where expected

Shared fields to check:

- `solver`
- `solverId`
- `contractVersion`
- `material`
- `outputFile`
- `outputDir`
- `inputConfig`

### Surface-point example

Check that:

- `results.peakTe_C` and `results.peakTl_C` are present
- `results.TeqVals_C` and `results.TresidVals_C` are present
- the primary text output file exists in the configured output folder

### Depth-profile example

Check that:

- `results.nPulses` is present
- `results.wallTime_s` is present
- the solver writes its output text file to the configured output folder
- the returned temperature and inversion arrays still exist if downstream tools use them

### Radial-profile example

Check that:

- `results.rGrid_um` is present
- `results.finalRadialProfile_C` is present
- the configured solve mode still matches the example expectation

### Scanning-beam example

Check that:

- `results.peakT_history` is present
- `results.outputFile` is set
- generated files stay within the repo `outputs/` tree

## Manual MATLAB Spot Check

After running an example, a minimal sanity check can be done in MATLAB with statements like:

```matlab
isfield(results, 'solverId')
isfield(results, 'outputFile')
exist(results.outputFile, 'file')
```

For pulse-based solvers, also check:

```matlab
isfield(results, 'nPulses')
```

## What This Guide Does Not Guarantee

This guide does not prove:

- full physical correctness for every material and laser regime
- manuscript-level validation against experiments
- numerical optimality of every discretization choice

It is a repo-maintenance guide, not a substitute for scientific validation.

## When To Extend This Guide

Extend this guide when:

- a new public example becomes a supported entry point
- a solver's returned-result contract changes
- runtime checks become available and you want stronger regression coverage