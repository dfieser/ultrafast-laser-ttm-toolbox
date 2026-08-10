# Examples

This folder contains small, editable entry-point scripts for common starting workflows. Each example adds `src/` to the MATLAB path, builds a compact config struct, and calls one of the reusable solvers.

## Recommended order

1. `Example_Surface_Point_Baseline.m`
   Start here if you want the simplest pulse-accumulation workflow.
2. `Example_Depth_Profile_Baseline.m`
   Use this for the main 1D depth-resolved solver.
3. `Example_Radial_Profile_Baseline.m`
   Use this when the main question is radial temperature distribution.
4. `Example_Scanning_Beam_Baseline.m`
   Use this for a moving-beam process example.

`Scanning_Beam_Single_Run.m` remains as a legacy alias so older notes and commands still work.

## Editing pattern

Each script is intentionally short:

- edit the values in the `cfg` or `params` struct
- run the script
- inspect the returned `results` struct and files written to `outputs/`

If you need repeated studies rather than one-off examples, use the reusable runners under `scripts/batch/`.
