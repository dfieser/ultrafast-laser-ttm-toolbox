# Agent Guidance

This file gives repository-specific guidance to coding agents and automated assistants working in this MATLAB repository.

## Primary Goal

Preserve the repository as a clean, reusable public-facing solver codebase. Prefer changes that improve clarity, portability, and usability for a general research audience without reintroducing manuscript-specific sprawl.

Agent use is a supported feature of this repository. Changes should make the codebase easier for an automated assistant to understand, navigate, and modify safely.

## Source Of Truth

- `src/` is the main code surface for reusable solver logic.
- `examples/` contains editable single-run walkthroughs and should stay lightweight.
- `scripts/batch/` contains reusable orchestration scripts, not manuscript-specific project automation.
- `docs/` explains repository scope, context, and curation intent.
- `outputs/` is for generated artifacts and should not accumulate committed result data.

## Recommended Starting Points For Agents

- Start in `examples/` when the task is about user workflow, tutorial polish, or the fastest runnable entry point.
- Start in `src/` when the task is about solver behavior, config interfaces, returned results, or output semantics.
- Start in `scripts/batch/` when the task is about repeated studies, multi-case orchestration, or summary generation.
- Start in `docs/` or `README.md` when the task is about discoverability, positioning, scope, or onboarding.

## Stable Edit Surfaces

- `src/Depth_Profile_Solver.m`, `src/Radial_Profile_Solver.m`, `src/Inversion_Quantifier.m`, `src/Single_Pulse_Visualizer.m`, `src/Surface_Point_Solver.m`, and `src/Scanning_Beam_Solver.m` are supported solver surfaces.
- `src/RepoUtils.m` is the shared utility surface for batch-script helpers.
- `examples/*.m` should remain thin wrappers or starter cases, not solver forks.
- `examples/Example_Scanning_Beam_Baseline.m` is the canonical scanning example; `examples/Scanning_Beam_Single_Run.m` is a legacy alias.
- `scripts/batch/*.m` should reuse shared patterns rather than embedding duplicated helper logic.

## Shared Result Contract

The main solver entry points now expose a shared minimal contract for agent-friendly automation. Preserve these fields when practical:

- `solver`
- `solverId`
- `contractVersion`
- `material`
- `outputFile`
- `outputDir`
- `inputConfig`

Where already present, also preserve `nPulses` and `wallTime_s`.

The machine-readable companion note for this contract lives in `docs/RESULT_CONTRACT_SCHEMA.json`.

## What To Preserve

- relative-path portability
- separation between reusable code and generated outputs
- separation between reusable solver code and editable examples
- generalized wording in documentation where possible
- the curation decision to avoid manuscript-bound repo structure
- version consistency: the `VERSION` file, `CITATION.cff`, and `.zenodo.json` must state the same version; bumping `VERSION` on `main` triggers a GitHub release and a new Zenodo version DOI

## What To Avoid

- reintroducing archived backups, duplicate source copies, or result dumps
- turning the repository back into a paper-specific validation bundle
- hardcoded local paths
- unnecessary large-scale rewrites that change solver behavior without need
- formatting-only edits across unrelated files

## Editing Priorities

When making changes, prefer this order:

1. fix portability or path issues
2. improve public-facing documentation and discoverability
3. clarify configuration surfaces and default behavior
4. make minimal solver-code fixes tied to a concrete problem

## Validation Expectations

- Prefer validating one touched file or one touched slice at a time.
- For documentation changes, check cross-references and file naming consistency.
- For example scripts, ensure they stay thin and call reusable code rather than duplicating solver internals.
- For solver changes, prefer preserving returned results and relative-path behavior.
- Use `docs/VALIDATION_GUIDE.md` as the default lightweight manual validation workflow when runtime execution is available.

## Documentation Expectations

If you add or reorganize code, update the corresponding documentation when needed:

- update `README.md` for user-facing entry points or workflow changes
- update `docs/PROJECT_CONTEXT.md` if repository scope or intent changes
- update `docs/REPOSITORY_NOTES.md` if curation choices or publishing guidance changes
- update `AGENTS.md` if the expected repository structure for future automation changes

## Batch Script Expectations

Batch scripts should remain reusable examples of parameter-study automation. Do not add one-off experiment scripts to `scripts/batch/` unless they are clearly generalized and documented as reusable.

Example scripts should live under `examples/` rather than `src/` when they exist mainly to be edited and run directly by a user.

Canonical starter examples should stay short and should prefer building a config struct and calling a solver rather than copying solver internals.

## Output Handling

Keep generated artifacts under `outputs/`. Do not commit transient run results, figure dumps, or machine-specific output snapshots unless the task explicitly requires a small curated example artifact.

## If You Need To Expand The Repo

Prefer additions that make the repo easier for an outside user to understand:

- example configurations
- concise usage documentation
- clearer parameter conventions
- lightweight validation examples that are generic rather than manuscript-bound
- agent-facing repo maps, workflows, or maintenance guidance when they improve automation reliability

If a proposed change mainly serves one historical study, it probably does not belong in this public-facing repository.