# Scripts

This folder contains reusable repository-level entry points that are broader than a single editable example.

## Included entry points

- `Verify_Public_Repo_Smoke.m` - lightweight runtime smoke test for the supported public-facing solver surfaces
- `batch/` - reusable multi-case study runners and parameter sweeps

## Recommended use

- Start with `examples/` when you want a single editable case.
- Use `Verify_Public_Repo_Smoke.m` when you want a quick regression check before publishing or after refactoring.
- Use `scripts/batch/` only after a single-case workflow is behaving the way you expect.

## Notes

- Batch runners should stay reusable and should not become manuscript-specific automation.
- Generated results from scripts belong under `outputs/` and should remain untracked.