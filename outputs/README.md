# Outputs

This folder is the default destination for generated text files, figures, and batch-study results.

## Version-control policy

- keep this folder in the repository so scripts have a predictable writable location
- do not commit generated run results, smoke-test artifacts, or machine-specific output snapshots
- preserve `.gitkeep` and this file so the folder remains visible in a clean clone

If you want to inspect current solver behavior, generate fresh outputs locally by running the examples or `scripts/Verify_Public_Repo_Smoke.m`.