# Contributing

Thanks for your interest in improving this repository.

## Good contributions

- bug fixes tied to a clear issue
- portability improvements
- example and documentation clarification
- small solver-interface improvements that preserve current behavior
- validation additions that stay generic rather than manuscript-specific

## Before opening a larger change

- prefer starting from the canonical examples in `examples/`
- keep reusable solver logic in `src/`
- keep generated outputs out of version control
- preserve the shared result-contract fields where practical

## Suggested workflow

1. make a focused change
2. run the relevant example or `scripts/Verify_Public_Repo_Smoke.m` when possible
3. update `README.md` or `docs/` if user-facing behavior changed
4. keep pull requests small enough that solver behavior changes are easy to review

## Scope guidance

This repository stays centered on reusable pulsed-laser thermal modeling workflows. Changes that mainly serve one historical study, or one manuscript-specific validation campaign, generally belong outside the public-facing copy.