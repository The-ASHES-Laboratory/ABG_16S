# Contributing

This repository is part of an active research program in the ASHES Laboratory, NC State University. Contributions are welcome via pull request.

## Issues

Use GitHub Issues for:

- Bug reports (with a minimal reproducible example: QIIME 2 version, command, and the relevant log output)
- Questions about analysis parameters, DADA2 truncation lengths, or rarefaction depth choices
- Requests to extend the pipeline to additional 16S variable regions or alternative classifiers
- Documentation gaps

Tag issues with one of `bug`, `enhancement`, `documentation`, `parameters`, `qiime2-version`.

## Pull requests

1. Fork and create a topic branch named `feature/short-description`.
2. Confirm that your changes run end-to-end on the QIIME 2 versions documented in the README (currently tested with 2022.2 and amplicon-2024.10).
3. Capture the active QIIME 2 environment (`qiime info`) in every new wrapper script so reproducibility tracking is preserved.
4. Update `CHANGELOG.md` with a brief entry under `[Unreleased]`.
5. Open a PR against `main`.

## Coding style

- Bash wrappers: `set -euo pipefail` at the top, CLI parsing via `getopts` or `argparse-bash`, full `--help` support.
- Python: 3.10+, Black formatting, `argparse` for CLI.
- R: tidyverse style, `argparse` for CLI.
- All scripts capture their parameters and the QIIME 2 environment to a log file in the output directory.

## Reproducibility

- New scripts must accept all inputs via CLI flags; no hardcoded paths or sample IDs.
- Random seeds (DADA2, FastTree, classifier training) must be settable via `--seed`.
- New analyses that depend on QIIME 2 plugins or external tools should document the version range tested.

## QIIME 2 version compatibility

This pipeline targets the QIIME 2 amplicon distribution. When QIIME 2 releases break upstream behavior (for example, the DADA2 → FastTree → midpoint rooting handoff in `07_build_phylogeny.sh`), open an issue tagged `qiime2-version` with the version, the failure mode, and a proposed fix.

## Community-engaged research note

This software supports research conducted in collaboration with descendant-community partners at the New York African Burial Ground. Contributions that propose extending the pipeline to additional descendant lineages or burial contexts should include a brief note on the community-engagement framework that will support the proposed extension.
