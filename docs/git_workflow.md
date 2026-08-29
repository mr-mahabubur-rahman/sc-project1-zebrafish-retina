# Git and GitHub

Repository name is yours to choose; the README examples assume
`sc-project1-zebrafish-retina` under the account `mr-mahabubur-rahman`.

## Check what would be committed BEFORE staging

```bash
git status --short
git status --ignored --short | head -40     # confirm data/ and results/ are ignored
du -sh data results                          # what you are keeping out
```

Never committed: `data/` (large, redistributable from GEO GSE202212) and
`results/*.h5ad` (regenerable; individual files exceed GitHub's 100 MB limit).
GitHub warns above 50 MB and rejects above 100 MB, and a large blob stays in
history even after deletion.

## First push

```bash
git init
git branch -M main
git add .gitignore README.md LICENSE pyproject.toml requirements.txt
git add scripts/ tools/ docs/ notebooks/ tables/ figures/
git status --short          # read this list before committing
git commit -m "Scanpy workflow for zebrafish retina regeneration scRNA-seq"

git remote add origin https://github.com/mr-mahabubur-rahman/sc-project1-zebrafish-retina.git
git push -u origin main
```

## Milestone tags

```bash
# Milestone 1 -- due 31 Aug 2026
git add -A && git commit -m "Draft report and complete pipeline"
git tag -a v1.0-peerreview -m "Milestone 1: draft report and pipeline for peer review"
git push origin main --tags        # a plain `git push` does NOT send tags

# Milestone 3 -- due 9 Sep 2026
git add -A && git commit -m "Final report incorporating peer review"
git tag -a v2.0-final -m "Milestone 2/3: final submission after rebuttal"
git push origin main --tags
```

Verify the tag is actually on the remote — a missing tag is a missed deliverable:

```bash
git ls-remote --tags origin
```

## Commands not to run

- `git push --force` — rewrites remote history; can destroy a reviewer's reference
  point mid-review. Only with a specific reason and after the review window.
- `git reset --hard` before checking `git status` — discards uncommitted work.
- `git clean -fdx` — deletes ignored files, i.e. your `data/` folder.
- `git rm -r --cached .` used casually — can unstage everything and produce a
  confusing diff.

## If a large file was committed by accident

Do not just `git rm` it: the blob stays in history and the push still fails.

```bash
git rm --cached results/04_clustered.h5ad
git commit -m "Remove checkpoint committed in error"
# if already pushed, the blob must be purged from history with git-filter-repo,
# which rewrites commits -- coordinate with anyone who has cloned first
```

Prevention: `git status --short` before every `git add`, and keep `.gitignore`
correct from the first commit.
