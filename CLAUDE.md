# MatchY — Claude Code instructions

## Workflow

**Commit before making changes.**
Before starting any new task or set of changes, create a git commit of the current working state (if there are uncommitted changes). This ensures every logical change is captured separately and makes it easy to revert or review individual pieces of work.

## Build & deploy

Run the build script from the repo root:

```
.\build.ps1
```

Artifacts are written to `C:\cargo-target\matchy\release\bundle\`.
