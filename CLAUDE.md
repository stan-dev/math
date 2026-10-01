@AGENTS.md

## Claude Code setup

- **Reuse skill.** Invoke the `stan-math-reuse` skill before writing any new function, overload, trait or helper under `stan/`. It walks through the idiom guide, the `LSP` tool and real call sites, then asks for a REUSE / EXTEND / NEW verdict.
- **Code intelligence (clangd).**
1. `./runClangd.py cdb` writes `compile_commands.json`, which is gitignored. Rerun it after adding headers or changing `make/local`.
2. `/plugin install clangd-lsp@claude-plugins-official` gives you the `LSP` tool (definitions, find-references) across the templates. If it is not available you should ask for it.
3. `./runClangd.py --clean` removes `compile_commands.json` and clangd's index in `.cache/clangd/`.
4. clangd indexes on startup and when a file is opened. After switching branches or pulling, restart the session so it re-indexes, and rerun `./runClangd.py cdb` first if headers were added.
5. `.claude/settings.json` sets `CLANGD_FLAGS=-j=2`. clangd's `-j` sizes three pools (background index, AST builds, preamble builds), so this keeps it to about 6 busy threads. For clangd in other editors, set the same variable or pass `-j=2`.
