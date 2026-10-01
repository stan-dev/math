---
name: stan-math-reuse
description: Use BEFORE writing any new function, overload, type trait, require_* alias, helper struct, or loop over Eigen matrices / std::vector in stan/math. Finds existing Stan Math symbols (require_*, is_*, value_type_t, scalar_type_t, return_type_t, plain_type_t, ref_type_t, to_ref, make_holder, arena_t, reverse_pass_callback, make_callback_var, apply_scalar_unary, apply_vector_unary, partials_propagator, check_*, value_of, log_sum_exp) and decides whether to reuse, extend, or write something new.
when_to_use: "add a function", "implement", "write a helper", "new trait", "SFINAE", "enable_if", "vectorize", "loop over", "gradient for", "new distribution", "_lpdf", before creating any new .hpp under stan/math.
allowed-tools: Read Grep Glob LSP Bash(grep *) Bash(git grep *) Bash(ls *)
---

# Stan Math reuse check

clangd setup: !`test -f "${CLAUDE_PROJECT_DIR}/compile_commands.json" && echo "compile_commands.json present" || echo CDB_MISSING`

Run this procedure **per task, immediately before writing code**, even if you
explored the repo earlier in the session. Agents routinely rewrite helpers
they read a few turns ago; the fix is to make the decision explicitly.

## Procedure

1. **Classify the need.** Write one line:
   `Need: <behavior>; inputs <types>; module <prim|rev|fwd|opencl>; kind <trait | SFINAE constraint | math function | vectorization | reverse-mode plumbing | error check | numerics>`

2. **Idiom guide.** Grep `doxygen/contributor_help_pages/idioms.md` with 2–3
   synonyms for the behavior, e.g.
   `grep -n -i -E "evaluate|expression" doxygen/contributor_help_pages/idioms.md`

3. **Look up candidates with the `LSP` tool.** Every operation needs a
   `filePath`, `line` and `character`; for `workspaceSymbol` any header at
   line 1, character 1 works.
   - `workspaceSymbol` with a name fragment finds every match across modules,
     e.g. `require_eigen` lists all the related aliases and `log_sum_exp`
     lists the prim, rev, fwd and opencl overloads.
   - `documentSymbol` on a candidate header lists every overload with its
     parameter types; `hover` on one shows its full doc comment.
   - Search by behavior in the doc comments:
     `grep -rli "log of the sum" stan/math/prim/fun`
   - If the LSP tool is missing or returns nothing, follow
     [grep-fallback.md](grep-fallback.md). Results stay empty until clangd
     has loaded its index, which happens after the first file is opened.
   - If the status above says `CDB_MISSING`, tell the user once: "Run
     `./runClangd.py cdb` so clangd has the right compile flags."

4. **Check real usage.** For the top 1–3 candidates:
   - Run LSP findReferences, or `git grep -l <sym> -- stan/math | head` if LSP
     is unavailable.
   - Read one real call site, so you use the symbol the way the library does.

5. **Decide, and write the decision down before any code:**

   ```
   Considered: `sym` (header) — why it fits / doesn't
   Verdict: REUSE `sym` | EXTEND `sym` (new overload/specialization in <header>) | NEW
   ```

   A NEW verdict must name the closest existing symbol and the concrete
   behavior it lacks.

6. **If EXTEND or NEW,** follow [conventions.md](conventions.md): layering,
   `check_*`, `to_ref`, arena types in rev, and registering the header.

## When to use a subagent

- Do the lookup inline; it is a handful of LSP queries and greps.
- Spawn an Explore subagent only when:
  - the need spans three or more modules, or
  - there is no LSP and the grep output would flood the context.
- Ask the subagent for `symbol | header | fit` lines, not source code.
