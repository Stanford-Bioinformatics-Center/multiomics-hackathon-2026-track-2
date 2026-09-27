# legacy/

Set-aside artifacts that are **not part of the live application** but are kept for
reference. Nothing here is read by the API (`apps/api`), the web app (`apps/web`), or the
engine at runtime.

## `legacy/site/`

The static report gallery previously committed at `hackathon/tool/site/`. It is a
**generated** artifact — the engine's `motrpac_probe site` CLI command re-renders it into
`hackathon/tool/site/` from the example run folders (see
`hackathon/tool/scripts_regenerate_gallery.md`). The committed copy here records a stale
build (an earlier commit) and is retained only for reference; regenerate a fresh gallery
with the CLI rather than relying on these files.

Moved here (via `git mv`, history preserved) so the live `hackathon/tool/` tree carries
only source, not a large regenerable output. See DECISIONS.md ADR-0026.
