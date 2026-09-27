## Summary

The Track 2 "Omic Discordance Explained" work. The **MoTrPAC Explorer** is the default view
(any molecule list in → every matching MoTrPAC comparison out, with a fused RNA↔protein
co-view). Alongside it: a read-only discordance API over committed demo outputs, a live query
builder with a mapping-preview confirmation gate, generalized upload + a shared results table,
retired mock views, a README rewritten to the hackathon's eight-section scaffold anchored to
Track 2, a fresh-clone bootstrap, and automated verification wired into CI. This PR also relocates
the stale, regenerable `hackathon/tool/site` gallery to `legacy/site` (ADR-0026).

Honest framing is preserved throughout: results are cross-cohort association findings, **not** a
PAH treatment-efficacy claim; cohorts are analyzed separately and never merged; the directed
headline (male rat SKM-GN protein, 8 weeks) is **q = 0.0584 and not significant** under the frozen
52-column BH family. All statistics originate in the engines; the API/web layers only serve, shape,
and render.

## What was tested

The verification battery is documented in VERIFICATION.md. Last recorded green state: engine
**123 passed / 34 skipped**, `apps/api` **81 passed**, web **38 passed** + `pnpm run build` OK,
httpx API smoke (headline q = 0.0584 ± 0.0001; discordance per-class counts 1/0/285/5642;
metabolomics null result), Playwright e2e (every live view renders real data; browser
point→evidence→source chain), and the golden-hash guard on `analysis.py` (**OK** —
`1c1089723ec7c9d67c4c0f7ab51262717f80d95322c7da8df60187a06cab7372`). Note: the full local battery
was not re-run immediately before this PR due to a deadline; CI runs it on this PR.

## Deferred items

- Live recompute of the discordance catalog (the adapter serves committed demo outputs only).
- Option-2 engine live metabolite scoring.
- Team actions before final submission: insert Explorer/Discordance screenshots (README §7) and
  confirm per-person contributor roles (README §8).
- Other follow-ons tracked in HANDOFF_NEXT_FEATURES.md.
