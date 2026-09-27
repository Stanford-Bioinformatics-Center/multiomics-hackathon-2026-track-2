// # Feature: finalize-discordance-mvp, Property 15
//
// Property 15 — retired-view redirect resolution (Requirement 6.2).
//
// Direct navigation to a retired mock design view (architecture/01, technical/02, workflow/03,
// dashboard/04, slide/08) — or to any unknown/garbage identifier — must resolve to the last-used
// LIVE view when one exists this session, otherwise to the default view (Explorer after the motrpac-explorer merge). A retired
// identifier is NEVER the resolved result; the resolver only ever yields a live, navigable view.
//
// A live identifier always resolves to itself, regardless of the last-used value.
//
// This exercises the pure `resolveView(requested, lastLiveView)` resolver directly (no DOM needed):
// it is deterministic and side-effect free, so fast-check can drive ≥100 iterations cheaply.

import { describe, expect, test } from "vitest";
import fc from "fast-check";
import { resolveView, RETIRED_VIEWS, NAV_VIEWS, DEFAULT_VIEW, isLiveView, type View } from "./App";

// ─────────────────────────── Fixed identifier sets ───────────────────────────
const LIVE_IDS: View[] = NAV_VIEWS.map((entry) => entry.id);
const LIVE_ID_SET = new Set<string>(LIVE_IDS);
const RETIRED_ID_SET = new Set<string>(RETIRED_VIEWS);

// ─────────────────────────── Arbitraries ───────────────────────────
const retiredId = fc.constantFrom(...RETIRED_VIEWS);
const liveId = fc.constantFrom(...LIVE_IDS);

// Arbitrary strings that could be typed into a URL. `fc.string()` can occasionally produce a value
// that happens to equal a live or retired id; that is fine — the assertions below cover every case.
const garbageString = fc.string();

// `requested`: retired ids ∪ live ids ∪ arbitrary strings.
const requestedArb = fc.oneof(retiredId, liveId, garbageString);

// `lastLiveView`: live ids ∪ retired ids ∪ null/undefined ∪ arbitrary strings.
const lastLiveArb = fc.oneof(
  liveId,
  retiredId,
  fc.constant(null),
  fc.constant(undefined),
  garbageString,
);

describe("Property 15: retired-view redirect resolution", () => {
  test("resolveView never returns a retired view and always returns a live view", () => {
    fc.assert(
      fc.property(requestedArb, lastLiveArb, (requested, lastLiveView) => {
        const resolved = resolveView(requested, lastLiveView);

        // Invariant: the result is ALWAYS a live, navigable view id and NEVER a retired one.
        expect(LIVE_ID_SET.has(resolved)).toBe(true);
        expect(RETIRED_ID_SET.has(resolved)).toBe(false);
        expect(isLiveView(resolved)).toBe(true);

        if (isLiveView(requested)) {
          // (2) A live requested id resolves to itself, regardless of the last-used value.
          expect(resolved).toBe(requested);
        } else {
          // (1) + (3) Retired or unknown requested id: resolves to last-used-if-live else default.
          if (isLiveView(lastLiveView)) {
            expect(resolved).toBe(lastLiveView);
          } else {
            expect(resolved).toBe(DEFAULT_VIEW);
          }
        }
      }),
      { numRuns: 100 },
    );
  });

  test("retired requested id resolves to last-used live view, else the default view — never retired", () => {
    fc.assert(
      fc.property(retiredId, lastLiveArb, (requested, lastLiveView) => {
        const resolved = resolveView(requested, lastLiveView);

        expect(RETIRED_ID_SET.has(resolved)).toBe(false);
        expect(LIVE_ID_SET.has(resolved)).toBe(true);

        if (isLiveView(lastLiveView)) {
          expect(resolved).toBe(lastLiveView);
        } else {
          expect(resolved).toBe(DEFAULT_VIEW);
        }
      }),
      { numRuns: 100 },
    );
  });

  test("live requested id resolves to itself regardless of last-used value", () => {
    fc.assert(
      fc.property(liveId, lastLiveArb, (requested, lastLiveView) => {
        const resolved = resolveView(requested, lastLiveView);
        expect(resolved).toBe(requested);
        expect(LIVE_ID_SET.has(resolved)).toBe(true);
      }),
      { numRuns: 100 },
    );
  });

  test("arbitrary unknown/garbage requested id behaves like a retired id", () => {
    fc.assert(
      fc.property(garbageString, lastLiveArb, (requested, lastLiveView) => {
        const resolved = resolveView(requested, lastLiveView);

        expect(RETIRED_ID_SET.has(resolved)).toBe(false);
        expect(LIVE_ID_SET.has(resolved)).toBe(true);

        if (isLiveView(requested)) {
          // A garbage string that happens to equal a live id resolves to itself.
          expect(resolved).toBe(requested);
        } else if (isLiveView(lastLiveView)) {
          expect(resolved).toBe(lastLiveView);
        } else {
          expect(resolved).toBe(DEFAULT_VIEW);
        }
      }),
      { numRuns: 100 },
    );
  });
});
