import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";
import path from "node:path";

/**
 * Vitest config for the web app's component/property tests (jsdom + Testing Library).
 *
 * Kept separate from `vite.config.ts` so the test runner does not pull in the Figma-Make
 * dev/build plugins (site-config, error-overlay replay, HMR fallbacks, kit route). Reusable by
 * every `*.test.tsx` / `*.property.test.tsx` suite (tasks 4.5, 4.6, and later fast-check suites).
 *
 * The existing Node built-in domain test (`src/domain/analysis.test.ts`) is intentionally NOT run
 * by Vitest — it stays on `pnpm run test:domain`. We exclude `*.test.ts` here so Vitest only owns the
 * jsdom/React `*.test.tsx` and `*.property.test.tsx` suites.
 */
export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: { "@": path.resolve(__dirname, "./src") },
  },
  test: {
    environment: "jsdom",
    globals: true,
    include: ["src/**/*.test.tsx", "src/**/*.property.test.tsx"],
    css: false,
  },
});
