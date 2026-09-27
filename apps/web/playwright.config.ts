import { defineConfig, devices } from "@playwright/test";

/**
 * Playwright end-to-end configuration for the Exercise Signature Explorer web app (R9 AC4, AC5).
 *
 * This boots BOTH processes the app needs — the FastAPI API and the Vite/React web dev server — via
 * the `webServer` array below, mirroring the exact commands and port model documented in
 * `RUN_LOCAL.md`:
 *
 *   - API  (uvicorn):  `uvicorn motrpac_probe_service.app:app --host 127.0.0.1 --port <API_PORT>`
 *                      run from `apps/api`, using the Python env that has `apps/api[api]` + the
 *                      mprobe store installed. RUN_LOCAL.md picks 8765 as an arbitrary free port.
 *   - web  (Vite):     `pnpm run dev` from `apps/web`, with `VITE_API_BASE_URL` pointed at the API
 *                      (without it the live views cannot fetch data) and `PORT` set so Vite serves
 *                      on a known origin. RUN_LOCAL.md uses 8443 by default.
 *
 * Ports and the API launch command are overridable via env so the same config works locally, in CI,
 * and when a port is already taken (RUN_LOCAL.md's "pick another" guidance). The repo standardizes on
 * pnpm (R8 AC6), so the web server command uses pnpm.
 *
 * The API needs the mprobe store (`mprobe store fetch`) and the editable installs from RUN_LOCAL.md's
 * prerequisites. `reuseExistingServer` lets a developer boot the two processes by hand (per
 * RUN_LOCAL.md) and just run `pnpm test:e2e` against them.
 */

const API_HOST = process.env.E2E_API_HOST ?? "127.0.0.1";
const API_PORT = process.env.E2E_API_PORT ?? "8765";
const WEB_PORT = process.env.E2E_WEB_PORT ?? "8443";

const API_BASE_URL = `http://${API_HOST}:${API_PORT}/api`;
const WEB_BASE_URL = `http://127.0.0.1:${WEB_PORT}`;

/**
 * Command that starts the API. Overridable via `E2E_API_COMMAND` because the interpreter/venv path
 * differs between a local checkout (RUN_LOCAL.md uses `hackathon/tool/.venv`) and CI (which installs
 * the package into the job's Python). The default assumes `uvicorn` is on PATH for the active env,
 * matching RUN_LOCAL.md once its prerequisite installs + `source .venv/bin/activate` have run.
 */
const API_COMMAND =
  process.env.E2E_API_COMMAND ??
  `uvicorn motrpac_probe_service.app:app --host ${API_HOST} --port ${API_PORT}`;

// The API is launched from apps/api; the web dev server from this directory (apps/web).
const API_CWD = process.env.E2E_API_CWD ?? "../api";

export default defineConfig({
  testDir: "./e2e",
  // Boot can be slow: the API loads the mprobe store and the first comparison run executes the
  // engine. Give both the servers and individual tests generous budgets.
  timeout: 120_000,
  expect: { timeout: 30_000 },
  fullyParallel: false,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 1 : 0,
  workers: 1,
  reporter: process.env.CI ? [["list"], ["html", { open: "never" }]] : [["list"]],

  use: {
    baseURL: WEB_BASE_URL,
    trace: "on-first-retry",
    screenshot: "only-on-failure",
  },

  projects: [
    {
      name: "chromium",
      use: { ...devices["Desktop Chrome"] },
    },
  ],

  // Two servers: the API first (the web app fetches from it), then Vite. Playwright waits for each
  // URL to answer before running tests. `reuseExistingServer` is honored outside CI so a hand-booted
  // pair (RUN_LOCAL.md) is reused rather than double-started.
  webServer: [
    {
      command: API_COMMAND,
      cwd: API_CWD,
      url: `${API_BASE_URL}/health`,
      timeout: 120_000,
      reuseExistingServer: !process.env.CI,
      stdout: "pipe",
      stderr: "pipe",
    },
    {
      command: "pnpm run dev",
      cwd: ".",
      url: WEB_BASE_URL,
      timeout: 120_000,
      reuseExistingServer: !process.env.CI,
      stdout: "pipe",
      stderr: "pipe",
      env: {
        // Point the frontend at the API booted above; without this the live views show a
        // "backend not reachable" state (RUN_LOCAL.md).
        VITE_API_BASE_URL: API_BASE_URL,
        PORT: WEB_PORT,
      },
    },
  ],
});
