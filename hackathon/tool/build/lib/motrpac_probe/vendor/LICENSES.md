# Vendored JavaScript

These files are inlined into reports that contain interactive (Vega-Lite) items,
so that reports work offline with no CDN dependency. They are unmodified
minified builds downloaded from jsDelivr on 2026-09-26.

| File | Package | Version | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| vega.min.js | vega | 5.33.1 | 515242 | 463f3db6a40b20e9747b4ed38f37ed0add508838f9141b1cf8366784b07b30c8 |
| vega-lite.min.js | vega-lite | 5.23.0 | 252198 | 58c27358e26f2d319cf62f45bc17a4c8362f08645001df2ec8d341eee4097c7f |
| vega-embed.min.js | vega-embed | 6.29.0 | 60630 | 12d02acfbe3ec59ef9a37dd4822a2e04e2961b5bbb671bbe661d2221715b99da |

Sources:

- https://cdn.jsdelivr.net/npm/vega@5.33.1/build/vega.min.js
- https://cdn.jsdelivr.net/npm/vega-lite@5.23.0/build/vega-lite.min.js
- https://cdn.jsdelivr.net/npm/vega-embed@6.29.0/build/vega-embed.min.js

To update: download the new builds to this directory, then update this table
(`sha256sum *.js`).

## Licence

All three packages are released under the BSD 3-Clause License by the
University of Washington Interactive Data Lab (https://github.com/vega).

    Copyright (c) 2015-2023, University of Washington Interactive Data Lab
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its
       contributors may be used to endorse or promote products derived from
       this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

Note for packaging: pyproject.toml `package-data` currently lists only
`templates/*`; add `vendor/*` so non-editable installs ship these files.
