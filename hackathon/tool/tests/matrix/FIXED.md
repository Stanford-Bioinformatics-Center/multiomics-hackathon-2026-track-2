- **Crash when RNA–protein ρ is exactly ±1** (cases `size5_secreted_half`, `absent_from_proteomics`, first matrix
  run): the Fisher-z CI bound is clipped at ±0.999999, so the error bar went slightly negative and matplotlib raised.
  Fixed in `figures.discordance` and `discord.fig_agreement` (error bars clipped at 0).
- **Empty tables when no ranking unit has enough measured genes** (`empty_overlap_tissue`): section F rendered
  empty tables, which the validator rejects. Fixed in `run.write`: an empty table is replaced by a one-row
  "no rows: …" table. The headline shows "no measured genes" cells, not an error.
- **Negative control above the class null** (found by the Phase 1 suite, confirmed here): the class null matched
  only on MitoCarta / complex / secreted flags, so a high-abundance random set sat above the 98th percentile in
  6 of 52 columns. The class null now also matches on abundance tertile (`nulls.column_null`).
