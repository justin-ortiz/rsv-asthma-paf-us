# RSV-attributable childhood asthma

R code for the estimated burden of childhood asthma at age 6
attributable to infant medically attended RSV lower respiratory tract
infection in the United States.

    rsv_asthma_model.R        model functions
    cache_inputs.R            reads an incidence workbook, writes inputs.rds
    run_analysis.R            assumptions and analysis
    compare_calibrations.R    prints the base case both ways

The incidence workbook is not in this repository. Download
`RSV_incidence_model_simulation_outputs_2026-08-28.zip` from
https://doi.org/10.5281/zenodo.22797654, unzip it, and put
`RSV_incidence_50pct_coverage_base_waning.xlsx` in the same folder as
the R files. Run `cache_inputs.R` once, then `run_analysis.R`. Needs
`readxl`.

Assumptions sit at the top of `run_analysis.R`. To add an effect
estimate, add a row to the `effects` table.

Asthma prevalence in children without infant RSV LRTI is solved from
national prevalence and the odds ratio. `prevalence_unexposed` is the
closed form described in the Methods. `prevalence_unexposed_published`
is the two-pass version the first analysis used, kept so the first
results can be reproduced. `compare_calibrations.R` prints both.

Trials are matched across scenarios within one workbook and can be
subtracted trial by trial. They are not matched across workbooks.

Incidence data: Hutton D, Ortiz J. Monte Carlo simulation outputs from
a US infant RSV lower respiratory tract infection incidence model
under three nirsevimab scenarios. Zenodo, 2026.
doi:10.5281/zenodo.22797654

The original analysis code for this study is at
github.com/MCFitz/RSV-asthma-US.
