# RSV-attributable childhood asthma

R code for the estimated burden of childhood asthma at age 6
attributable to infant medically attended RSV lower respiratory tract
infection in the United States.

    rsv_asthma_model.R        model functions
    cache_inputs.R            reads an incidence workbook, writes inputs.rds
    run_analysis.R            assumptions and analysis
    compare_calibrations.R    prints the base case both ways

The incidence workbooks are not in this repository. Download
`RSV_incidence_model_simulation_outputs_2026-08-28.zip` from
https://doi.org/10.5281/zenodo.22797654, unzip it, and put these three
workbooks in the same folder as the R files:

    RSV_incidence_50pct_coverage_base_waning.xlsx
    RSV_incidence_50pct_coverage_postlicensure_waning.xlsx
    RSV_incidence_90pct_coverage_base_waning.xlsx

Run `cache_inputs.R` once for each workbook. Lines 14 and 15 name the
workbook and the file it writes:

    50% coverage, base waning          inputs.rds
    50% coverage, post-licensure       inputs_postlicensure.rds
    90% coverage, base waning          inputs_90pct.rds

Then run `run_analysis.R`. Needs `readxl`. It prints the results and
writes them to CSV files.

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
