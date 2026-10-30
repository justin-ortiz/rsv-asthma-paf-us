# RSV-attributable childhood asthma

R code for the estimated burden of childhood asthma at age 6
attributable to infant medically attended RSV lower respiratory tract
infection in the United States.

    rsv_asthma_model.R        model functions
    cache_inputs.R            reads an incidence workbook, writes an .rds file
    run_analysis.R            assumptions and analysis
    compare_calibrations.R    compares the two-pass and closed-form calibrations

## Data

The incidence workbooks are not in this repository. Download
`RSV_incidence_model_simulation_outputs_2026-08-28.zip` from
https://doi.org/10.5281/zenodo.22797654 (version 1.0, CC BY 4.0) and
unzip it. The three workbooks we used are:

    RSV_incidence_50pct_coverage_base_waning.xlsx
    RSV_incidence_50pct_coverage_postlicensure_waning.xlsx
    RSV_incidence_90pct_coverage_base_waning.xlsx

In the zip each name ends in `.xlsx.xlsx`. Remove one `.xlsx` from each
file name.

## To run

1. Download this repository as a zip file from the green Code button
   and unzip it.
2. Put the three workbooks in the same folder as the R files.
3. In RStudio, open `cache_inputs.R`. Set the working folder with
   Session, Set Working Directory, To Source File Location.
4. Install the reader once: `install.packages("readxl")`
5. Run `cache_inputs.R` once for each workbook. Lines 14 and 15 name the
   workbook and the file it writes. Edit both lines, save, and run.

        50% coverage, base waning          inputs.rds
        50% coverage, post-licensure       inputs_postlicensure.rds
        90% coverage, base waning          inputs_90pct.rds

6. Run `run_analysis.R`. It prints the results and writes seven CSV files
   to the working folder. A new run overwrites them. The files in this repository are the results we report.

Assumptions sit at the top of `run_analysis.R`. To add an effect
estimate, add a row to the `effects` table.

## Method

Asthma prevalence in children without infant RSV LRTI is solved from
national prevalence and the odds ratio. `prevalence_unexposed` is the
closed form described in the Methods. `compare_calibrations.R` compares
it with the two-pass calculation used in the first analysis, at the same
base-case inputs.

Trials are matched across scenarios within one workbook, so scenarios
can be subtracted trial by trial. Trials are not matched across
workbooks.

Incidence data: Hutton D, Ortiz J. Monte Carlo simulation outputs from
a US infant RSV lower respiratory tract infection incidence model
under three nirsevimab scenarios. Zenodo, 2026.
doi:10.5281/zenodo.22797654

The method follows the original analysis code for this study by Meagan Fitzpatrick, github.com/MCFitz/RSV-asthma-US.
