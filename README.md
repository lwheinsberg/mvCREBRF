Code Repository for Multivariate Analysis
================
Jerry Z. Zhang, Lacey W. Heinsberg, and Daniel E. Weeks<br/>Department
of Human Genetics<br/>University of Pittsburgh



# Sample analysis code

Copyright 2022, University of Pittsburgh. All Rights Reserved.  
License: CC BY-SA 3.0
([link](https://creativecommons.org/licenses/by-sa/3.0/))

This repository contains example R markdown analysis code used for
“Multivariate analysis of a missense variant in CREBRF reveals
association with measures of regional and total adiposity in people of
Polynesian ancestries” 
([published](https://onlinelibrary.wiley.com/doi/10.1002/gepi.22508) and
[pre-print](https://www.medrxiv.org/content/10.1101/2022.09.08.22279720v1)).

The purpose of this paper was to improve our understanding of a missense
variant in CREBRF, rs373863828, that is paradoxically associated with
increased body mass index, but a decreased odds of diabetes in
Polynesian groups. Specifically, we examined the association between
rs373863828 and a panel of correlated anthropometric and lipid profile
phenotypes (body mass index, weight, height, HDL cholesterol, net
triglycerides, total cholesterol, fat-free mass, fat mass, waist-hip
ratio, abdominal circumference, and hip circumference) using Bayesian
multivariate and network analyses.

The purpose of this repository is to more fully document the details of
our analyses through annotated code.

While the Samoan data analyzed in our paper are available via dbGAP
([link](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000914.v1.p1)),
a synthetic simulated data set has been created for immediate use with
the example analysis code. The synthetic data set was created using
[`synthpop`](https://cran.r-project.org/web/packages/synthpop/index.html)
to emulate correlation structures found within one of the cohorts
studied in the paper. The data are wholly synthetic and contain no
observations from the actual study.

If you are adapting this code to perform the analyses in your own data
set, you will first need to pre-process your data to assure the data
set:

1.  Contains only unrelated participants;

2.  Uses additive coding of genotype(s) (based on the number of minor
    allele copies) as 0, 1, 2;

3.  Includes all correlated phenotypes of interest; and

4.  Includes all covariates that will need to be regressed out (e.g.,
    age, gender, ancestry principal components, etc).

The example analysis code is presented as two R Markdown files: (1)
[`01_phenotypes_mvBIMBAM_run.Rmd`](./01_phenotypes_mvBIMBAM_run.Rmd) and
(2) [`02_mvBayesNet_ReadRDS.Rmd`](./02_mvBayesNet_ReadRDS.Rmd).
Alternatively, we have also provided typeset
[`.md`](./01_phenotypes_mvBIMBAM_run.md) and
[`.pdf`](./01_phenotypes_mvBIMBAM_run.pdf) versions of the files. Note
that the `.md` versions contain results and are nicely readable online
within the GitHub website.

In `01_phenotypes_mvBIMBAM_run.Rmd`, we perform multivariate association
analysis between rs373863828 and phenotypes of interest using
[mvBIMBAM](https://github.com/heejungshim/mvBIMBAM), a Bayesian approach
for joint genetic association analysis of multiple related phenotypes.

In `02_mvBayesNet_ReadRDS.Rmd`, we further study and plot the structure
of the relationships between rs373863828 and phenotypes of interest
using [bnlearn](https://www.bnlearn.com/), an R package developed to
perform Bayesian network analyses. Note, the code within
`02_mvBayesNet_ReadRDS.Rmd` was adapted from
<http://www.bnlearn.com/research/genetics14/> under the Creative Commons
Attribution-Share Alike License
(<https://creativecommons.org/licenses/by-sa/3.0/>) which accompanies
the paper “Multiple Quantitative Trait Analysis Using Bayesian Networks”
by Scutari, Howell, Balding, Mackay (Genetics, 2014).

# Dependencies

The first R Markdown, `01_phenotypes_mvBIMBAM_run.Rmd`, relies on
calling [mvBIMBAM](https://github.com/heejungshim/mvBIMBAM) using
`base::system()`. This requires `bimbam` to be added to `$PATH` on
Unix-like systems or `%PATH%` on Windows based systems as described at
<https://github.com/heejungshim/mvBIMBAM>.

A compiler will be needed to compile mvBIMBAM for usage on your local
machine.

The code also relies on R package dependencies listed within the “Load
libraries” sections of both R Markdown files.

# Execution

R Studio can be used to execute the .Rmd code using the `Knit` function
from the quick bar.

Alternatively, this code can be executed without RStudio relying on a
standalone Pandoc installation. With Pandoc installed,
`rmarkdown::render()` can be used to execute and knit the .Rmd code.

The synthetic data created for use with this code has been uploaded in
this repository in the file
[`git_synth_data.rds`](./git_synth_data.rds).

The R Markdown `01_phenotypes_mvBIMBAM_run.Rmd` should be run first, as
it creates a quantile normalized, covariate-adjusted, and outlier-free
version of the synthetic data set and computes Bayes factors that are
needed as input in the `02_mvBayesNet_ReadRDS.Rmd` R Markdown file.

The R Markdown `02_mvBayesNet_ReadRDS.Rmd` should be run second.

Note that in `02_mvBayesNet_ReadRDS.Rmd`, the code for
`f_run_plot_graph` assumes there are 10 processors available for
parallel processing (ncluster = 10). Depending on the capabilities of a
user’s machine, this may need to be adjusted.

# Example run and interpretation

Example analysis code has been carefully annotated in the R Markdown
files

`01_phenotypes_mvBIMBAM_run.Rmd`

and

`02_mvBayesNet_ReadRDS.Rmd`.

These files also include notes about interpretation of the results.

Compiled versions of these, complete with output and figures, can be found in these Markdown files:

[01_phenotypes_mvBIMBAM_run.md](https://github.com/lwheinsberg/mvCREBRF/blob/master/01_phenotypes_mvBIMBAM_run.md)

and

[02_mvBayesNet_ReadRDS.md](https://github.com/lwheinsberg/mvCREBRF/blob/master/02_mvBayesNet_ReadRDS.md)

or in the corresponding PDF files.

# Contact information

If you have any questions or comments, please feel free to contact us.

Daniel E. Weeks: <weeks@pitt.edu>
