library(tidyverse)
library(flextable)



# automatically create a bib database for R packages ----
knitr::write_bib(c(
  .packages(), 'bookdown', 'knitr', 'rmarkdown',
  'cluster', 'vegan', 'Rtsne', "dbscan"
), 'assets/bib/packages.bib')
