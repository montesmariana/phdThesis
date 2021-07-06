library(tidyverse)
library(flextable)
library(kableExtra)
library(colorblindr)


# automatically create a bib database for R packages ----
knitr::write_bib(c(
  .packages(), 'bookdown', 'knitr', 'rmarkdown',
  'cluster', 'vegan', 'Rtsne', 'dbscan', 'umap'
), 'assets/bib/packages.bib')

# utils ----
na2zero <- function(x) {
  x[is.na(x)] <- "-"
  x
}

# load data ----
d <- readRDS(here::here("assets", "data.rds"))
medoid_data <- read_tsv(here::here("assets", "medoid_data.tsv"), col_types = cols())

# plots ----
plotWithCws <- function(lname, mnum){
  m <- d[[lname]]$medoidCoords[mnum]
  cw_data <- medoid_data %>% 
    filter(lemma == lname, model == names(m)) %>% 
    select(cluster, cloud_type, top_Cw, topFscore)
  plotCoords(m[[1]], cw_data)
}

plotCoords <- function(medoid, cw_data = tibble()) {
  df <- medoid$coords %>% 
    mutate(
      cluster = if_else(cluster == "0", NA_character_, as.character(cluster))
    ) %>% filter(!sense %in% c("none", "remove"))
  if (nrow(cw_data) > 0) {
    cw_data <- cw_data %>%
      mutate(
        cw = paste0(top_Cw, " (", round(topFscore, 2), ")"),
        cluster = as.character(cluster)
      ) %>% 
      select(cluster, cw) %>% deframe()
    df <- df %>% mutate(cluster = cw_data[cluster])
  } else {
    df <- df %>% mutate(cluster = fct_reorder(cluster, as.numeric(cluster)))
  }
  
  ggplot(df) +
    geom_point(aes(x = model.x, y = model.y, color = cluster, shape = sense, alpha = eps),
               size = 3) +
    theme_void() +
    scale_color_OkabeIto(darken = 0.1, use_black = TRUE, na.value = "#9b9c9f") +
    coord_fixed() +
    scale_alpha(range = c(1, 0), guide = "none") +
    scale_shape(guide = "none") #+
  # theme(legend.position = "bottom")
}
