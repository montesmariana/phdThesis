library(tidyverse)
library(flextable)
library(kableExtra)
library(colorblindr)
library(sna)
library(GGally)
library(gghalves)
library(ggpubr)
library(cowplot)
library(here)


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

cloud_foto <- function(fname){
  knitr::include_graphics(here("assets", "img", paste0(fname, ".jpg")))
}

html2md <- function(txt) {
  str_replace_all(txt, "<strong>([^<]+)</strong>", "*\\1*") %>%
    str_replace_all("<span class='target'>([^<]+)</span>", "**\\1**") %>% 
    str_remove_all("<sup>[^<]+</sup>")
}

# load data ----
cloud_order <- c("Cumulus", "Stratocumulus", "Cirrus", "Cumulonimbus", "Cirrostratus")
d <- readRDS(here("assets", "data.rds"))
medoid_data <- read_tsv(here("assets", "medoid_data.tsv"), col_types = cols()) %>% 
  mutate(cloud_type = fct_relevel(cloud_type, cloud_order))
cloud_data <- read_tsv(here("assets", "classification.tsv"), col_types = cols()) %>% 
  mutate(cloud_type = fct_relevel(cloud_type, cloud_order))

# plots ----

complex_cor <- function(df, x_coord, y_coord, color_var, xlab = "", ylab = "", colorlab = "", add_abline = FALSE){
  baseplot <- df %>% 
    ggplot(aes(x = {{ x_coord }}, y = {{ y_coord }}, color = {{ color_var }})) +
    geom_point(alpha = 0.5, size = 3) +
    theme_pubr() +
    scale_colour_viridis_d() +
    labs(x = xlab, y = ylab, color = colorlab)
  if (add_abline) baseplot <- baseplot + geom_abline()
  
  xplot <- axis_canvas(baseplot, axis = "x") +
    geom_boxplot(data = df, aes(x = {{ x_coord }}, fill = {{ color_var }}), color = "black", alpha = 0.5) +
    scale_fill_viridis_d()
  yplot <- axis_canvas(baseplot, axis = "y", coord_flip = TRUE) +
    geom_boxplot(data = df, aes(x = {{ y_coord }}, fill = {{ color_var }}), color = "black", alpha = 0.5) +
    scale_fill_viridis_d() +
    coord_flip()
  
  baseplot %>% 
    insert_xaxis_grob(xplot, position = "top") %>%
    insert_yaxis_grob(yplot, position = "right") %>%
    ggdraw()
}

plotBasic <- function(mdata){
  m <- d[[mdata$lemma]]$medoidCoords[[mdata$model]]
  m$coords <- mutate(m$coords, sense = "sense")
  plotCoords(m) + theme(legend.position = 'none')
}

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

network_cws <- function(coords, n = 150) {
  stopifnot(require("sna"))
  cws_count <- coords$cws %>% flatten_chr %>% table %>% sort %>% tail(n)
  all_cws <- names(cws_count)
  cw_freq <- unname(cws_count)
  cws_mtx <- map(all_cws, function(a) {
    A <- filter(coords, map_lgl(cws, has_element, a))
    map_dbl(all_cws, ~nrow(filter(A, map_lgl(cws, has_element, .x))))
  })
  cws_mtx <- flatten_dbl(cws_mtx) %>% matrix(nrow = length(all_cws), dimnames = list(all_cws, all_cws))
  net <- network::network(cws_mtx, ignore.eval = FALSE, names.eval = "weights")
  net %v% "freq" = rowSums(cws_mtx)
  cw_freq <- map_dbl(all_cws, ~max(cws_mtx[.x,]))
  net %v% "cw_freq" = case_when(cw_freq > 10 ~ 6, cw_freq > 1 ~ 4, TRUE ~2)
  set.seed(8541)
  GGally::ggnet2(net, size = "freq", label = TRUE, label.size ="cw_freq",
                 edge.size = "weights", edge.alpha = 0.3) +
    guides(size = "none")
}

# Examples ----
source_names <- c(
  'nrc_handelsblad' = 'NRC Handelsblad',
  'de_standaard' = 'De Standaard',
  'de_morgen' = 'De Morgen',
  'het_laatste_nieuws' = 'Het Laatste Nieuws',
  'volkskrant' = 'Volkskrant',
  'algemeen_dagblad' = 'Algemeen Dagblad',
  'parool' = 'Het Parool',
  'het_nieuwsblad' = 'Het Nieuwsblad'
)
 
sampleCloud <- function(cloudtype) {
  set.seed(8541)
  medoid_data %>% 
    filter(
      (cloudtype %in% maincat & maincat == cloudtype) |
        (!cloudtype %in% maincat & str_detect(maincat, cloudtype)),
      !Hail
    ) %>% slice_sample(n = 1)
}

sampleCtxt <- function(l, mnum = 1, clusn = NULL, stag = NULL, cw = NULL, seed = 8541, n = 1) {
  set.seed(seed)
  m <- d[[l]]$medoidCoords[[mnum]]$coords
  if (is.null(clusn)) clusn <- unique(m$cluster)
  if (is.null(stag)) stag <- unique(m$sense)
  m <- m %>% filter(
    cluster %in% clusn,
    sense %in% stag
  )
  if (!is.null(cw)) m <- m %>% filter(map_lgl(cws, has_element, cw))
  m %>% sample_n(size = 1) %>% 
    separate(`_id`, into = c('lemma', 'pos', 'source', 'line'), sep = '/') %>% 
    mutate(
      numbers = str_replace(source, "[a-z_]+_([0-9_]+)", "\\1") %>% str_replace("_01_", "_"),
      source = source_names[str_remove(source, paste0("_", numbers))]
    ) %>% 
    separate(numbers, into = c("date", "artn"), remove = F) %>% 
    mutate(date = as.Date(date, format = "%Y%m%d"),
           txt = str_glue("{ctxt} (*{source}*, {date}, Art. {artn})")) %>% 
    pull(txt) %>% html2md()
}

# Computations ----

countCTypes <- function(df) {
  df %>% count(cloud_type, Hail) %>% 
    pivot_wider(names_from = Hail, values_from = n) %>% 
    mutate(n = str_glue("{`FALSE`+`TRUE`} ({`TRUE`})")) %>% 
    select(cloud_type, n)
}

nouns <- c("schaal", "staal", "stof", "spot", "blik", "horde", "hoop")
adjs <- c("hachelijk", "heet", "gemeen", "dof", "geldig", "heilzaam", "hemels", "grijs",
          "hoopvol", "goedkoop", "hoekig", "gekleurd", "geestig")
verbs <- c("heffen", "harden", "herstellen", "herroepen", "huldigen", "haken", "herhalen",
           "herstructureren", "herinneren", "diskwalificeren", "helpen", "haten")

by_lemma <- cloud_data %>% 
    filter(cluster != "0" | maincat == "Cirrostratus") %>% 
    select(lemma, model, cloud_type) %>% distinct() %>% 
    count(lemma, cloud_type) %>% 
    mutate(n = if_else(lemma %in% nouns, n/212, n/200)) %>% 
    pivot_wider(names_from = cloud_type, values_from = n, values_fill = 0) %>% 
    select(lemma, all_of(cloud_order)) %>% 
    arrange(desc(Cumulus))

mean_n_clusters <- cloud_data %>%
  filter(cluster != '0' | maincat == 'Cirrostratus') %>%
  count(lemma, model) %>% group_by(lemma) %>%
  summarize(mean_n_clusters = mean(n)) %>% arrange(mean_n_clusters) %>% 
  deframe %>% round(3)
