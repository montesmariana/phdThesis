library(tidyverse)
library(kableExtra)
library(colorblindr)
library(sna)
library(GGally)
library(gghalves)
library(ggpubr)
library(cowplot)
library(here)
library(knitr)
library(extrafont)
library(jsonlite)
library(readxl)
library(network)

# automatically create a bib database for R packages ----
write_bib(c(
  .packages(), 'bookdown', 'knitr', 'rmarkdown',
  'cluster', 'vegan', 'Rtsne', 'dbscan', 'umap',
  'semvar', 'entropy', 'shiny', 'colorblindr',
  'plotly'
), 'assets/bib/packages.bib')

loadfonts(device = "pdf")
# utils ----
hmean <- function(...){
  return(1/mean(1/c(...), na.rm = TRUE))
}
na2zero <- function(x) {
  x[is.na(x)] <- "-"
  x
}

pmi2positive <- function(x) {
  x[is.na(x)] <- 0
  x[x < 0] <- 0
  x
}

cloud_foto <- function(fname, extension = "jpg"){
  include_graphics(
    here("assets", "img", paste0(fname, ".", extension)),
    auto_pdf = T,
    dpi = 600)
}

html2md <- function(txt) {
  str_replace_all(txt, "<strong>([^<]+)</strong>", "*\\1*") %>%
    str_replace_all("<span class='target'>([^<]+)</span>", "**\\1**") %>% 
    str_remove_all("<sup>[^<]+</sup>")
}

sc <- function(txt){
  str_glue('<span style="font-variant:small-caps;">{txt}</span>')
}

latIt <- function(txt){
  if (is.na(txt)) "" else str_replace(r"((\textit{X}))", "X", txt)
}

appPvalue <- function(p) {
  if (p < 0.001) {
    return ("p-value < 0.001")
  } else if (p < 0.01) {
    return ("p-value < 0.01")
  } else if (p < 0.05) {
    return ("p-value < 0.05")
  } else {
    template <- r"(p-value $\approx$ X)"
    return(str_replace(template, "X", as.character(round(p, 3))))
  }
}

reportStat <- function(stat){
  estimate <- stat$estimate
  p.value <- stat$p.value
  measure <- str_replace(r"($\greek$)", "greek", names(estimate))
  str_glue("{measure} = {round(estimate, 2)}, {appPvalue(p.value)}")
}
scUpper <- function(x) {
  up <- str_extract(x, "[A-Z]+")
  if (is.na(up)) x else str_replace(x, up, sc(tolower(up)))
  }

nameModel <- function(mname){
  parts <- str_split(mname, pattern = "\\.")[[1]]
  if (length(parts) == 4) parts <- c("lemma", parts)
  foc <- str_remove(parts[[2]], "BOW") %>% str_remove("LEMMA") %>% 
    # str_replace("([A-Z]+)([^A-Z]+)", paste0(sc(tolower("\\1")), "\\2"))
    scUpper
  pmi <- str_replace(parts[[3]], "PPMI", sc("ppmi"))
  length <- sc(tolower(str_remove(parts[[4]], 'LENGTH')))
  socpos <- str_remove(parts[[5]], 'SOCPOS')
  paste(foc, pmi, paste0(length, socpos), sep="-")
}

# load data ----
d <- readRDS(here("assets", "data.rds"))
lrel <- read_json(here("assets", "lemmarel.json"))
definitions <- read_excel(here("assets", "definitions.xlsx")) %>%
  filter(lemma != 'spoor') %>% 
  select(lemma, sense = code,
         enl = example, een = example_translation,
         dnl = definition, den = definition_translation,
         rel_freq = my_dist) %>%
  arrange(lemma, sense) %>%
  mutate(
    sense = if_else(
      str_starts(den, "\\d"),
      str_extract(den, "\\d(\\.\\d)?"),
      str_remove(sense, paste0(lemma, "_"))
    ),
    dnl = str_remove(dnl, '\\d(\\.\\d)? '),
    den = str_remove(den, '\\d(\\.\\d)? '),
    enl = map_chr(enl, function(x) if (is.na(x)) "" else latIt(x)),
    een = map_chr(een, function(x) if (is.na(x)) "" else latIt(x)),
    Dutch = if_else(is.na(dnl), str_glue(""), str_glue("{dnl} {enl}")),
    English = str_glue("{den} {een}")
    ) %>%
  select(lemma, sense, Dutch, English) 

heilzaam <- read_tsv(here("assets", "heilzaam.cws.tsv"), col_types = c("sense_4" = "c")) %>% 
  rename("BOW" = conc_distance, "CW" = cw_type) %>% 
  mutate(path = str_replace_all(path, "->", r"( $\\rightarrow$ )"),
         path = str_remove(path, "#"))

lemmas <- read_tsv(here("assets", "lemmas.tsv"), col_types = cols()) %>% 
  select(pos, lemma = type, frequency, batches) %>% 
  filter(lemma %in% names(d)) %>% 
  mutate(
    pos = if_else(pos == "adj", "adjectives", paste0(pos, "s")),
    pos = fct_relevel(pos, c('nouns', 'adjectives', 'verbs')),
    senses = map(lemma, ~prop.table(table(d[[.x]]$senses$sense))),
    n_senses = map_dbl(lemma, ~nrow(count(d[[.x]]$senses, sense)))) %>% 
  arrange(pos, batches, frequency)
confidence <- read_tsv(here("assets", "confidences.tsv"), col_types = cols()) %>% 
  filter(! lemma %in% c("spoor", "herkennen")) 
kappas <- read_tsv(here("assets", "kappas.tsv"), col_types = cols()) %>% 
  filter(! lemma %in% c("spoor", "herkennen")) 
final_agreement <- read_tsv(here("assets", "final.agreement.tsv"), col_types = cols()) %>% 
  filter(! type %in% c("spoor", "herkennen")) %>% 
  mutate(
    final_action = fct_relevel(final_action, c("Same", "Other", "Remove")),
    majority_agreement = fct_relevel(majority_agreement, c("Full", "Majority", "No agreement"))
  )

geen <- readRDS(here("assets", "geen_data.rds"))
cloud_order <- c("Cumulus", "Stratocumulus", "Cirrus", "Cumulonimbus", "Cirrostratus")
qlvlcorp <- readRDS(here("assets", "qlvlnews.summary.rds"))
medoid_data <- read_tsv(here("assets", "medoid_data.tsv"), col_types = cols()) %>% 
  mutate(cloud_type = fct_relevel(cloud_type, cloud_order))
cloud_data <- read_tsv(here("assets", "classification.tsv"), col_types = cols()) %>% 
  mutate(cloud_type = fct_relevel(cloud_type, cloud_order))
examples <- read_tsv("assets/examples_translated.tsv", col_types = cols())


to_show <- c("heet", "stof", "dof", "huldigen", "haten", "hoop")
popular_medoid <- "BOWbound5lex.PPMIselection.LENGTHFOC.SOCPOSall"

# plots ----

complex_cor <- function(df,
                        x_coord, y_coord, color_var,
                        xlab = "", ylab = "", colorlab = "",
                        add_abline = FALSE){
  ncats <- nrow(count(df, {{ color_var }}))
  baseplot <- df %>% 
    ggplot(aes(x = {{ x_coord }}, y = {{ y_coord }}, color = {{ color_var }})) +
    geom_point(alpha = 0.5, size = 3) +
    theme_pubr(base_family = "Bookman Old Style") +
    (if (ncats > 3) scale_colour_viridis_d(guide = guide_legend(
      direction = "horizontal",
      title.position = "top"
    )) else scale_colour_grey(start = 0.8, end = 0, guide = guide_legend(
      direction = "horizontal",
      title.position = "top"
    ))) +
    labs(x = xlab, y = ylab, color = colorlab)
  if (add_abline) baseplot <- baseplot + geom_abline()
  
  xplot <- axis_canvas(baseplot, axis = "x") +
    geom_boxplot(data = df, aes(x = {{ x_coord }}, fill = {{ color_var }}), color = "black", alpha = 0.5) +
    if (ncats > 3) scale_fill_viridis_d() else scale_fill_grey(start = 0.8, end = 0)
  yplot <- axis_canvas(baseplot, axis = "y", coord_flip = TRUE) +
    geom_boxplot(data = df, aes(x = {{ y_coord }}, fill = {{ color_var }}), color = "black", alpha = 0.5) +
    (if (ncats > 3) scale_fill_viridis_d() else scale_fill_grey(start = 0.8, end = 0)) +
    coord_flip()
  
  baseplot %>% 
    insert_xaxis_grob(xplot, position = "top") %>%
    insert_yaxis_grob(yplot, position = "right") %>%
    ggdraw()
}

plotBasic <- function(mdata){
  m <- if ("lemma" %in% names(mdata)) d[[mdata$lemma]]$medoidCoords[[mdata$model]] else mdata
  m <- if ("coords" %in% names(m)) m else list(coords = m)
  m$coords <- mutate(m$coords, sense = "sense")
  plotCoords(m, colorguide = "none")
}

plotWithCws <- function(lname, mnum, colorguide = "right"){
  m <- d[[lname]]$medoidCoords[mnum]
  cw_data <- medoid_data %>% 
    filter(lemma == lname, model == names(m)) %>% 
    select(cluster, cloud_type, top_Cw, topFscore)
  plotCoords(m[[1]], cw_data, colorguide = colorguide)
}

plotCoords <- function(medoid, cw_data = tibble(), colorguide = "right") {
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
  plotCloud(df, cluster, sense, eps, colorguide = colorguide)
  
}

plotSenses <- function(l, mname){
  if (!str_starts(mname, l)) mname <- str_glue("{l}.{mname}")
  df <- d[[l]]$medoidCoords[[mname]]$coords %>% 
    mutate(sense = str_remove(sense, paste0(l, ".")))
  plotCloud(df, sense, "shape", 0, "bottom")
}
plotNude <- function(df) {
  plotCloud(df, colorguide = "none")
}
plotCloud <- function(df, colorvar = "", shapevar = "", alphavar = 1, colorguide = "right"){
  g <- ggplot(df) +
    geom_point(aes(x = model.x, y = model.y,
                   color = {{ colorvar }},
                   shape = {{ shapevar }},
                   alpha = {{ alphavar }}),
               size = 3) +
    theme_void() +
    coord_fixed() +
    scale_shape(guide = "none") +
    theme(legend.position = colorguide)
  g <- if (g$labels$alpha != "alpha") g + scale_alpha(range = c(1, 0), guide = "none") else g + scale_alpha(range = c(0.5, 0.5), guide = "none")
  g <- if (g$labels$colour != "colour") g + scale_color_OkabeIto(darken = 0.1, use_black = TRUE, na.value = "#9b9c9f") else g + scale_colour_grey()
  g
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
  net <- network(cws_mtx, ignore.eval = FALSE, names.eval = "weights")
  net %v% "freq" = rowSums(cws_mtx)
  cw_freq <- map_dbl(all_cws, ~max(cws_mtx[.x,]))
  net %v% "cw_freq" = case_when(cw_freq > 10 ~ 6, cw_freq > 1 ~ 4, TRUE ~2)
  set.seed(8541)
  ggnet2(net, size = "freq", label = TRUE, label.size ="cw_freq",
                 edge.size = "weights", edge.alpha = 0.1) +
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
readDutch <- function(ex) {
  filter(examples, code == ex) %>% 
    mutate(ctxt = str_glue("{dutch} (*{source}*, {date}, Art. {artn})")) %>% 
    pull(ctxt)
}
readTranslation <- function(ex) {
  en <- filter(examples, code == ex) %>% 
    pull(english)
  return(str_glue("'{en}'"))
}

readExample <- function(ex) {
  filter(examples, code == ex) %>% 
    mutate(ctxt = str_glue("{dutch}   <br><br>\\s\\s{english}   <br> \n <br>    (*{source}*, {date}, Art. {artn})")) %>% 
    pull(ctxt)
}
 
sampleCloud <- function(cloudtype) {
  set.seed(8541)
  cloud_selection <- medoid_data %>% 
    filter(
      (cloudtype %in% maincat & maincat == cloudtype) |
        (!cloudtype %in% maincat & str_detect(maincat, cloudtype)),
      !Hail
    )
  if (cloudtype == "Cumulonimbus") slice(arrange(cloud_selection, desc(rel_size)), 1) else slice_sample(cloud_selection, n = 1)
}
sampleHail <- function() {
  set.seed(8541)
  medoid_data %>%
    filter(maincat != "Cirrostratus") %>% 
    group_by(lemma, model) %>% 
    summarize(hail_p = sum(Hail)/n(), n = n()) %>% 
    ungroup() %>% 
    filter(hail_p == max(hail_p)) %>% 
    slice_sample(n = 1)
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
      numbers = str_replace(source, "[a-z_]+_([0-9_]+)", "\\1"),
      source = source_names[str_remove(source, paste0("_", numbers))],
      numbers = str_replace(numbers, "_01_", "_")
    ) %>% 
    separate(numbers, into = c("date", "artn"), remove = F) %>% 
    mutate(date = as.Date(date, format = "%Y%m%d"),
           txt = str_glue("{ctxt} (*{source}*, {date}, Art. {artn})")) %>% 
    pull(txt) %>% html2md()
}

getBest <- function(l) {
  medoid_data %>% filter(lemma == l, cloud_type != "Cirrostratus") %>% 
    mutate(model = str_remove(model, paste0(lemma, "."))) %>% 
    select(model, entropy) %>% group_by(model) %>% 
    mutate(mean_entropy = mean(entropy)) %>% 
    arrange(mean_entropy) %>% 
    select(model, mean_entropy) %>% 
    ungroup() %>% 
    slice(1) %>% 
    deframe()
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

# Tables ----
showDefs <- function(lemmas, caption) {
  subdefs <- filter(definitions, lemma %in% lemmas)
  kbl(select(subdefs, Dutch, sense, English), escape = F, booktabs = T, caption = caption,
        longtable = T, linesep = "\\addlinespace") %>% 
    kable_paper(font_size = 7, latex_options = c("repeat_header"), full_width = T,
                repeat_header_method = "replace") %>% 
    pack_rows(index = table(subdefs$lemma),
              latex_align = "c",
              # latex_gap_space = "0.5em",
              indent = F) %>% 
    column_spec(2, width = "1.5em") %>%
    column_spec(c(1, 3), width = "19em") %>% 
    row_spec(0, align = "c", bold = T, font_size = 9)
}

countCues <- function(colname, n = 10){
  heilzaam %>% count({{ colname }}, majority_sense) %>% 
    arrange(desc(n), desc(majority_sense)) %>% 
    pivot_wider(names_from = majority_sense, values_from = n, values_fill = 0) %>% 
    head(n)
}
