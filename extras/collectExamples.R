sampleCtxt2 <- function(l, mnum = 1, clusn = NULL, stag = NULL, cw = NULL, seed = 8541, n = 1) {
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
           ctxt = html2md(ctxt)) %>% 
    select(source, date, artn, ctxt)
}
collect_examples <- function(file){
  examples <- grep("sampleCtxt", readLines(file), value = TRUE)
  if (length(examples) == 0) return(tibble())
  transformed <- str_extract(examples, "sampleCtxt[^`]+") %>% 
    str_replace("sampleCtxt", "sampleCtxt2")
  map_dfr(transformed, ~eval(parse(text = .x))) %>% 
    mutate(code = str_extract(examples, "@[^)]+")) %>% 
    select(code, everything())
}
bind_rows(
  collect_examples("semantic_interpretation.Rmd"),
  collect_examples("no_optimal_solution.Rmd")) %>% 
  write_tsv("assets/examples.tsv")

# Read translations

