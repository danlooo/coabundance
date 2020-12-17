cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))


print.coabundance <- function(x) {
  cat_subtle(str_glue("# A coabundance object of method {x$method}: {x$graph %>% gsize()} interactions and {x$graph %>% gorder()} nodes \n\n"))
  cat_subtle("# Nodes:\n")
  print(x$graph %>% activate(nodes) %>% as_tibble(), n = 5)
  cat_subtle("# Interactions:\n")
  print(x$graph %>% activate(edges) %>% as_tibble() %>% arrange(-estimate), n = 5)
}
