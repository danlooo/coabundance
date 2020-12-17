cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))


print.coabundance <- function(x) {
  n_interactions <- igraph::gsize(x$graph)
  n_nodes <- igraph::gorder(x$graph)

  cat_subtle(str_glue("# A coabundance object of method {x$method}: {n_interactions} interactions and {n_nodes} nodes \n\n"))
  cat_subtle("# Nodes:\n")
  print(x$graph %>% activate(nodes) %>% as_tibble(), n = 5)
  cat_subtle("# Interactions:\n")
  print(x$graph %>% activate(edges) %>% as_tibble() %>% arrange(-estimate), n = 5)
}
