#' Create a coabundance object
#' @param cor_res result object originally describing the relationships (e.g. object of class `rcorr`)
#' @param edges tibble with columns from and to describing edges
#' @param method character of correlation method used
#' @export
coabundance <- function(cor_res, edges, nodes = NULL, method = NULL, ...) {
  if (!is.null(nodes)) {
    taxa <- edges$from %>%
      union(edges$to) %>%
      unique()
    nodes <- nodes %>% dplyr::filter(taxon %in% taxa)
  }

  edges <- edges %>% arrange(from, to)
  graph <- tidygraph::tbl_graph(edges = edges, nodes = nodes, directed = FALSE)

  res <- list(graph = graph, result = cor_res, method = method)
  class(res) <- "coabundance"
  res
}

cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))


print.coabundance <- function(x) {
  n_interactions <- igraph::gsize(x$graph)
  n_nodes <- igraph::gorder(x$graph)
  
  cat_subtle(str_glue("A coabundance object:\n   method: {x$method}\n   interactions: {n_interactions}\n   nodes: {n_nodes}\n"))
}


expand_nodes <- function(nodes) {
  tidyr::expand_grid(from = nodes, to = nodes) %>%
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>%
                                             sort() %>%
                                             paste0(collapse = ""))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::filter(from != to) %>%
    dplyr::ungroup() %>%
    dplyr::select(-comp) 
}

as_coabundance.spiec_easi_sparcc_res <- function(cor_res, ...) {
  taxa <- cor_res$boot$data %>% colnames()

  edges <-
    expand_nodes(taxa) %>%
    dplyr::mutate(
      estimate = sparcc_pval$cors,
      p.value = sparcc_pval$pvals,
      q.value = p.adjust(p.value, method = "fdr")
    )

  coabundance(cor_res = cor_res, edges = edges, method = "sparcc", ...)
}

as_coabundance.rcorr <- function(cor_res, nodes = NULL, method = "rcorr") {
  edges <- cor_res %>%
    broom::tidy() %>%
    dplyr::rename(from = column1, to = column2) %>%
    dplyr::mutate(q.value = p.adjust(p.value, method = "fdr"))

  coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method)
}

as_coabundance.tbl_df <- function(cor_res, nodes = NULL, method = NULL) {
  if (!all(c("from", "to", "estimate") %in% colnames(cor_res))) {
    stop("cor_res must have at least columns from and to")
  }

  coabundance(cor_res = cor_res, edges = cor_res, nodes = nodes, method = method)
}

as_coabundance.pulsar.refit <- function(cor_res, nodes = NULL, method = NULL) {
  if (method != "mb") stop("Only method mb implemented in objects of type pulsar.refit")

  used_taxa <- cor_res$est$data %>% colnames()

  edges <-
    cor_res %>%
    SpiecEasi::getOptCov() %>%
    as.matrix() %>%
    tibble::as_tibble(rownames = "from") %>%
    tidyr::pivot_longer(-one_of("from"), names_to = "to", values_to = "estimate") %>%
    dplyr::mutate(to = to %>% str_remove("^V")) %>%
    readr::type_convert(col_types = cols(from = col_integer(), to = col_integer(), estimate = col_double()))

  graph <-
    cor_res %>%
    SpiecEasi::getRefit() %>%
    SpiecEasi::adj2igraph() %>%
    tidygraph::as_tbl_graph() %>%
    tidygraph::mutate(taxon = used_taxa) %>%
    tidygraph::activate(edges) %>%
    dplyr::left_join(edges, by = c("from", "to")) %>%
    tidygraph::activate(nodes)

  cur_nodes <-
    graph %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble()
  
  zero_edges <-
    expand_nodes(used_taxa) %>%
    mutate(estimate = 0)

  edges <-
    graph %>%
    tidygraph::activate(edges) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(cur_nodes %>% dplyr::rename(from_taxon = taxon), by = c("from" = "name")) %>%
    dplyr::left_join(cur_nodes %>% dplyr::rename(to_taxon = taxon), by = c("to" = "name")) %>%
    dplyr::select(from = from_taxon, to = to_taxon, estimate) %>%
    bind_rows(zero_edges) %>%
    group_by(from, to) %>%
    arrange(- abs(estimate)) %>%
    slice(1) %>%
    ungroup()

  coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method)
}

as_coabundance.coabundance <- function(x) {
  x
}

as_coabundance.default <- function(cor_res, edges, nodes = NULL, method = NULL, ...) {
  if (!is.null(nodes)) {
    taxa <- edges$from %>%
      union(edges$to) %>%
      unique()
    nodes <- nodes %>% dplyr::filter(taxon %in% taxa)
  }

  edges <- edges %>% arrange(from, to)
  graph <- tidygraph::tbl_graph(edges = edges, nodes = nodes, directed = FALSE)

  res <- list(graph = graph, result = cor_res, method = method)
  class(res) <- "coabundance"
  res
}

#' Convert coabundance objects
#' 
#' @param x object to transform to class coabundance. Supported classes: `pulsar.refit`, `tbl_df`, `rcorr`, `spiec_easi_sparcc_res`
#' @export
as_coabundance <- function(x, ...) {
  UseMethod("as_coabundance")
}
