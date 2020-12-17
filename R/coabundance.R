#' Create a coabundance object
#' @param cor_res result object originally describing the relationships (e.g. object of class `rcorr`)
#' @param edges tibble with columns from and to describing edges
#' @param method character of correlation method used
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

as_coabundance.spiec_easi_sparcc_res <- function(cor_res, ...) {
  taxa <- cor_res$boot$data %>% colnames()

  edges <-
    tidyr::expand_grid(from = taxa, to = taxa) %>%
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>%
      sort() %>%
      paste0(collapse = ""))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::filter(from != to) %>%
    dplyr::ungroup() %>%
    dplyr::select(-comp) %>%
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
    getOptBeta() %>%
    as.matrix() %>%
    as_tibble(rownames = "from") %>%
    pivot_longer(-one_of("from"), names_to = "to", values_to = "estimate") %>%
    mutate(to = to %>% str_remove("^V")) %>%
    readr::type_convert(col_types = cols(from = col_integer(), to = col_integer(), estimate = col_double()))

  graph <-
    cor_res %>%
    getRefit() %>%
    adj2igraph() %>%
    as_tbl_graph() %>%
    mutate(taxon = used_taxa) %>%
    activate(edges) %>%
    left_join(edges, by = c("from", "to")) %>%
    activate(nodes)

  cur_nodes <-
    graph %>%
    activate(nodes) %>%
    as_tibble()

  edges <-
    graph %>%
    activate(edges) %>%
    as_tibble() %>%
    left_join(cur_nodes %>% dplyr::rename(from_taxon = taxon), by = c("from" = "name")) %>%
    left_join(cur_nodes %>% dplyr::rename(to_taxon = taxon), by = c("to" = "name")) %>%
    dplyr::select(from = from_taxon, to = to_taxon, estimate)

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
as_coabundance <- function(x, ...) {
  UseMethod("as_coabundance")
}
