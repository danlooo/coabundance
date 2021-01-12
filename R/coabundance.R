#' Create a coabundance object
#' @param cor_res result object originally describing the relationships (e.g. object of class `rcorr`)
#' @param edges tibble with columns from and to describing edges
#' @param method character of correlation method used
#' @export
coabundance <- function(cor_res, edges, nodes = NULL, method = NULL, max_pval = NULL, min_abs_estimate = NULL, ...) {
  if (!is.null(max_pval)) {
    if ("p.value" %in% colnames(edges)) {
      edges <-
        edges %>%
        filter(p.value <= max_pval)
    } else {
      warning("Ignore option max_pval: Not applicable")
    }
  }

  if (!is.null(min_abs_estimate)) {
    edges <-
      edges %>%
      filter(abs(estimate) >= min_abs_estimate)
  }

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

  res <- res %>% topologize()

  res
}

#' @export
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
      estimate = cor_res$pval$cors,
      p.value = cor_res$pval$pvals,
      q.value = p.adjust(p.value, method = "fdr")
    )

  coabundance(cor_res = cor_res, edges = edges, method = "sparcc", ...)
}

as_coabundance.rcorr <- function(cor_res, nodes = NULL, method = "rcorr", ...) {
  edges <- cor_res %>%
    broom::tidy() %>%
    dplyr::rename(from = column1, to = column2) %>%
    dplyr::mutate(q.value = p.adjust(p.value, method = "fdr"))

  coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method, ...)
}

as_coabundance.tbl_df <- function(cor_res, nodes = NULL, method = NULL, ...) {
  if (!all(c("from", "to", "estimate") %in% colnames(cor_res))) {
    stop("cor_res must have at least columns from and to")
  }

  coabundance(cor_res = cor_res, edges = cor_res, nodes = nodes, method = method, ...)
}

as_coabundance.pulsar.refit <- function(cor_res, nodes = NULL, method = NULL, ...) {
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
    tidygraph::activate(edges) %>%
    left_join(edges, by = c("from", "to")) %>%
    tidygraph::activate(nodes)

  cur_nodes <-
    graph %>%
    tidygraph::activate(nodes) %>%
    as_tibble()

  edges <-
    graph %>%
    tidygraph::activate(edges) %>%
    as_tibble() %>%
    left_join(cur_nodes %>% dplyr::rename(from_taxon = taxon), by = c("from" = "name")) %>%
    left_join(cur_nodes %>% dplyr::rename(to_taxon = taxon), by = c("to" = "name")) %>%
    dplyr::select(from = from_taxon, to = to_taxon, estimate)

  coabundance(cor_res = cor_res, edges = edges, nodes = nodes, method = method, ...)
}

as_coabundance.coabundance <- function(x, nodes = NULL, method = NULL, ...) {
  if (is.null(method)) {
    method <- x$method
  }

  if (is.null(nodes)) {
    return(x)
  }

  edges <-
    x$graph %>%
    tidygraph::activate(edges) %>%
    as_tibble()

  nodes <-
    x$graph %>%
    tidygraph::activate(nodes) %>%
    as_tibble() %>%
    rename(taxon = name) %>%
    left_join(nodes, by = "taxon")


  coabundance(cor_res = x$result, edges = edges, nodes = nodes, method = method, ...)
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
  res %>% as_coabundance.coabundance(...)
}

#' Convert coabundance objects
#' @export
as_coabundance <- function(x, ...) {
  UseMethod("as_coabundance")
}

cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))

#' @export
print.coabundance <- function(x) {
  n_interactions <- igraph::gsize(x$graph)
  n_nodes <- igraph::gorder(x$graph)

  cat_subtle(str_glue("# A coabundance object of method {x$method}: {n_interactions} interactions and {n_nodes} nodes \n\n"))
  cat_subtle("# Nodes:\n")
  print(x$graph %>% tidygraph::activate(nodes) %>% as_tibble(), n = 5)
  cat_subtle("# Interactions:\n")
  print(x$graph %>% tidygraph::activate(edges) %>% as_tibble() %>% arrange(-estimate), n = 5)
}

#' Adds topology data about nodes
#' @export
topologize <- function(x) {
  graph <- x$graph
  orig_state <- graph %>% tidygraph::active()

  graph <-
    graph %>%
    tidygraph::activate(nodes) %>%
    mutate(
      degree = tidygraph::centrality_degree(),
      component = tidygraph::group_components(),
      closeness = tidygraph::centrality_closeness(),
      betweeness = tidygraph::centrality_betweenness()
    )

  graph <- graph %>% tidygraph::activate(!!orig_state)
  x$graph <- graph
  x
}

#' Filter a coabundance object
#'
#' @param x coabundance object
#' @param max_pval maximum p value to keep an edge. Will be ignores if the p.value is not available
#' @param min_abs_estimate minimal absolute value of the estimate to keep an edge. This is to filter edges with a low effect size.
#' @param remove_isolated_nodes TRUE if nodes not being part of any edge after filtering should be removed, FALSE otherwise.
#' @param recalculate_topology TRUE if node topology metrics e.g. centrality scores should be recalculated after filtering, FALSE otherwise.
#' @export
filter.coabundance <- function(x, max_pval = 0.05, min_abs_estimate = NULL, remove_isolated_nodes = TRUE, recalculate_topology = TRUE) {
  graph <- x$graph
  orig_state <- graph %>% tidygraph::active()

  edge_colnames <- graph %>%
    tidygraph::activate(edges) %>%
    as_tibble() %>%
    colnames()

  if (!is.null(max_pval)) {
    if (!"p.value" %in% edge_colnames) {
      warning("Ignore option max_pval: Not applicable")
    } else {
      graph <-
        graph %>%
        tidygraph::activate(edges) %>%
        tidygraph::filter(p.value <= max_pval)
    }
  }

  if (!is.null(min_abs_estimate)) {
    graph <-
      graph %>%
      tidygraph::activate(edges) %>%
      filter(estimate >= min_abs_estimate)
  }

  if (remove_isolated_nodes) {
    graph <-
      graph %>%
      tidygraph::activate(nodes) %>%
      filter(!tidygraph::node_is_isolated())
  }

  graph <- graph %>% tidygraph::activate(!!orig_state)
  x$graph <- graph

  if (recalculate_topology) {
    x <-
      x %>%
      topologize()
  }

  x
}
