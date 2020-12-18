library(tidyverse)

create_example_data <- function(n_samples = 20, n_taxa = 10, seed = 1337) {
  set.seed(seed)
  data <- rnbinom(n_samples * n_taxa, 10, 0.95) %>% matrix(nrow = n_samples)
  
  rownames(data) <- n_samples %>%
    seq() %>%
    paste0("sample_", .)
  
  colnames(data) <- n_taxa %>%
    seq() %>%
    paste0("taxon_", .)
  
  return(data)
}

expect_graph <- function(res) {
  expect_equal(res$graph %>% igraph::gsize(), n_taxa_pairs)
  expect_equal(res$graph %>% igraph::gorder(), n_taxa)
}

n_samples <- 20
n_taxa <- 10
data <- create_example_data(n_samples = n_samples, n_taxa = n_taxa)
n_taxa_pairs <- n_taxa*(n_taxa-1)/2

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("correlate pearson", {
  res <- correlate_pearson(data)

  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "pearson")
  expect_graph(res)
})

test_that("correlate spearman", {
  res <- correlate_spearman(data)

  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "spearman")
  expect_graph(res)
})

test_that("correlate mb", {
  options(mc.cores = 5)
  res <- correlate_mb(data)
  
  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "mb")
  expect_graph(res)
})

test_that("correlate mb fail", {
  # prevent other methods to be set in this explicit mb function
  expect_error(correlate_mb(data, method = "glasso"))
})

test_that("correlate fastspar", {
  res <- correlate_fastspar(data, bootstraps = 20, iterations = 10)

  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "sparcc")
  expect_graph(res)
})

test_that("correlate sparcc", {
  res <- correlate_spiec_easi_sparcc(data, bootstraps = 20, iterations = 10)

  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "sparcc")
  expect_graph(res)
})

test_that("parse pearson", {
  res <- Hmisc::rcorr(data, type = "pearson") %>% as_coabundance(method = "pearson")

  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "pearson")
  expect_graph(res)
})

test_that("parse spearman", {
  res <- Hmisc::rcorr(data, type = "spearman") %>% as_coabundance(method = "spearman")

  expect_s3_class(res, "coabundance")
  expect_equal(res$method, "spearman")
  expect_graph(res)
})
