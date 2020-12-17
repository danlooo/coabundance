#' Do fastspar coabundance
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
correlate_fastspar <- function(data, iterations = 50, exclude_iterations = 10, bootstraps = 200, threads = getOption("mc.cores")) {
  fastspar_bin_dir_path <-"/analysis/miniconda3/bin"
  Sys.setenv(PATH = paste0(fastspar_bin_dir_path, ":", Sys.getenv("PATH")))
  
  system <- function(...) base::system(ignore.stdout = TRUE, ignore.stderr = TRUE, ...)
  
  # sanity checks
  if(class(data) != "matrix") stop("data must be of type matrix")
  if(str_glue("fastspar --version") %>% system() != 0) {
    stop("Command fastspar not found")
  }
  
  dir <- tempfile(pattern = "fastspar")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  
  data_path <- str_glue("{dir}/data.tsv")
  
  data %>%
    t() %>%
    as_tibble(rownames = "#OTU ID") %>%
    write_tsv(data_path)
  
  paste(
    "fastspar",
    "--yes",
    "--iterations", iterations,
    "--exclude_iterations", exclude_iterations,
    "--otu_table", data_path,
    "--correlation", paste0(dir, "/cor.tsv"),
    "--covariance",  paste0(dir, "median_covariance.tsv"),
    "--threads", threads,
    sep = " "
  ) %>%
    system()
  
  
  paste0(dir, "/bootstraps_counts") %>% dir.create()
  paste0(dir, "/bootstraps_cor") %>% dir.create()
  
  paste(
    "fastspar_bootstrap",
    "--otu_table", data_path,
    "--number", bootstraps,
    "--prefix", paste0(dir, "/bootstraps_counts", "/data"),
    sep = " "
  ) %>%
    system()
  
  paste(
    "parallel",
    "--jobs", threads,
    
    "fastspar",
    "--yes",
    "--otu_table {}",
    "--correlation", paste0(dir, "/bootstraps_cor/cor_{/}"),
    "--covariance", paste0(dir, "/bootstraps_cor/cov_{/}"),
    "--iterations", iterations,
    "--exclude_iterations", exclude_iterations,
    ":::",
    paste0(dir, "/bootstraps_counts/*"),
    sep = " "
  ) %>%
    system()
  
  paste(
    "fastspar_pvalues",
    "--otu_table", data_path,
    "--correlation", paste0(dir, "/cor.tsv"),
    "--prefix", paste0(dir, "/bootstraps_cor/cor_data_"),
    "--permutations", bootstraps,
    "--outfile", paste0(dir, "/pvals.tsv"),
    sep = " "
  ) %>%
    system()
  
  pval_tbl <-
    paste0(dir, "/pvals.tsv") %>%
    read_tsv(col_types = cols(.default = "c")) %>%
    dplyr::rename(from = `#OTU ID`) %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "p.value")
  
  cor_tbl <-
    paste0(dir, "/cor.tsv") %>%
    readr::read_tsv(col_types = cols(.default = "c")) %>%
    dplyr::rename(from = `#OTU ID`) %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "estimate")
  
  paste0("rm -rf ", dir) %>% system()
  
  res <-
    pval_tbl %>%
    dplyr::inner_join(cor_tbl, by = c("from", "to")) %>%
    # only keep triangle
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>% sort() %>% paste0(collapse = "-"))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::select(-comp) %>%
    readr::type_convert() %>%
    dplyr::ungroup() %>%
    mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(from, to, p.value, q.value, estimate)
  
  res %>% as_coabundance(method = "")
}
