#' Coabundance Analysis
#' 
#' @details 
#'  The option `method` (`sparcc` by default) defines the method use to calculate the interaction between pairs of coabundant taxa:
#'  
#'  * `pearson`: Pearson correlation coefficient
#'  * `spearman`: Spearman's rank correlation coefficient
#'  * `mb`: Inverse Covariance based on \insertCite{mb}{coabundance} as implemented in \insertCite{spiec_easi}{coabundance}
#'  * `sparcc`: SparCC correlation based on \insertCite{sparcc}{coabundance} as implemented in \insertCite{spiec_easi}{coabundance}
#'  or \insertCite{fastspar}{coabundance}
#'  
#' @references 
#' 
#' \insertRef{mb}{coabundance}
#' 
#' \insertRef{spiec_easi}{coabundance}
#' 
#' \insertRef{sparcc}{coabundance}
#' 
#' \insertRef{fastspar}{coabundance}
#' 
#' @param data matrix or data frame with abundance count data
#' @param method character of coabundance method. One of 'sparcc', 'mb', 'pearson', or 'spearman'
#' @param ... arguments passed to selected function (one of \link[coabundance]{correlate_sparcc},
#' \link[coabundance]{correlate_mb}, \link[coabundance]{correlate_spearman}, or \link[coabundance]{correlate_pearson})
#' @export
correlate <- function(data, method = "sparcc", ...) {
  switch(method,
         "sparcc" = correlate_sparcc(data = data, ...),
         "mb" = correlate_mb(data = data, ...),
         "pearson" = correlate_pearson(data = data, ...),
         "spearman" = correlate_spearman(data = data, ...),
         stop(str_glue("method {method} is not implemented!"))
    )
}

#' Coabundance Analysis Using SparCC as Implemented in Fastspar
#'
#' This implementation is way faster than `correlate_spiec_easi_sparcc` but requires linux and the external shell command `fastspar`.
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
#' @export
correlate_fastspar <- function(data, iterations = 50, exclude_iterations = 10, bootstraps = 200, threads = getOption("mc.cores")) {
  system <- function(...) base::system(ignore.stdout = TRUE, ignore.stderr = TRUE, ...)

  threads <- min(parallel::detectCores(), threads)

  # sanity checks
  if (class(data) != "matrix") stop("data must be of type matrix")
  if (str_glue("fastspar --version") %>% system() != 0) {
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
    "--covariance", paste0(dir, "median_covariance.tsv"),
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
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>%
      sort() %>%
      paste0(collapse = "-"))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::select(-comp) %>%
    readr::type_convert() %>%
    dplyr::ungroup() %>%
    dplyr::filter(from != to) %>%
    dplyr::mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(from, to, p.value, q.value, estimate)

  res %>% as_coabundance(method = "sparcc")
}

#' Coabundance analysis using mb as implemented in SpiecEasi
#'
#' This is a wrapper arround function `SpiecEasi::spiec.easi` with argument `method = "mb"`.
#'
#' @references
#'
#' \insertRef{spiec_easi}{coabundance}
#'
#' \insertRef{mb}{coabundance}
#'
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
#' @param pulsar.params list of options passed to \link[SpiecEasi]{pulsar.params}
#' @param ... further options passed to \link[SpiecEasi]{spiec.easi}
#' @export
correlate_mb <- function(
                         data,
                         pulsar.params = list(
                           thresh = 0.05,
                           subsample.ratio = 0.8,
                           ncores = getOption("mc.cores"),
                           rep.num = 20,
                           seed = 1337
                         ), ...) {
  params <- list(...)
  
  if(! is.null(params[["method"]])) {
    stop("Argument method as part of ... must not be set. This method uses mb explicitly.")
  }
  
  if(is.null(pulsar.params$ncores)) {
    warning("Option mc.cores is not set. Defaulting to one thread.")
    pulsar.params$ncores <- 1
  }
  
  list(
    data = data,
    method = "mb",
    pulsar.params = pulsar.params,
  ) %>%
    c(params) %>%
    do.call(SpiecEasi::spiec.easi) %>%
    as_cobundance(method = "mb")
}

#' Coabundance analysis using SparCC
#' @param implementation Character indicating the implementation of SparCC algorithm to use. One of "fastspar", "spiec_easi"
#' @param ... further arguments passed to the sparcc correlation functions
#' @export
correlate_sparcc <- function(implementation = "spiec_easi", ...) {
  switch(implementation,
    "fastspar" = correlate_fastspar(...),
    "spiec_easi" = correlate_spiec_easi_sparcc(...),
    stop(stringr::str_glue("{implementation} must be one of fastspar, auto, or spiec_easi"))
  )
}

#' Coabundance analysis using SparCC as implemented in SpiecEasi
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
#' @export
correlate_spiec_easi_sparcc <- function(data, iterations = 10, bootstraps = 200, th = 0.1, threads = getOption("mc.cores")) {
  sparcc_boot <-
    SpiecEasi::sparccboot(
      data = data,
      R = bootstraps,
      ncpus = threads,
      sparcc.params = list(
        iter = iterations,
        inner_iter = iterations,
        th = th
      )
    )

  sparcc_pval <-
    sparcc_boot %>%
    SpiecEasi::pval.sparccboot(sided = "both")

  list(boot = sparcc_boot, pval = sparcc_pval) %>%
    structure(class = "spiec_easi_sparcc_res") %>%
    as_coabundance(cor_res = .)
}


#' Correlate Pearson
#' @export
correlate_pearson <- function(data) {
  data %>%
    Hmisc::rcorr(type = "pearson") %>%
    as_coabundance(method = "pearson")
}


#' Correlate Spearman
#' @export
correlate_spearman <- function(data) {
  data %>%
    Hmisc::rcorr(type = "spearman") %>%
    as_coabundance(method = "spearman")
}