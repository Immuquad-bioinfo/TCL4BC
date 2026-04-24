suppressMessages(library(quadrat))
suppressMessages(library(tidyverse))
args <- commandArgs(T)

txt_path <- args[1]
to_path <- args[2]
cal_divindx <- function(data = NULL, sfopath = NULL, MARGIN = 1, base = exp(1)){

  ## FUNCTION---------------------------------------------------------------------------------
  easy_divindx <- function(x, index = "shannon", MARGIN = 1, base = exp(1)) {
    x <- drop(as.matrix(x))
    if (any(x < 0, na.rm = TRUE))
      stop("input data must be non-negative")
    INDICES <- c("shannon", "simpson", "invsimpson")
    index <- match.arg(index, INDICES)
    if (length(dim(x)) > 1) {
      total <- apply(x, MARGIN, sum)
      x <- sweep(x, MARGIN, total, "/")
    }
    else {
      x <- x/(total <- sum(x))
    }
    if (index == "shannon")
      x <- -x * log(x, base)
    else x <- x * x
    if (length(dim(x)) > 1)
      H <- apply(x, MARGIN, sum, na.rm = TRUE)
    else H <- sum(x, na.rm = TRUE)
    if (index == "simpson")
      H <- 1 - H
    else if (index == "invsimpson")
      H <- 1/H
    if (any(is.na(total)))
      is.na(H) <- is.na(total)
    H
  }


  ## RUN--------------------------------------------------------------------------------------
  finame <- names(data)
  file_out <- file.path(sfopath, "stats_divIndx.csv")
  indices_matrix <- c()
  for (var in seq_along(data)) {
    dat <- data[[var]]
    # 1.singleton index
    singleton_per <- sum(dat$count[which(dat$count == 1)]) / sum(dat$count)
    # 2.Shannon index
    sha_ind <- easy_divindx(dat$count, index = "shannon", MARGIN = MARGIN, base = base)
    # 3.Pielou's index
    piel_ind <- sha_ind / log(length(dat$count), base = base)
    # 4. Clonality index
    clona <- 1 - piel_ind
    # 5.Simpson index
    sim_ind <- easy_divindx(dat$count, index = "simpson", MARGIN = MARGIN, base = base)
    # 6.invert-Simson index
    invsim_ind <- easy_divindx(dat$count, index = "invsimpson", MARGIN = MARGIN, base = base)
    # 7.hvj Index
    dat <- dat %>%
      dplyr::mutate(
        vj_comb = stringi::stri_c(vseg, "_", jseg)
      )
    dat_vj <- aggregate(dat$count, by = list(dat$vj_comb), sum)
    colnames(dat_vj) <- c("vj_comb", "count")
    hvj_ind <- easy_divindx(dat_vj[, "count"])
    # 8.TCR convergence
    st_cdr3aa_vseg <- dat %>%
      mutate(cdr3aa_vseg = paste0(cdr3aa, "_", vseg),
             cdr3aa_vseg_count = 1) %>%
      select(cdr3aa_vseg, cdr3aa_vseg_count) %>%
      aggregate(cdr3aa_vseg_count~., ., sum) %>%
      arrange(desc(cdr3aa_vseg_count)) %>%
      filter(cdr3aa_vseg_count > 1)
    convergence <- dat %>%
      mutate(cdr3aa_vseg = paste0(cdr3aa, "_", vseg)) %>%
      left_join(., st_cdr3aa_vseg, c("cdr3aa_vseg" = "cdr3aa_vseg")) %>%
      filter(is.na(cdr3aa_vseg_count) == F) %>%
      pull(freq) %>%
      sum()
    # Stats
    temp <- tibble(
      "FileName"         = finame[var],
      "Clonality"        = clona,
      "Pielous"          = piel_ind,
      "Shannon Index"    = sha_ind,
      "Invsimpson Index" = invsim_ind,
      "Simpson Index"    = sim_ind,
      "Hvj Index"        = hvj_ind,
      "Singleton"        = singleton_per,
      "Convergence"      = convergence
    )
    indices_matrix <- rbind(indices_matrix, temp)
  }
  write.csv(data.frame(indices_matrix), file = file_out, quote = FALSE, row.names = FALSE)
  return(indices_matrix)

}
dfs <- read_data(txt_path)
cal_divindx(dfs, sfopath = to_path)


