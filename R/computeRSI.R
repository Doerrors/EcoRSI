validate_EcoRSI_data <- function(data, site_col, index_col, umd_label, re_label){

  # 1. Check UMD & RE exist
  sites_present <- unique(na.omit(data[[site_col]]))

  if(!(umd_label %in% sites_present))
    stop("UMD site not found in data")

  if(!(re_label %in% sites_present))
    stop("RE site not found in data")

  # 2. Check missing site/index
  if(any(is.na(data[[site_col]])))
    stop("Missing values found in site column")

  if(any(is.na(data[[index_col]])))
    stop("Missing values found in index column")

  # 3. Equal rows per site
  site_counts <- table(data[[site_col]])
  if(length(unique(site_counts)) != 1){
    stop(paste("Unequal rows per site:",
               paste(names(site_counts), site_counts, collapse = ", ")))
  }

  # 4. Index repetition check
  expected_index <- sort(unique(data[[index_col]]))

  index_check <- tapply(data[[index_col]], data[[site_col]], function(x){
    identical(sort(unique(x)), expected_index)
  })

  if(!all(index_check)){
    stop("Index values are not repeating correctly across sites")
  }

  # 5. Duplicate index inside site
  dup_check <- tapply(data[[index_col]], data[[site_col]], function(x){
    any(duplicated(x))
  })

  if(any(dup_check)){
    stop("Duplicate index values found within a site")
  }

  TRUE
}


#' Compute Restoration State Index (RSI)
#'
#' Computes the Restoration State Index (RSI) as a composite,
#' reference ecosystem-based framework to quantify ecological
#' recovery trajectories. RSI integrates multidimensional indicators
#' across functional domains (e.g. vegetation, fauna, soil nutrients,
#' heavy metals) into a unified measure of restoration status relative
#' to a reference ecosystem.
#'
#' By default, equal weights are applied across all domains
#' (w_k = 1/K), ensuring cross-site and cross-study comparability.
#' Custom domain weights can be supplied when the ecological context
#' demands context-specific prioritization (e.g. soil-focused mine
#' restoration vs. biodiversity-focused forest restoration). When
#' custom weights are used, RSI scores are context-specific and
#' may not be directly comparable across studies using default
#' equal weights.
#'
#' @param data A data frame containing ecological data.
#' @param umd_label Character. Label identifying the degraded baseline
#'   site (Unrestored Mine Dump or equivalent) used to compute
#'   lower bounds (L) and upper bounds (U) for normalization.
#' @param re_label Character. Label identifying the Reference Ecosystem
#'   site. All indicators are normalized relative to this site.
#' @param domain_list A named list where each element is a character
#'   vector of column names belonging to that domain.
#'   Example:
#'   \code{list(
#'     Vegetation  = "floral_richness",
#'     Fauna       = "faunal_richness",
#'     SoilNut     = c("N", "P", "K", "OC"),
#'     HeavyMetal  = c("Mn", "Fe", "pH", "Pb", "Cr")
#'   )}
#' @param positive_indicators Character vector of indicator column names
#'   where higher values reflect improved ecosystem condition
#'   (e.g. species richness, soil organic carbon, nutrients).
#'   Normalized as: z+ = (X - L) / (RE - L)
#' @param negative_indicators Character vector of indicator column names
#'   where higher values reflect greater degradation
#'   (e.g. heavy metals, salinity).
#'   Normalized as: z- = (U - X) / (U - RE)
#' @param index_col Character. Name of the column containing the
#'   time or sampling index (must be identical across all sites).
#' @param site_col Character. Name of the column identifying sites.
#' @param domain_weights Optional named numeric vector of weights for
#'   each domain. Names must exactly match names in \code{domain_list}.
#'   Weights must sum to 1. Default is \code{NULL} (equal weights
#'   applied automatically: w_k = 1/K).
#'   Example for soil-focused mine restoration:
#'   \code{c(Vegetation = 0.10, Fauna = 0.10,
#'            SoilNut = 0.40, HeavyMetal = 0.40)}
#'
#' @return A data frame with the following columns:
#'   \item{Index}{Sampling index value}
#'   \item{Site}{Site label}
#'   \item{One column per domain}{Domain-level sub-index (I_k),
#'     named as the domain names supplied in \code{domain_list}}
#'   \item{RSI}{Final Restoration State Index (weighted mean of
#'     domain sub-indices). Bounded between 0 and 1. Values near 1
#'     indicate strong similarity to the reference ecosystem.}
#'
#' @details
#' \strong{Computation steps:}
#' \enumerate{
#'   \item Bounds (L, U) are computed from the degraded baseline site
#'     (umd_label) across all time steps.
#'   \item Each indicator is normalized per time step relative to the
#'     matched Reference Ecosystem value at the same index.
#'     Positive: z+ = clamp(0,1, (X - L) / (RE - L))
#'     Negative: z- = clamp(0,1, (U - X) / (U - RE))
#'   \item Domain sub-indices (I_k) are computed as the arithmetic
#'     mean of all normalized indicator scores within each domain.
#'   \item Final RSI = sum(w_k * I_k) across all K domains,
#'     where w_k = 1/K by default (equal weighting).
#' }
#'
#' \strong{Equal weighting rationale:}
#' Equal domain weights are the default because RSI measures holistic
#' ecosystem convergence toward a reference state. In a fully restored
#' ecosystem, all domains must simultaneously converge. Prioritizing
#' one domain would imply that partial recovery is sufficient, which
#' contradicts the ecological definition of restoration success.
#' Equal weighting also ensures RSI scores are scale-invariant and
#' comparable across sites and studies.
#'
#' \strong{Custom weighting:}
#' A warning is issued when custom weights are supplied, reminding the
#' user that resulting RSI scores are context-specific.
#'
#' @export
computeRSI <- function(data,
                       umd_label,
                       re_label,
                       domain_list,
                       positive_indicators,
                       negative_indicators,
                       index_col,
                       site_col,
                       domain_weights = NULL) {

  # -------------------------------------------------------
  # 1. Input Validation
  # -------------------------------------------------------

  if (!is.data.frame(data))
    stop("`data` must be a data frame.")

  if (!is.character(umd_label) || !is.character(re_label))
    stop("`umd_label` and `re_label` must be character values.")

  if (!is.list(domain_list) || is.null(names(domain_list)))
    stop("`domain_list` must be a named list of character vectors.")

  if (!is.character(positive_indicators) || !is.character(negative_indicators))
    stop("Indicators must be provided as character vectors.")

  if (!is.character(index_col) || !is.character(site_col))
    stop("`index_col` and `site_col` must be character values.")

  # All indicators pooled
  indicators <- c(positive_indicators, negative_indicators)

  # Check every indicator in domain_list is declared as positive or negative
  all_domain_indicators <- unlist(domain_list, use.names = FALSE)

  undeclared <- setdiff(all_domain_indicators, indicators)
  if (length(undeclared) > 0)
    stop(paste(
      "Indicators in domain_list not found in positive_indicators or negative_indicators:",
      paste(undeclared, collapse = ", ")
    ))

  unassigned <- setdiff(indicators, all_domain_indicators)
  if (length(unassigned) > 0)
    stop(paste(
      "Indicators declared in positive/negative_indicators but not assigned to any domain:",
      paste(unassigned, collapse = ", ")
    ))

  # Check required columns exist in data
  required_cols <- c(indicators, index_col, site_col)
  missing_cols  <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0)
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))

  # -------------------------------------------------------
  # Validate domain_weights if supplied
  # -------------------------------------------------------

  K <- length(domain_list)

  if (!is.null(domain_weights)) {

    if (!is.numeric(domain_weights))
      stop("`domain_weights` must be a named numeric vector.")

    if (is.null(names(domain_weights)))
      stop("`domain_weights` must be named. Names must match names in `domain_list`.")

    missing_domains <- setdiff(names(domain_list), names(domain_weights))
    extra_domains   <- setdiff(names(domain_weights), names(domain_list))

    if (length(missing_domains) > 0)
      stop(paste(
        "These domains are in domain_list but missing from domain_weights:",
        paste(missing_domains, collapse = ", ")
      ))

    if (length(extra_domains) > 0)
      stop(paste(
        "These domains are in domain_weights but not in domain_list:",
        paste(extra_domains, collapse = ", ")
      ))

    if (abs(sum(domain_weights) - 1) > 1e-6)
      stop(paste0(
        "`domain_weights` must sum to 1. Current sum = ",
        round(sum(domain_weights), 6)
      ))

    # Reorder weights to match domain_list order
    domain_weights <- domain_weights[names(domain_list)]

    warning(
      "Custom domain weights applied. RSI scores are context-specific and ",
      "may not be directly comparable across studies using default equal weights."
    )

  } else {
    # Default: equal weights
    domain_weights        <- rep(1 / K, K)
    names(domain_weights) <- names(domain_list)
  }

  # -------------------------------------------------------
  # Data Structure Validation
  # -------------------------------------------------------

  validate_EcoRSI_data(data, site_col, index_col, umd_label, re_label)

  # -------------------------------------------------------
  # 2. Separate degraded baseline (UMD) and Reference Ecosystem
  # -------------------------------------------------------

  umd       <- data[data[[site_col]] == umd_label, , drop = FALSE]
  reference <- data[data[[site_col]] == re_label,  , drop = FALSE]

  if (nrow(umd) == 0)       stop("No rows found for UMD label.")
  if (nrow(reference) == 0) stop("No rows found for RE label.")

  # -------------------------------------------------------
  # 3. Compute Bounds from degraded baseline (UMD)
  #    L = lower bound (min of UMD) — used for positive indicators
  #    U = upper bound (max of UMD) — used for negative indicators
  # -------------------------------------------------------

  lower_bounds <- sapply(indicators, function(col) min(umd[[col]], na.rm = TRUE))
  upper_bounds <- sapply(indicators, function(col) max(umd[[col]], na.rm = TRUE))

  # -------------------------------------------------------
  # 4. Normalize All Indicators
  #    Positive: z+ = clamp(0,1, (X - L) / (RE - L))
  #    Negative: z- = clamp(0,1, (U - X) / (U - RE))
  # -------------------------------------------------------

  norm_matrix           <- matrix(NA_real_, nrow = nrow(data), ncol = length(indicators))
  colnames(norm_matrix) <- indicators   # named by indicator for domain lookup below

  for (j in seq_along(indicators)) {

    col <- indicators[j]

    for (i in seq_len(nrow(data))) {

      X   <- data[[col]][i]
      idx <- data[[index_col]][i]

      RE_row <- reference[reference[[index_col]] == idx, , drop = FALSE]

      if (nrow(RE_row) == 0 || is.na(X)) next

      RE <- RE_row[[col]][1]
      if (is.na(RE)) next

      if (col %in% positive_indicators) {
        L     <- lower_bounds[col]
        value <- if ((RE - L) != 0) (X - L) / (RE - L) else 0
      } else {
        U     <- upper_bounds[col]
        value <- if ((U - RE) != 0) (U - X) / (U - RE) else 0
      }

      # Clamp to [0, 1]
      norm_matrix[i, j] <- max(0, min(1, value))
    }
  }

  # -------------------------------------------------------
  # 5. Compute Domain Sub-Indices (I_k)
  #    I_k = arithmetic mean of normalized scores within domain k
  #    This is the step that was missing in the original code.
  # -------------------------------------------------------

  domain_subindex_matrix           <- matrix(NA_real_, nrow = nrow(data), ncol = K)
  colnames(domain_subindex_matrix) <- names(domain_list)

  for (k in seq_along(domain_list)) {

    domain_indicators <- domain_list[[k]]

    # Columns in norm_matrix that belong to this domain
    domain_cols <- which(colnames(norm_matrix) %in% domain_indicators)

    if (length(domain_cols) == 0) next

    if (length(domain_cols) == 1) {
      sub_scores <- norm_matrix[, domain_cols, drop = TRUE]
    } else {
      sub_scores <- norm_matrix[, domain_cols, drop = FALSE]
    }

    # I_k = mean of normalized indicator scores within this domain
    domain_subindex_matrix[, k] <- if (is.matrix(sub_scores)) {
      ifelse(
        rowSums(!is.na(sub_scores)) == 0,
        NA_real_,
        rowMeans(sub_scores, na.rm = TRUE)
      )
    } else {
      sub_scores  # single indicator domain: sub-index = the score itself
    }
  }

  # -------------------------------------------------------
  # 6. Compute Final RSI
  #    RSI = sum(w_k * I_k)   across all K domains
  # -------------------------------------------------------

  RSI <- rep(NA_real_, nrow(data))

  for (i in seq_len(nrow(data))) {

    row_subindices <- domain_subindex_matrix[i, ]
    valid          <- !is.na(row_subindices)

    if (!any(valid)) next

    # Renormalize weights for domains with valid scores
    # (handles edge case where a domain has all-NA scores)
    w_valid <- domain_weights[valid]
    w_valid <- w_valid / sum(w_valid)

    RSI[i] <- sum(w_valid * row_subindices[valid])
  }

  # -------------------------------------------------------
  # 7. Assemble Output
  #    Returns: Index, Site, one column per domain sub-index, RSI
  # -------------------------------------------------------

  result <- data.frame(
    Index = data[[index_col]],
    Site  = data[[site_col]],
    domain_subindex_matrix,
    RSI   = RSI,
    check.names = FALSE
  )

  return(result)
}
