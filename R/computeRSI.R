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

#' Compute Relative Sustainability Index (RSI)
#'
#' Computes the Relative Sustainability Index (RSI) using
#' positive and negative ecological indicators based on
#' UMD and RE site references.
#'
#' @param data A data frame containing ecological data.
#' @param umd_label Character value identifying the UMD site.
#' @param re_label Character value identifying the RE site.
#' @param positive_indicators Character vector of positive indicator column names.
#' @param negative_indicators Character vector of negative indicator column names.
#' @param index_col Character name of index column.
#' @param site_col Character name of site column.
#'
#' @return A data frame containing Index, Site, and RSI values.
#' @export
computeRSI <- function(data,
                       umd_label,
                       re_label,
                       positive_indicators,
                       negative_indicators,
                       index_col,
                       site_col) {
  
  # -----------------------------
  # 1. Input Validation
  # -----------------------------
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  
  if (!is.character(umd_label) || !is.character(re_label)) {
    stop("`umd_label` and `re_label` must be character values.")
  }
  
  if (!is.character(positive_indicators) || !is.character(negative_indicators)) {
    stop("Indicators must be provided as character vectors.")
  }
  
  if (!is.character(index_col) || !is.character(site_col)) {
    stop("`index_col` and `site_col` must be character values.")
  }
  
  indicators <- c(positive_indicators, negative_indicators)
  
  required_cols <- c(indicators, index_col, site_col)
  
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  # -----------------------------
  # Data Structure Validation
  # -----------------------------
  validate_EcoRSI_data(
    data,
    site_col,
    index_col,
    umd_label,
    re_label
  )
  
  
  # -----------------------------
  # 2. Separate UMD and RE
  # -----------------------------
  
  umd <- data[data[[site_col]] == umd_label, , drop = FALSE]
  reference <- data[data[[site_col]] == re_label, , drop = FALSE]
  
  if (nrow(umd) == 0) {
    stop("No rows found for UMD label.")
  }
  
  if (nrow(reference) == 0) {
    stop("No rows found for RE label.")
  }
  
  # -----------------------------
  # 3. Compute Bounds from UMD
  # -----------------------------
  
  lower_bounds <- sapply(indicators, function(col) {
    min(umd[[col]], na.rm = TRUE)
  })
  
  upper_bounds <- sapply(indicators, function(col) {
    max(umd[[col]], na.rm = TRUE)
  })
  
  # -----------------------------
  # 4. Normalize Indicators
  # -----------------------------
  
  norm_matrix <- matrix(NA_real_, nrow = nrow(data), ncol = length(indicators))
  colnames(norm_matrix) <- paste0(indicators, "_norm")
  
  for (j in seq_along(indicators)) {
    
    col <- indicators[j]
    
    for (i in seq_len(nrow(data))) {
      
      X <- data[[col]][i]
      idx <- data[[index_col]][i]
      
      RE_row <- reference[reference[[index_col]] == idx, , drop = FALSE]
      
      if (nrow(RE_row) == 0 || is.na(X)) {
        next
      }
      
      RE <- RE_row[[col]][1]
      
      if (is.na(RE)) {
        next
      }
      
      if (col %in% positive_indicators) {
        
        L <- lower_bounds[col]
        
        value <- if ((RE - L) != 0) {
          (X - L) / (RE - L)
        } else {
          0
        }
        
      } else {
        
        U <- upper_bounds[col]
        
        value <- if ((U - RE) != 0) {
          (U - X) / (U - RE)
        } else {
          0
        }
      }
      
      # Clamp between 0 and 1
      value <- max(0, min(1, value))
      
      norm_matrix[i, j] <- value
    }
  }
  
  # -----------------------------
  # 5. Compute Final RSI
  # -----------------------------
  
  RSI <- ifelse(
    rowSums(!is.na(norm_matrix)) == 0,
    NA,
    rowMeans(norm_matrix, na.rm = TRUE)
  )
  
  result <- data.frame(
    Index = data[[index_col]],
    Site = data[[site_col]],
    RSI = RSI
  )
  
  return(result)
}