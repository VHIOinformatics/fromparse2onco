#' Parse and Combine Variant Data
#'
#' Import and Filter by Samples
#'
#' @param minimum_samples_mutated Minimum number of samples that must be mutated to be included in the filtered matrix.
#' @return Matrix filtered by number of samples mutated
#' @import readxl
#' @import dplyr
#' @import stringr
#' @import purrr
#' @examples
#' \dontrun{
#' filtered_matrix <- importandfilterbysamples(minimum_samples_mutated = 5)
#'}
#' @export

importandfilterbysamples <- function(minimum_samples_mutated){
  # Read the oncoplot matrix
  maf.matrix <- as.matrix(read.table("onco_matrix.txt", header = TRUE, sep = '\t', quote = ""))

  # Calculate mutation ratio
  mut_ratio <- apply(maf.matrix, 1, function(x) {y = (x != ''); sum(y)})

  # Filter the matrix by minimum samples mutated
  maf.matrix.filt <- maf.matrix[names(mut_ratio)[mut_ratio >= minimum_samples_mutated], ]

  # Replace specific strings in the filtered matrix
  maf.matrix.filt <- gsub('Frame_Shift', 'Frameshift', maf.matrix.filt)
  colnames(maf.matrix.filt) <- gsub("\\.","-", colnames(maf.matrix.filt))

  # Return the filtered matrix
  return(maf.matrix.filt)
}
