#' Create oncomatrix from MAF Data Frame
#'
#' This function takes a MAF dataframe, filters it and creates an oncomatrix. It also prints summary plots and creates a TMB table. It does so by executing functions in maftools package.
#'
#' @param maf_df A data frame in MAF format.
#' @param remove A logical value indicating whether to remove duplicated variants. Equivalent to removeDuplicatedVariants in maftools::read.maf function. Default is TRUE.
#' @param flags A logical value indicating whether to remove FLAG genes. Equivalent to rmFlags in maftools::read.maf function. Default is FALSE.
#' @param minimalMutations Minimum number of samples mutated in a gene to be drawn. Equivalent to minMut in maftools::oncoplot function. Default is 2.
#' @param topgenes Number of top genes to be drawn. Equivalent to top in maftools::oncoplot function. Default = 20000.
#' @param nonSyn Vector of Variant Classifications to keep. Equivalent to vc_nonSyn in maftools::read.maf function. Default is c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation").
#' 
#' @param captureSize The size (in Mb) of the sequenced region used to compute the TMB (TMB = mutations / captureSize). Default is 50 (whole exome, as in maftools::tmb). Use the real size of your panel if the data come from a targeted panel.
#' @param draw_oncoplot A logical value. maftools::oncoplot is only used here to write the oncomatrix (onco_matrix.txt); if FALSE the plot is drawn on a null device (nothing is shown) which saves time when topgenes is large. Default is TRUE.
#'
#' @return It prints summary plots and saves the oncomatrix (onco_matrix.txt) and the TMB table (TMB_table.txt) in the working directory. Invisibly, it returns a list with the maf object (`maf`) and the TMB table (`tmb`).
#' 
#' @import dplyr
#' @import maftools
#' @import stringr
#' @importFrom utils write.table
#' 
#' @examples
#' prepareForOncoplot(filtered_df, minimalMutations=3)
#'
#' @export
prepareForOncoplot <- function(maf_df, remove=TRUE, flags=FALSE, minimalMutations = 2, topgenes = 20000, nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"), captureSize = 50, draw_oncoplot = TRUE) {
  # Read and summarize MAF file
  maf_object <- maftools::read.maf(maf = maf_df, removeDuplicatedVariants = remove, rmFlags = flags, vc_nonSyn = nonSyn)
  maftools::plotmafSummary(maf = maf_object,
                           addStat = 'median',
                           titvRaw = FALSE,
                           showBarcodes = TRUE)

  # oncoplot() writes onco_matrix.txt (writeMatrix = TRUE) in the working directory.
  # Drawing thousands of genes is slow, so the plot can be sent to a null device.
  write_matrix <- function() {
    maftools::oncoplot(maf = maf_object,
                       minMut = minimalMutations,
                       showTumorSampleBarcodes = TRUE,
                       top = topgenes,
                       removeNonMutated = TRUE,
                       writeMatrix = TRUE)
  }
  if (draw_oncoplot) {
    write_matrix()
  } else {
    grDevices::pdf(NULL)
    tryCatch(write_matrix(), finally = grDevices::dev.off())
  }

  tmb <- maftools::tmb(maf_object, captureSize = captureSize)
  tmb_df <- tmb[, c("Tumor_Sample_Barcode", "total_perMB")]
  tmb_df <- tmb_df[order(tmb_df$total_perMB), ]
  write.table(tmb_df, file = "TMB_table.txt", row.names = FALSE, sep = "\t")

  invisible(list(maf = maf_object, tmb = tmb_df))
}
