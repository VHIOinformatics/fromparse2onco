#' Create oncomatrix from MAF file
#'
#' This function reads a MAF file, filters it and creates an oncomatrix. It also prints summary plots and creates a TMB table. It does so by executing functions in maftools package.
#'
#' @param maf_df A data frame in MAF format.
#' @param remove A logical value indicating whether to remove duplicated variants. Equivalent to removeDuplicatedVariants in maftools::read.maf function. Default is FALSE.
#' @param flags A logical value indicating whether to remove FLAG genes. Equivalent to rmFlags in maftools::read.maf function. Default is FALSE.
#' @param minimalMutations Minimum number of samples mutated in a gene to be drawn. Equivalent to minMut in maftools::oncoplot function. Deafult is 2.
#' @param topgenes Number of top genes to be drawn. Equivalent to top in maftools::oncoplot function. Default = 20000.
#' @param nonSyn Vector of Variant Classifications to keep. Equivalent to vc_nonSyn in maftools::read.maf function. Default is c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation").
#' 
#' @return Summary plots. It also saves the oncomatrix (onco_matrix.txt) and the TMB table (TMB_table.txt) in the working directory.
#' 
#' @import dplyr
#' @import maftools
#' @import stringr
#' 
#' @examples
#' # Example usage:
#' # modclassvartmb(filtered_df, remove=TRUE, minimalMutations=3)
#'
#' @export
modclassvartmb <- function(maf_df, remove=FALSE, flags=FALSE, minimalMutations = 2, topgenes = 20000, nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")) {
  # Read and summarize MAF file
  maf_object <- maftools::read.maf(maf = maf_df, removeDuplicatedVariants = remove, rmFlags=flags,vc_nonSyn = nonSyn)
  maftools::plotmafSummary(maf = maf_object,
                           addStat = 'median',
                           titvRaw = FALSE,
                           showBarcodes = TRUE)

  # Should we print the number of unique genes in the input matrix?
  #print(length(unique(maf_df$Hugo_Symbol)))

  # Create an oncoplot
  maftools::oncoplot(maf = maf_object,
                     minMut = minimalMutations,
                     showTumorSampleBarcodes = TRUE,
                     top = topgenes,
                     removeNonMutated = TRUE,
                     writeMatrix = TRUE)

  tmb <- maftools::tmb(maf_object)
  tmb_df <- tmb[,c("Tumor_Sample_Barcode", "total_perMB")]
  tmb_df <- tmb_df[order(tmb_df$total_perMB), ]
  write.table(tmb_df, file = "TMB_table.txt", row.names = FALSE, sep = "\t")
  #Return maf_object or nothing?
}




