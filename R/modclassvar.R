#' Modify and Classify Variants, Calculate TMB
#'
#' This function modifies the classification of variants and generates summary plots.
#'
#' @param to_class A data frame containing the variants to be classified.
#' @param remove A logical value indicating whether to remove duplicated variants. Default is FALSE.
#' @param minimum_samples_mutated An integer specifying the minimum number of samples mutated to filter the matrix. Default is 1.
#'
#' @return A filtered matrix of the variant classifications.
#' @import dplyr
#' @import maftools
#' @import stringr
#' @examples
#' # Example usage:
#' # filtered_matrix <- modclassvar(variant_data, remove=TRUE, minimum_samples_mutated=5)
#'
#' @export
modclassvartmb <- function(to_class, remove=FALSE, minimum_samples_mutated=1,flags=FALSE, minimalMutations = 2, topgenes = 20000, nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")) {
  # Read and summarize MAF file
  maf_object <- maftools::read.maf(maf = to_class, removeDuplicatedVariants = remove, rmFlags=flags,vc_nonSyn = nonSyn)
  maftools::plotmafSummary(maf = maf_object,
                           addStat = 'median',
                           titvRaw = FALSE,
                           showBarcodes = TRUE)

  # Print the number of unique genes
  print(length(unique(to_class$Hugo_Symbol)))

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
  return(tmb_df)
  # ¿Deberíamos poner el tmb como return o el el objeto maf?
}




