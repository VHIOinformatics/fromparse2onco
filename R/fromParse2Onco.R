#' Convert output of parse_variants program into an oncoplot
#'
#' This function reads one or more Excel files produced by parse_variants program and makes and oncoplot.
#'
#' @param path_to_parse A vector of paths to the excel files to be read.
#' @param tumor_only A logical value indicating if variant calling was performed in tumor only mode. Default is FALSE, indicating it was performed in paired mode.
#' 
#' @param filter_column A character vector specifying the values to filter the `FILTER` column by. Default is "PASS".
#' @param VAF_tumor The minimum variant allele frequency in the tumor sample (from 0 to 1). Default is 0.
#' @param VAF_control The maximum variant allele frequency in the control sample (from 0 to 1). Default is 0.
#' @param total_tumor_reads The minimum number of reads in the tumor sample. Default is 0.
#' @param alt_tumor_reads The minimum number of alternative reads in the tumor sample. Default is 0.
#' @param cgi A logical value indicating if the CGI Summary column should be used to filter keeping labels included in cgi_list.
#' @param oncokb A logical value indicating if the OncoKB column should be used to filter keeping labels included in oncokb_list.
#' @param cgi_list Vector of CGI Summary labels to keep.
#' @param oncokb_list Vector of OncoKB labels to keep.
#' @param annott Vector of Annotation Impact labels from SnpEff to keep.
#' @param tumor_samples_out Vector of tumor samples names to remove.
#' @param control_samples_out Vector of control samples names to remove (this makes sense when the same tumor sample has been analysed using different control samples as germline control).
#' 
#' @param remove A logical value indicating whether to remove duplicated variants. Equivalent to removeDuplicatedVariants in maftools::read.maf function. Default is TRUE.
#' @param flags A logical value indicating whether to remove FLAG genes. Equivalent to rmFlags in maftools::read.maf function. Default is FALSE.
#' @param minimalMutations Minimum number of samples mutated in a gene to be drawn. Equivalent to minMut in maftools::oncoplot function. Deafult is 2.
#' @param topgenes Number of top genes to be drawn. Equivalent to top in maftools::oncoplot function. Default = 20000.
#' @param nonSyn Vector of Variant Classifications to keep. Equivalent to vc_nonSyn in maftools::read.maf function. Default is c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation").
#'
#' @param Missense_color Color for missense variants. Default is "#2a9134".
#' @param Nonsense_color Color for nonsense variants. Default is "#ffca3a".
#' @param Nonstop_color Color for nonstop variants. Default is "#000000".
#' @param FrameDel_color Color for frame shift deletions. Default is "blue".
#' @param FrameIns_color Color for frame shift insertions. Default is "purple".
#' @param In_Frame_Ins_color Color for in-frame insertions. Default is "lightblue".
#' @param In_Frame_Del_color Color for in-frame deletions. Default is "plum1".
#' @param Translation_Start_Site_color Color for translation start sites. Default is "#ff0a54".
#' @param Splice_site_color Color for splicing sites. Default is "darkorange".
#' @param Multihit_color Color for multi-hit genes. Default is "#dab49d".
#' 
#' @return An oncoplot object which is also saved as a png file.
#'
#' @import readxl
#' @import dplyr
#' @import stringr
#' @import purrr
#' @import maftools
#' @import grid
#' @import ComplexHeatmap
#'
#' @examples
#' oncobuddy <-fromParse2Onco(path_to_parse = c("/path/to/file1","/path/to/file2"), VAF_tumor=0.05, minimalMutations=3, Missense_color="#FF5733")
#' 
#' @export

fromParse2Onco <- function(path_to_parse,tumor_only = FALSE,
                           filter_column=c("PASS"), VAF_tumor=0, VAF_control=0, total_tumor_reads=0, alt_tumor_reads=0, cgi=FALSE, oncokb=FALSE, cgi_list=c("oncogenic (predicted)", "oncogenic (predicted and annotated)", "oncogenic (annotated)"), oncokb_list=c("Likely Oncogenic", "Oncogenic"), annott=c("HIGH","MODERATE","MODIFIER"), tumor_samples_out = NULL, control_samples_out = NULL,
                           remove=TRUE, flags=FALSE, minimalMutations = 2, topgenes = 20000, nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),
                           Missense_color="#2a9134", Nonsense_color="#ffca3a", Nonstop_color="#000000", FrameDel_color="blue", FrameIns_color="purple", In_Frame_Ins_color="lightblue", In_Frame_Del_color="plum1", Translation_Start_Site_color="#ff0a54", Splice_site_color="darkorange", Multihit_color="#dab49d")
                           {
  
  # Read parse_variants output and convert to MAF
  variants_df <- fromParse2MAF(path_to_parse = path_to_parse, tumor_only, oncokb, cgi)
  
  # Filter MAF data frame
  filtered_df <- filterMAF(variants_df, tumor_only, filter_column, VAF_tumor, VAF_control, total_tumor_reads, alt_tumor_reads, cgi, oncokb, cgi_list, oncokb_list, annott, tumor_samples_out, control_samples_out)
  
  # Create oncomatrix from MAF data frame
  prepareForOncoplot(filtered_df, remove, flags, minimalMutations, topgenes, nonSyn)
  
  # Make oncoplot
  oncoplot <- makeOncoplot(Missense_color, Nonsense_color, Nonstop_color, FrameDel_color, FrameIns_color, In_Frame_Ins_color, In_Frame_Del_color, Translation_Start_Site_color, Splice_site_color, Multihit_color)
  
  return(oncoplot)
}

