#' Convert output of parseVCF program into an oncoplot
#'
#' This function reads one or more Excel files produced by parseVCF program and makes and oncoplot.
#'
#' @param path_to_parse A vector of paths to the excel files to be read.
#' @param tumor_only A logical value indicating if variant calling was performed in tumor only mode. Default is FALSE, indicating it was performed in paired mode.
#' @param cosmic A logical value indicating if the Cosmic columns should be used to filter (variants with both Cosmic columns empty are removed). Default is FALSE.
#' @param cosmicpref A string with the version prefix of the Cosmic columns in the input files. Default is "COSMIC99".
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
#' @param show_row_names Whether to show gene names. Default is TRUE.
#' @param show_pct Whether to show percentage of samples mutated per gene. Default is TRUE.
#' @param output Output file name (with png extension). Default is "oncoplot.png".
#' @param absent Vector of sample names (as they appear in the plot) that must have NO mutation in the genes shown. Used to compare a sample or a group of samples against the rest. Default is NULL.
#' @param topx Number of genes (rows) shown in the final oncoplot. Default is 50.
#' @param captureSize Size (in Mb) of the sequenced region used to compute the TMB. Default is 50 (exome); use the size of your panel for targeted panels.
#' @param draw_oncoplot Whether the intermediate maftools oncoplot (only needed to write the oncomatrix) is drawn on screen. FALSE saves time with many genes. Default is TRUE.
#' 
#' @return An oncoplot object which is also saved as a png file.
#'
#'
#' @examples
#' oncobuddy <-fromParse2Onco(path_to_parse = c("/path/to/file1","/path/to/file2"), VAF_tumor=0.05, minimalMutations=3, Missense_color="#FF5733", output ="my_oncoplot.png")
#' 
#' @export
fromParse2Onco <- function(path_to_parse, tumor_only = FALSE,
                           filter_column=c("PASS"), VAF_tumor=0, VAF_control=0, total_tumor_reads=0, alt_tumor_reads=0, cgi=FALSE, oncokb=FALSE, cosmic = FALSE, cosmicpref="COSMIC99", cgi_list=.default_cgi_list, oncokb_list=c("Likely Oncogenic", "Oncogenic"), annott=c("HIGH","MODERATE","MODIFIER"), tumor_samples_out = NULL, control_samples_out = NULL,
                           remove=TRUE, flags=FALSE, minimalMutations = 2, topgenes = 20000, nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),
                           Missense_color="#2a9134", Nonsense_color="#ffca3a", Nonstop_color="#000000", FrameDel_color="blue", FrameIns_color="purple", In_Frame_Ins_color="lightblue", In_Frame_Del_color="plum1", Translation_Start_Site_color="#ff0a54", Splice_site_color="darkorange", Multihit_color="#dab49d", show_row_names=TRUE, show_pct=TRUE, output="oncoplot.png", absent = NULL, topx = 50,
                           captureSize = 50, draw_oncoplot = TRUE) {

  # Read parse_variants output and convert to MAF
  variants_df <- fromParse2MAF(path_to_parse = path_to_parse, tumor_only = tumor_only,
                               oncokb = oncokb, cgi = cgi, cosmic = cosmic, cosmicpref = cosmicpref)

  # Filter MAF data frame
  filtered_df <- filterMAF(to_filter = variants_df, tumor_only = tumor_only, filter_column = filter_column,
                           VAF_tumor = VAF_tumor, VAF_control = VAF_control,
                           total_tumor_reads = total_tumor_reads, alt_tumor_reads = alt_tumor_reads,
                           cgi = cgi, oncokb = oncokb, cosmic = cosmic,
                           cgi_list = cgi_list, oncokb_list = oncokb_list, annott = annott,
                           tumor_samples_out = tumor_samples_out, control_samples_out = control_samples_out)

  # Create oncomatrix from MAF data frame
  prepareForOncoplot(maf_df = filtered_df, remove = remove, flags = flags,
                     minimalMutations = minimalMutations, topgenes = topgenes, nonSyn = nonSyn,
                     captureSize = captureSize, draw_oncoplot = draw_oncoplot)

  # Make oncoplot
  oncoplot <- makeOncoplot(Missense_color = Missense_color, Nonsense_color = Nonsense_color,
                           Nonstop_color = Nonstop_color, FrameDel_color = FrameDel_color,
                           FrameIns_color = FrameIns_color, In_Frame_Ins_color = In_Frame_Ins_color,
                           In_Frame_Del_color = In_Frame_Del_color,
                           Translation_Start_Site_color = Translation_Start_Site_color,
                           Splice_site_color = Splice_site_color, Multihit_color = Multihit_color,
                           show_row_names = show_row_names, show_pct = show_pct,
                           output = output, topx = topx, absent = absent)

  return(oncoplot)
}
