#' Make Oncoplot
#'
#' This function imports an onco_matrix.txt file in the working directory and makes an oncoplot using the oncoPrint function in ComplexHeatmap package.
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
#'
#' @return An oncoplot object which is also saved as a png file.
#'
#' @import readxl
#' @import dplyr
#' @import stringr
#' @import purrr
#' @import grid
#' @import ComplexHeatmap
#'
#' @examples
#' oncoplot <- makeOncoplot(Missense_color="#FF5733", Nonsense_color="#33FF57", output = "my_oncoplot.png")
#'
#' @export

makeOncoplot <- function(Missense_color="#2a9134", Nonsense_color="#ffca3a", Nonstop_color="#000000", FrameDel_color="blue", FrameIns_color="purple", In_Frame_Ins_color="lightblue", In_Frame_Del_color="plum1", Translation_Start_Site_color="#ff0a54", Splice_site_color="darkorange", Multihit_color="#dab49d", show_row_names = TRUE, show_pct = TRUE, output="oncoplot.png") {

  # Read the oncoplot matrix
  onco.matrix <- as.matrix(read.table("onco_matrix.txt", header = TRUE, sep = '\t', quote = ""))
  
  # Replace specific strings in matrix
  onco.matrix <- gsub('Frame_Shift', 'Frameshift', onco.matrix)
  colnames(onco.matrix) <- gsub("\\.","-", colnames(onco.matrix))
  
  # Define colors for each type of mutation
  col <- c(Missense_Mutation = Missense_color,
           Nonsense_Mutation = Nonsense_color,
           Nonstop_Mutation = Nonstop_color,
           Frameshift_Del = FrameDel_color,
           Frameshift_Ins = FrameIns_color,
           In_Frame_Ins = In_Frame_Ins_color,
           In_Frame_Del = In_Frame_Del_color,
           Translation_Start_Site = Translation_Start_Site_color,
           Splice_Site = Splice_site_color,
           Multi_Hit = Multihit_color,
	   Intron = "maroon",
           IGR = "aquamarine",
           `3'Flank` = "lightsalmon3",
           `5'Flank` = "yellow2",
           RNA = "yellowgreen",
           `3'UTR` = "firebrick2",
           `5'UTR` = "mediumpurple1")

  # Assign shapes for each type of mutation
  alter_fun <- list(
    background = alter_graphic("rect", fill = "#CCCCCC"),
    Missense_Mutation = alter_graphic("rect",
                                      width = 1,
                                      height = 1,
                                      fill = col["Missense_Mutation"]),
    Nonsense_Mutation = alter_graphic("rect",
                                      width = 1,
                                      height = 1,
                                      fill = col["Nonsense_Mutation"]),
    Nonstop_Mutation = alter_graphic("rect",
                                     width = 1,
                                     height = 1,
                                     fill = col["Nonstop_Mutation"]),
    Multi_Hit = alter_graphic("rect",
                              width = 1,
                              height = 1,
                              fill = col["Multi_Hit"]),
    Frameshift_Del = alter_graphic("rect",
                                   width = 1,
                                   height = 1,
                                   fill = col["Frameshift_Del"]),
    Frameshift_Ins = alter_graphic("rect",
                                   width = 1,
                                   height = 1,
                                   fill = col["Frameshift_Ins"]),
    Translation_Start_Site = alter_graphic("rect",
                                           width = 1,
                                           height = 1,
                                           fill = col["Translation_Start_Site"]),
    Splice_Site = alter_graphic("rect",
                                width = 1,
                                height = 1,
                                fill = col["Splice_Site"]),
    In_Frame_Ins = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["In_Frame_Ins"]),
    In_Frame_Del = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["In_Frame_Del"]),
    Intron = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["Intron"]),
    IGR = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["IGR"]),
    `3'Flank` = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["3'Flank"]),
    RNA = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["RNA"]),
    `5'Flank` = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["5'Flank"]),
    `3'UTR` = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["3'UTR"]),
    `5'UTR` = alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["5'UTR"]))
  
  #Execute oncoPrint
  p <- ComplexHeatmap::oncoPrint(mat = onco.matrix, col = col, 
                                 alter_fun = alter_fun, alter_fun_is_vectorized = FALSE, 
                                 show_row_names = show_row_names,
                                 pct_side = "right",
                                 row_names_gp = gpar(fontsize = 10,fontface = "italic"),
                                 show_pct = show_pct,
                                 column_names_side = c("bottom"), 
                                 show_column_names = TRUE,
                                 show_heatmap_legend = TRUE,
                                 column_order = order(apply(onco.matrix,2,function(x){length(which(x != ""))} ), decreasing = TRUE),
                                 row_names_side = "left")
  png(output, width = 2000, height = 1200, res = 150)
  draw(p)
  dev.off()
  
  return(p)
}

