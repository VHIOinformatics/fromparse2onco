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
#' oncoplot <- makeOncoplot(Missense_color="#FF5733", Nonsense_color="#33FF57")
#'
#' @export

makeOncoplot <- function(Missense_color="#2a9134", Nonsense_color="#ffca3a", Nonstop_color="#000000", FrameDel_color="blue", FrameIns_color="purple", In_Frame_Ins_color="lightblue", In_Frame_Del_color="plum1", Translation_Start_Site_color="#ff0a54", Splice_site_color="darkorange", Multihit_color="#dab49d") {

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
           Multi_Hit = Multihit_color)

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
                                 fill = col["In_Frame_Del"]))
  
  #Execute oncoPrint
  p <- ComplexHeatmap::oncoPrint(mat = onco.matrix, col = col, 
                                 alter_fun = alter_fun, alter_fun_is_vectorized = FALSE, 
                                 show_row_names = TRUE,
                                 pct_side = "right",
                                 row_names_gp = gpar(fontsize = 10,fontface = "italic"),
                                 show_pct = TRUE,
                                 column_names_side = c("bottom"), 
                                 show_column_names = TRUE,
                                 show_heatmap_legend = TRUE,
                                 column_order = order(apply(onco.matrix,2,function(x){length(which(x != ""))} ), decreasing = TRUE),
                                 row_names_side = "left")
  png("oncoplot.png", width = 2000, height = 1200, res = 150)
  draw(p)
  dev.off()
  
  return(p)
}

