#' Configure Oncoplot Colors and Shapes
#'
#' This function configures the colors and shapes used in oncoplots for various types of mutations.
#'
#' @param Missense_color Color for missense variants. Default is "#2a9134".
#' @param Nonsense_color Color for nonsense variants. Default is "#ffca3a".
#' @param Nonstop_color Color for nonstop variants. Default is "#000000".
#' @param FrameDel_color Color for frame shift deletions. Default is "blue".
#' @param FrameIns_color Color for frame shift insertions. Default is "purple".
#' @param In_Frame_Ins_color Color for in-frame insertions. Default is "lightblue".
#' @param In_Frame_Del_color Color for in-frame deletions. Default is "plum1".
#' @param Tranlsation_Start_Site_color Color for translation start sites. Default is "#ff0a54".
#' @param Splice_site_color Color for splicing sites. Default is "darkorange".
#' @param Multihit_color Color for multi-hit genes. Default is "#dab49d".
#'
#' @return A list containing the color configuration and alteration function for the oncoplot.
#'
#' @import grid
#' @import ComplexHeatmap
#'
#' @examples
#' #color_config <- confoplot(Missense_color="#FF5733", Nonsense_color="#33FF57")
#'
#' @export
confoplot <- function(Missense_color="#2a9134", Nonsense_color="#ffca3a", Nonstop_color="#000000", FrameDel_color="blue", FrameIns_color="purple", In_Frame_Ins_color="lightblue", In_Frame_Del_color="plum1", Tranlsation_Start_Site_color="#ff0a54", Splice_site_color="darkorange", Multihit_color="#dab49d") {

  # Define colors for each type of mutation
  col <- c(Missense_Mutation = Missense_color,
           Nonsense_Mutation = Nonsense_color,
           Nonstop_Mutation = Nonstop_color,
           Frameshift_Del = FrameDel_color,
           Frameshift_Ins = FrameIns_color,
           In_Frame_Ins = In_Frame_Ins_color,
           In_Frame_Del = In_Frame_Del_color,
           Translation_Start_Site = Tranlsation_Start_Site_color,
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

  return(list(col=col, alter_fun=alter_fun))
}

