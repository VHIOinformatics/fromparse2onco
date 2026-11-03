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
#' @param topx Number of rows of the desired matrix. Default is 50. Strongly recommended.
#' @param absent This parameter filters genes that are mutation-free in the specified samples but may be mutated in others. Use the sample names as they appear in the plot (after removing repeated names).
#' @param matrix_file Path to the oncomatrix written by prepareForOncoplot. Default is "onco_matrix.txt" (in the working directory).
#' @param width,height,res Size (pixels) and resolution of the png file. Default is 2000 x 1200 at 150 dpi.
#'
#' @return An oncoplot object which is also saved as a png file.
#'
#' @import grid
#' @import ComplexHeatmap
#' @importFrom utils head read.table
#' @importFrom grDevices png
#'
#' @examples
#' oncoplot <- makeOncoplot(Missense_color="#FF5733", Nonsense_color="#33FF57", output = "my_oncoplot.png")
#'
#' @export

makeOncoplot <- function(Missense_color="#2a9134", Nonsense_color="#ffca3a", Nonstop_color="#000000", FrameDel_color="blue", FrameIns_color="purple", In_Frame_Ins_color="lightblue", In_Frame_Del_color="plum1", Translation_Start_Site_color="#ff0a54", Splice_site_color="darkorange", Multihit_color="#dab49d", show_row_names = TRUE, show_pct = TRUE, output="oncoplot.png", topx = 50, absent = NULL, matrix_file = "onco_matrix.txt", width = 2000, height = 1200, res = 150) {

  if (!file.exists(matrix_file)) {
    stop("Oncomatrix file not found: '", matrix_file, "'. Run prepareForOncoplot() first ",
         "(it writes onco_matrix.txt in the working directory).", call. = FALSE)
  }

  # Read the oncoplot matrix.
  # check.names = FALSE keeps the sample names exactly as they are in the file
  # (by default read.table changes '-' into '.' and adds an 'X' to names that start
  # with a number, which had to be undone afterwards and also changed the dots
  # that really belonged to the sample name).
  # colClasses = "character": an all-empty column would otherwise be read as logical NA.
  onco.matrix <- as.matrix(read.table(matrix_file, header = TRUE, sep = "\t", quote = "",
                                      check.names = FALSE, colClasses = "character"))

  # Replace specific strings in matrix
  onco.matrix <- gsub("Frame_Shift", "Frameshift", onco.matrix)

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

  # One grey background + one rectangle per mutation type, built from `col`
  # (before: 17 hand-written alter_graphic() calls).
  # The ORDER of alter_fun decides how the alterations are stacked in the bar
  # plots and in which order they appear in the legend, so it is kept explicit.
  alter_order <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Multi_Hit",
                   "Frameshift_Del", "Frameshift_Ins", "Translation_Start_Site", "Splice_Site",
                   "In_Frame_Ins", "In_Frame_Del", "Intron", "IGR", "3'Flank", "RNA",
                   "5'Flank", "3'UTR", "5'UTR")
  alter_fun <- c(
    list(background = alter_graphic("rect", fill = "#CCCCCC")),
    lapply(col[alter_order], function(color) alter_graphic("rect", width = 1, height = 1, fill = color))
  )

  # Some sample names come repeated (e.g. "A_B_A_B"): keep only one copy ("A_B")
  fix_repeat <- function(samplename) {
    parts <- strsplit(samplename, "_")[[1]]
    n <- length(parts)
    if (n == 0 || n %% 2 != 0) return(samplename)
    half <- n / 2
    if (identical(parts[1:half], parts[(half + 1):n])) {
      return(paste(parts[(half + 1):n], collapse = "_"))
    }
    samplename
  }
  colnames(onco.matrix) <- vapply(colnames(onco.matrix), fix_repeat, character(1), USE.NAMES = FALSE)

  # absent: keep only the genes with NO mutation in the samples given in `absent`
  # (to compare a sample or a group of samples against the rest).
  if (!is.null(absent) && length(absent) > 0) {
    not_found <- setdiff(absent, colnames(onco.matrix))
    if (length(not_found) > 0) {
      stop("These samples in 'absent' are not in the oncomatrix: ",
           paste(not_found, collapse = ", "), call. = FALSE)
    }
    rows_to_keep <- rowSums(onco.matrix[, absent, drop = FALSE] == "") == length(absent)
    onco.matrix <- onco.matrix[rows_to_keep, , drop = FALSE]
    if (nrow(onco.matrix) == 0) {
      stop("No genes left after applying 'absent'.", call. = FALSE)
    }
  }

  # topx: keep only the topx most frequently mutated genes
  freq_onco_matrix <- rowSums(onco.matrix != "")
  topx_mat <- head(onco.matrix[order(freq_onco_matrix, decreasing = TRUE), , drop = FALSE], topx)

  # Execute oncoPrint
  p <- ComplexHeatmap::oncoPrint(mat = topx_mat, col = col,
                                 alter_fun = alter_fun,
                                 show_row_names = show_row_names,
                                 pct_side = "right",
                                 row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                                 show_pct = show_pct,
                                 column_names_side = c("bottom"),
                                 show_column_names = TRUE,
                                 show_heatmap_legend = TRUE,
                                 column_order = order(colSums(topx_mat != ""), decreasing = TRUE),
                                 row_names_side = "left")

  png(output, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)   # closes the file even if draw() fails
  draw(p)
  return(p)
}
