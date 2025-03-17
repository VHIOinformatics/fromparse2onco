#' Convert Parsed Data to Oncoplot
#'
#' @param path_to_parse Path to the parsed data file(s).
#' @param more_than_one Logical indicating whether there are multiple files to process.
#' @param tumor_only Logical indicating whether only tumor data should be considered.
#' @param pattern_excel Pattern for matching Excel files.
#' @param filter_column Filter column for variant filtering.
#' @param VAF_tumor Minimum VAF for tumor samples.
#' @param VAF_control Minimum VAF for control samples.
#' @param total_reads Minimum total reads.
#' @param tumor_reads Minimum tumor reads.
#' @param remove Logical indicating whether duplicated variants should be removed.
#' @param minimum_samples_mutated Minimum number of samples mutated for inclusion in the oncoplot.
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
#' @param minimalMutations Minimum number of mutations per gene to not be discarded. Default is 10.
#' @param topgenes Number of genes to show in the oncoplot
#'
#' @return An object to represent an oncoplot or an oncoplot
#' @export
#'
#' @import readxl
#' @import dplyr
#' @import maftools
#' @import ComplexHeatmap
#' @import stringr
#' @import ggplot2
#' @import tidyverse
#'
#'
#' @examples
#' oncobuddy <-fromparse2onco("/path/to/parse/",TRUE,TRUE,pattern_excel="pattern")
fromparse2onco <- function(path_to_parse, more_than_one, tumor_only, pattern_excel=FALSE, filter_column=c("PASS"), VAF_tumor=0, VAF_control=0,
                           total_reads=0,tumor_reads=0,remove=FALSE,minimum_samples_mutated=1,Missense_color="#2a9134",Nonsense_color="#ffca3a",
                           Nonstop_color="#000000",FrameDel_color="blue",FrameIns_color="purple",In_Frame_Ins_color="lightblue",In_Frame_Del_color="plum1",
                           Tranlsation_Start_Site_color="#ff0a54",Splice_site_color="darkorange",Multihit_color="#dab49d", flags=FALSE,
                           nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),
                           cgi_path=FALSE, oncokb_path=FALSE, cgi_list=c("oncogenic (predicted)", "oncogenic (predicted and annotated)"), oncokb_list=c("Likely Oncogenic", "Oncogenic"),
                           annott=c("HIGH","MODERATE","MODIFIER"),oncokb=FALSE, cgi=FALSE, minimalMutations = 10, topgenes = 50) {
  # Read the variants data
  variants_df <- NULL
  if (more_than_one) {
    print("Processing multiple files")
    results_parse_more_than_one <- list.files(path = path_to_parse, pattern = pattern_excel)

    # Use lapply to read and bind all files into one dataframe
    variants_df <- do.call(bind_rows, lapply(results_parse_more_than_one, function(i) {
      read_excel(file.path(path_to_parse, i))
    }))
  } else {
    print("Processing a single file")
    variants_df <- read_excel(path_to_parse)
  }

  # Process data based on tumor_only flag
  if (tumor_only) {
    # Columns to convert to numeric for tumor-only mode
    numeric_columns <- c("POS", "ADref", "ADalt", "DPtootal", "VAF")

    # Mapping of column names to new names for tumor-only mode
    rename_map <- list(
      "Tumor_Sample_Barcode" = "Sample",
      "Chromosome" = "CHROM",
      "Start_Position" = "POS",
      "Reference_Allele" = "REF",
      "Tumor_Seq_Allele2" = "ALT",
      "FILTER" = "FILTER",
      "Variant_Type" = "Feature_Type",
      "Hugo_Symbol" = "Gene_Name",
      "Variant_Classification" = "Annotation",
      "IMPACT" = "Annotation_Impact",
      "ExonicFunction_refGene" = "Transcript_BioType",
      "BNChange" = "HGVS.c",
      "AAChange_refGene" = "HGVS.p",
      "ReferenceReads_TumorSample" = "ADref",
      "AlternativeReads_TumorSample" = "ADalt",
      "TotalReads" = "DPtootal",
      "VAF" = "VAF"
    )
  } else {
    # Columns to convert to numeric for paired tumor-normal mode
    numeric_columns <- c(
      "POS", "Tumor_Reference_Reads", "Tumor_Alternative_Reads", "Tumor_Total_Reads", "Tumor_VAF_Alternative",
      "Control_Reference_Reads", "Control_Alternative_Reads", "Control_Total_Reads", "Control_VAF_Alternative"
    )

    # Mapping of column names to new names for paired tumor-normal mode
    rename_map <- list(
      "Tumor_Sample_Barcode" = "Tumor",
      "Chromosome" = "CHROM",
      "Start_Position" = "POS",
      "Reference_Allele" = "REF",
      "Tumor_Seq_Allele2" = "ALT",
      "FILTER" = "FILTER",
      "Variant_Type" = "Feature_Type",
      "Hugo_Symbol" = "Gene_Name",
      "Variant_Classification" = "Annotation",
      "IMPACT" = "Annotation_Impact",
      "ExonicFunction_refGene" = "Transcript_BioType",
      "BNChange" = "HGVS.c",
      "AAChange_refGene" = "HGVS.p",
      "TotalReads_Normal" = "Control_Total_Reads",
      "ReferenceReads_NormalSample" = "Control_Reference_Reads",
      "AlternativeReads_NormalSample" = "Control_Alternative_Reads",
      "VAF_normal" = "Control_VAF_Alternative",
      "TotalReads_Tumor" = "Tumor_Total_Reads",
      "ReferenceReads_TumorSample" = "Tumor_Reference_Reads",
      "AlternativeReads_TumorSample" = "Tumor_Alternative_Reads",
      "VAF" = "Tumor_VAF_Alternative"
    )
  }

  # Add optional mappings for CGI and OncoKB data if enabled
  if (cgi) {
    rename_map <- c(rename_map,
                    "CGI-SUMMARY" = "CGI-Oncogenic Summary",
                    "CGI-PREDICTION" = "CGI-Oncogenic Prediction")
  }
  if (oncokb) {
    rename_map <- c(rename_map, "OncoKB" = "OncoKB")
  }
  # Apply transformations:
  # 1. Convert specified columns to numeric (if they exist).
  # 2. Rename columns based on the rename_map.
  variants_df <- variants_df %>%
    mutate(across(all_of(numeric_columns), as.numeric)) %>%  # Convert numeric columns
    rename(!!!rename_map)  # Rename columns using rename_map


  # Calculate End_Position and modify Variant_Type column
  variants_df <- variants_df %>%
    mutate(dif_len = str_length(Tumor_Seq_Allele2) - str_length(Reference_Allele),
           End_Position = Start_Position + abs(dif_len) - 1,
           Variant_Type = case_when(
             dif_len == 0 & str_length(Tumor_Seq_Allele2) == 1 ~ 'SNP',
             dif_len > 0 ~ 'INS',
             dif_len < 0 ~ 'DEL',
             TRUE ~ Variant_Type
           ))
  # Initialize the filtered table with base conditions
  filter_conditions <- variants_df %>%
    filter(FILTER %in% filter_column,
           VAF >= VAF_tumor,
           AlternativeReads_TumorSample >= tumor_reads,
           IMPACT %in% annott)

  # Add additional conditions if considering only tumor
  # Define a function to apply additional filters
  apply_additional_filters <- function(df, tumor_only, cgi_path, oncokb_path) {
    if (tumor_only) {
      df <- df %>%
        filter(TotalReads >= total_reads)
    } else {
      df <- df %>%
        filter(
          VAF_normal >= VAF_control,
          TotalReads_Tumor >= total_reads
        )
    }
    if (cgi_path) {
      df <- df %>%
        filter(`CGI-SUMMARY` %in% cgi_list)
    }
    if (oncokb_path) {
      df <- df %>%
        filter(OncoKB %in% oncokb_list)
    }
    return(df)
  }

  # Replace specific strings in the Annotation column
  filter_conditions <- filter_conditions %>%
    mutate(
      Variant_Classification = str_replace(Variant_Classification,
                                           "prime_UTR_truncation&exon_loss_variant",
                                           "prime_UTR_truncation+exon_loss_variant"),
      # Create a new column 'disgreggation' by removing '&' and subsequent characters
      disgreggation = str_remove(Variant_Classification, "&.*")
    )

  # Define a named list of patterns and their corresponding classifications
  classification_patterns <- list(
    "Splice_Site" = "splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant|5_prime_UTR_truncation\\+exon_loss_variant|3_prime_UTR_truncation\\+exon_loss_variant",
    "Nonsense_Mutation" = "stop_gained",
    "Frame_Shift_Del" = "frameshift_variant",
    "Frame_Shift_Ins" = "frameshift_variant",
    "Nonstop_Mutation" = "stop_lost",
    "Translation_Start_Site" = "initiator_codon_variant|start_lost",
    "In_Frame_Ins" = "inframe_insertion",
    "In_Frame_Del" = "inframe_deletion",
    "Missense_Mutation" = "missense_variant|coding_sequence_variant|rare_amino_acid_variant",
    "Intron" = "transcript_amplification|intron_variant|intragenic_variant",
    "Silent" = "synonymous_variant|stop_retained_variant|NMD_transcript_variant|start_retained",
    "RNA" = "mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript|nc_transcript_variant",
    "5'UTR" = "5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant",
    "3'UTR" = "3_prime_UTR_variant",
    "IGR" = "TF_binding_site_variant|regulatory_region|intergenic",
    "5'Flank" = "upstream_gene_variant",
    "3'Flank" = "downstream_gene_variant",
    "Gene_Fusion" = "gene_fusion"
  )

  # Apply patterns to classify Variant_Classification
  filter_conditions <- filter_conditions %>%
    mutate(
      Variant_Classification = case_when(
        grepl(classification_patterns[["Splice_Site"]], disgreggation) ~ "Splice_Site",
        grepl(classification_patterns[["Nonsense_Mutation"]], disgreggation) ~ "Nonsense_Mutation",
        grepl(classification_patterns[["Frame_Shift_Del"]], disgreggation) & Variant_Type == "DEL" ~ "Frame_Shift_Del",
        grepl(classification_patterns[["Frame_Shift_Ins"]], disgreggation) & Variant_Type == "INS" ~ "Frame_Shift_Ins",
        grepl(classification_patterns[["Nonstop_Mutation"]], disgreggation) ~ "Nonstop_Mutation",
        grepl(classification_patterns[["Translation_Start_Site"]], disgreggation) ~ "Translation_Start_Site",
        grepl(classification_patterns[["In_Frame_Ins"]], disgreggation) & Variant_Type == "INS" ~ "In_Frame_Ins",
        grepl(classification_patterns[["In_Frame_Del"]], disgreggation) & Variant_Type == "DEL" ~ "In_Frame_Del",
        grepl(classification_patterns[["Missense_Mutation"]], disgreggation) ~ "Missense_Mutation",
        grepl(classification_patterns[["Intron"]], disgreggation) ~ "Intron",
        grepl(classification_patterns[["Silent"]], disgreggation) ~ "Silent",
        grepl(classification_patterns[["RNA"]], disgreggation) ~ "RNA",
        grepl(classification_patterns[["5'UTR"]], disgreggation) ~ "5'UTR",
        grepl(classification_patterns[["3'UTR"]], disgreggation) ~ "3'UTR",
        grepl(classification_patterns[["IGR"]], disgreggation) ~ "IGR",
        grepl(classification_patterns[["5'Flank"]], disgreggation) ~ "5'Flank",
        grepl(classification_patterns[["3'Flank"]], disgreggation) ~ "3'Flank",
        grepl(classification_patterns[["Gene_Fusion"]], disgreggation) ~ "Gene_Fusion",
        TRUE ~ Variant_Classification
      )
    ) %>%
    filter(!is.na(Variant_Classification))  # Remove rows with NA classifications


  # Read and summarize MAF file
  maf_object <- maftools::read.maf(maf = filter_conditions, removeDuplicatedVariants = remove, rmFlags=flags,vc_nonSyn = nonSyn)
  maftools::plotmafSummary(maf = maf_object,
                           addStat = 'median',
                           titvRaw = FALSE,
                           showBarcodes = TRUE)

  # Print the number of unique genes
  print(length(unique(filter_conditions$Hugo_Symbol)))

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
  # Read the oncoplot matrix
  maf.matrix <- as.matrix(read.table("onco_matrix.txt", header = TRUE, sep = '\t', quote = ""))

  # Calculate mutation ratio
  mut_ratio <- apply(maf.matrix, 1, function(x) {y = (x != ''); sum(y)})

  # Filter the matrix by minimum samples mutated
  maf.matrix.filt <- maf.matrix[names(mut_ratio)[mut_ratio >= minimum_samples_mutated], ]

  # Replace specific strings in the filtered matrix
  maf.matrix.filt <- gsub('Frame_Shift', 'Frameshift', maf.matrix.filt)
  colnames(maf.matrix.filt) <- gsub("\\.","-", colnames(maf.matrix.filt))


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

  # Define the mutation types
  mutation_types <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                      "Frameshift_Del", "Frameshift_Ins", "In_Frame_Ins",
                      "In_Frame_Del", "Translation_Start_Site", "Splice_Site",
                      "Multi_Hit")

  # Create a list of mutation shapes dynamically
  alter_fun <- list(
    background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
    Missense_Mutation = ComplexHeatmap::alter_graphic("rect",
                                      width = 1,
                                      height = 1,
                                      fill = col["Missense_Mutation"]),
    Nonsense_Mutation = ComplexHeatmap::alter_graphic("rect",
                                      width = 1,
                                      height = 1,
                                      fill = col["Nonsense_Mutation"]),
    Nonstop_Mutation = ComplexHeatmap::alter_graphic("rect",
                                     width = 1,
                                     height = 1,
                                     fill = col["Nonstop_Mutation"]),
    Multi_Hit = ComplexHeatmap::alter_graphic("rect",
                              width = 1,
                              height = 1,
                              fill = col["Multi_Hit"]),
    Frameshift_Del = ComplexHeatmap::alter_graphic("rect",
                                   width = 1,
                                   height = 1,
                                   fill = col["Frameshift_Del"]),
    Frameshift_Ins = ComplexHeatmap::alter_graphic("rect",
                                   width = 1,
                                   height = 1,
                                   fill = col["Frameshift_Ins"]),
    Translation_Start_Site = ComplexHeatmap::alter_graphic("rect",
                                           width = 1,
                                           height = 1,
                                           fill = col["Translation_Start_Site"]),
    Splice_Site = ComplexHeatmap::alter_graphic("rect",
                                width = 1,
                                height = 1,
                                fill = col["Splice_Site"]),
    In_Frame_Ins = ComplexHeatmap::alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["In_Frame_Ins"]),
    In_Frame_Del = ComplexHeatmap::alter_graphic("rect",
                                 width = 1,
                                 height = 1,
                                 fill = col["In_Frame_Del"]))

  # Convert the list into a named list, keeping the mutation types as names
  #names(alter_fun)[-1] <- mutation_types

  p <- ComplexHeatmap::oncoPrint(mat = maf.matrix.filt,
                            col = col,
                            alter_fun = alter_fun,
                            alter_fun_is_vectorized = FALSE,
                            show_row_names = TRUE,
                            pct_side = TRUE,
                            row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                            show_pct = TRUE,
                            column_names_side = c("bottom"),
                            column_order = newOrder,
                            show_column_names = TRUE,
                            show_heatmap_legend = TRUE)
  png("oncoplot.png", width = 2000, height = 1200, res = 150)
  draw(p)
  dev.off()
  return(p)
}
