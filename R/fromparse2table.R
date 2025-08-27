#' Parse and Combine Variant Data
#'
#' This function reads and processes variant data from one or more Excel files, filtering and formatting the data as specified.
#'
#' @param path_to_parse A character string specifying the path to the directory or file to be parsed.
#' @param more_than_one A logical value indicating if there are multiple files to parse. Default is FALSE.
#' @param tumor_only A logical value indicating if only tumor data should be processed. Default is FALSE.
#' @param pattern_excel A character string specifying the pattern to match file names if there are multiple files. Default is FALSE.
#'
#' @return A data frame containing the parsed and processed variant data.
#' @import readxl
#' @import dplyr
#' @import stringr
#' @import purrr
#' @examples
#' # Example usage:
#' # variants_df <- fromparse2table("/path/to/files", more_than_one=TRUE, tumor_only=TRUE, pattern_excel="*.xlsx")
#'
#' @export
fromparse2table <- function(path_to_parse, more_than_one, tumor_only, pattern_excel=FALSE, oncokb=FALSE, cgi=FALSE) {
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
    numeric_columns <- c("POS", "ADref", "ADalt", "DPtotal", "VAF")

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
      "TotalReads" = "DPtotal",
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

  return(variants_df)
}

