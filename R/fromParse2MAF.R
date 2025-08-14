#' Read parseVCF Output and Convert To MAF
#'
#' This function reads variant data from one or more Excel files produced by parseVCF program and processes it in order to match MAF format specifications.
#'
#' @param path_to_parse A vector of paths to the excel files to be read.
#' @param tumor_only A logical value indicating if variant calling was performed in tumor only mode. Default is FALSE, indicating it was performed in paired mode.
#' @param oncokb A logical value indicating whether OncoKB annotation was performed. Default is FALSE.
#' @param cgi A logical value indicating whether CGI annotation was performed. Default is FALSE.
#'
#' @return A data frame containing the processed variant data in MAF format.
#' @import readxl
#' @import dplyr
#' @import stringr
#' @import purrr
#'
#' @examples
#' variants_df <- fromParse2MAF(c("/path/to/file1","/path/to/file2"), tumor_only=TRUE)
#'
#' @export
fromParse2MAF <- function(path_to_parse,tumor_only = FALSE, oncokb=FALSE, cgi=FALSE) {
  # Read the variants data from one or multiple excel files
  variants_df <- do.call(rbind, lapply(path_to_parse, 
                                             function(i) {
                                               read_excel(file.path(i))
                                             }))


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
           #According to MAF format specifications (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) the end position is the highest numeric genomic position of the reported variant
           #The End_Position field is not actually used to make the oncoplot but it is required by the read.maf function when validating the required fields from MAF
           End_Position = Start_Position + str_length(Tumor_Seq_Allele2) - 1,
           Variant_Type = case_when(
             dif_len == 0 & str_length(Tumor_Seq_Allele2) == 1 ~ 'SNP',
             dif_len == 0 & str_length(Tumor_Seq_Allele2) == 2 ~ 'DNP',
             dif_len == 0 & str_length(Tumor_Seq_Allele2) == 3 ~ 'TNP',
             dif_len == 0 & str_length(Tumor_Seq_Allele2) > 3 ~ 'ONP',
             dif_len > 0 ~ 'INS',
             dif_len < 0 ~ 'DEL'
           ))
  
  # Replace specific strings in the Annotation (Variant_Classification) column
  variants_df <- variants_df %>%
    mutate(
      Variant_Classification = str_replace(Variant_Classification,
                                           "prime_UTR_truncation&exon_loss_variant",
                                           "prime_UTR_truncation+exon_loss_variant"),
      # Create a new column 'disgreggation' by removing '&' and subsequent characters
      disgreggation = str_remove(Variant_Classification, "&.*")
    )

  # Define a named list of patterns and their corresponding classifications
  classification_patterns <- list(
    "Splice_Site" = "splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant|5_prime_UTR_truncation\\+exon
_loss_variant|3_prime_UTR_truncation\\+exon_loss_variant",
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
  variants_df <- variants_df %>%
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
    ) 

  return(variants_df)
}

