#' Filter MAF Data Frame
#'
#' This function filters a data frame of genetic variants based on specified criteria (VAF, reads & flags). It also allows the removal of tumor or control samples.
#'
#' @param to_filter A data frame containing the variants to be filtered in MAF format.
#' @param tumor_only A logical value indicating if the variant calling was performed in tumor-only mode. Default is FALSE, indicating it was performed in paired mode.
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
#' @return A filtered data frame containing variants that meet the specified criteria.
#' @examples
#' # Example usage:
#' # filtered_df <- filterMAF(variants_df, tumor_only=TRUE, VAF_tumor=0.05, alt_tumor_reads=10)
#' @export
filterMAF <- function(to_filter, tumor_only=FALSE, filter_column=c("PASS"), VAF_tumor=0, VAF_control=0, total_tumor_reads=0, alt_tumor_reads=0, cgi=FALSE, oncokb=FALSE, cgi_list=c("oncogenic (predicted)", "oncogenic (predicted and annotated)", "oncogenic (annotated)"), oncokb_list=c("Likely Oncogenic", "Oncogenic"), annott=c("HIGH","MODERATE","MODIFIER"), tumor_samples_out = NULL, control_samples_out = NULL) {
  # Initialize the filtered table with base conditions
  filter_conditions <- to_filter %>%
    # Filter by the specified columns and minimum conditions for VAF and tumor alternative reads
    filter(FILTER %in% filter_column,
           VAF >= VAF_tumor,
           AlternativeReads_TumorSample >= alt_tumor_reads,
           IMPACT %in% annott)

  # Add additional conditions if considering only tumor
  if (tumor_only) {
    filter_conditions <- filter_conditions %>%
      # Filter by the total number of reads
      filter(TotalReads >= total_tumor_reads)

    # If cgi is true, filter by the CGI list
    if (cgi) {
      filter_conditions <- filter_conditions %>%
        filter(`CGI-SUMMARY` %in% cgi_list)
    }

    # If oncokb is true, filter by the OncoKB list
    if (oncokb) {
      filter_conditions <- filter_conditions %>%
        filter(OncoKB %in% oncokb_list)
    }

    # If tumor_samples_out is not NULL, filter them out
    if (length(tumor_samples_out) > 0) {
      filter_conditions <- filter_conditions %>%
        filter(! Tumor_Sample_Barcode %in% tumor_samples_out)
    }



  } else {
    # If not only tumor, add conditions for normal reads and control VAF
    filter_conditions <- filter_conditions %>%
      filter(VAF_normal <= VAF_control,
             TotalReads_Tumor >= total_tumor_reads)

    # If cgi is true, filter by the CGI list
    if (cgi) {
      filter_conditions <- filter_conditions %>%
        filter(`CGI-SUMMARY` %in% cgi_list)
    }

    # If oncokb is true, filter by the OncoKB list
    if (oncokb) {
      filter_conditions <- filter_conditions %>%
        filter(OncoKB %in% oncokb_list)
    }

    # If tumor_samples_out is not NULL, filter them out
    if (length(tumor_samples_out) > 0) {
      filter_conditions <- filter_conditions %>%
        filter(! Tumor_Sample_Barcode %in% tumor_samples_out)
    }

    # If control_samples_out is not NULL, filter them out
    if (length(control_samples_out) > 0) {
      filter_conditions <- filter_conditions %>%
        filter(! Control %in% control_samples_out)
    }



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
  # Return the filtered table
  return(filter_conditions)
}

