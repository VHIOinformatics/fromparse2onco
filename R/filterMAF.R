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
#'
#' @examples
#' filtered_df <- filterMAF(variants_df, tumor_only=TRUE, VAF_tumor=0.05, alt_tumor_reads=10)
#'
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
  # Remove rows with NA in Variant_Classification
  filter_conditions <- filter_conditions %>%
    filter(!is.na(Variant_Classification))  # Remove rows with NA classifications
  # Return the filtered table
  return(filter_conditions)
}

