# Default CGI labels kept by the filter. Shared with fromParse2Onco so both
# functions always use the same list (before, each one had a different one).
.default_cgi_list <- c("Oncogenic (predicted by BoostDM)",
                       "Oncogenic (predicted by RulesDM)",
                       "Oncogenic (annotated and predicted by BoostDM)",
                       "Potentially Oncogenic (predicted by RulesDM)")

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
#' @param cosmic A logical value indicating if the Cosmic columns should be used to filter. Default is FALSE
#'
#' @return A filtered data frame containing variants that meet the specified criteria.
#'
#' @examples
#' filtered_df <- filterMAF(variants_df, tumor_only=TRUE, VAF_tumor=0.05, alt_tumor_reads=10)
#'
#' @import dplyr
#'
#' @export
filterMAF <- function(to_filter, tumor_only=FALSE, filter_column=c("PASS"), VAF_tumor=0, VAF_control=0, total_tumor_reads=0, alt_tumor_reads=0, cgi=FALSE, oncokb=FALSE, cosmic=FALSE, cgi_list=.default_cgi_list, oncokb_list=c("Likely Oncogenic", "Oncogenic"), annott=c("HIGH","MODERATE","MODIFIER"), tumor_samples_out = NULL, control_samples_out = NULL) {

  # filter() silently drops the rows where the condition is NA. Tell the user
  # how many variants are going to disappear for that reason.
  numeric_cols <- if (tumor_only) c("VAF", "AlternativeReads_TumorSample", "TotalReads")
                  else c("VAF", "AlternativeReads_TumorSample", "TotalReads_Tumor", "VAF_normal")
  for (col in numeric_cols) {
    n_na <- sum(is.na(to_filter[[col]]))
    if (n_na > 0) message(n_na, " variants have NA in '", col, "' and will be removed by the filters.")
  }

  # Conditions common to both modes
  df <- filter(to_filter,
               FILTER %in% filter_column,
               VAF >= VAF_tumor,
               AlternativeReads_TumorSample >= alt_tumor_reads,
               IMPACT %in% annott)

  # Conditions that depend on the mode
  if (tumor_only) {
    df <- filter(df, TotalReads >= total_tumor_reads)
    if (length(control_samples_out) > 0) {
      warning("control_samples_out is ignored in tumor_only mode.", call. = FALSE)
    }
  } else {
    df <- filter(df,
                 VAF_normal <= VAF_control,
                 TotalReads_Tumor >= total_tumor_reads)
    if (length(control_samples_out) > 0) {
      df <- filter(df, !Control %in% control_samples_out)
    }
  }

  # Optional filters (identical in both modes, before they were written twice)
  if (cgi) {
    df <- filter(df, `CGI-SUMMARY` %in% cgi_list)
  }
  if (oncokb) {
    df <- filter(df, OncoKB %in% oncokb_list)
  }
  if (cosmic) {
    # fromParse2MAF always renames the Cosmic columns to COSMIC_Coding / COSMIC_nonCoding
    message("Applying Cosmic filter: removing rows with empty Cosmic columns")
    df <- filter(df,
                 (!is.na(COSMIC_Coding)    & COSMIC_Coding    != "") |
                 (!is.na(COSMIC_nonCoding) & COSMIC_nonCoding != ""))
  }
  if (length(tumor_samples_out) > 0) {
    df <- filter(df, !Tumor_Sample_Barcode %in% tumor_samples_out)
  }

  # Remove rows with NA in Variant_Classification
  df <- filter(df, !is.na(Variant_Classification))

  return(df)
}
