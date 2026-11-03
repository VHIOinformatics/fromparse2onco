#' Read parseVCF Output and Convert To MAF
#'
#' This function reads variant data from one or more Excel files produced by parseVCF program and processes it in order to match MAF format specifications.
#'
#' @param path_to_parse A vector of paths to the excel files to be read.
#' @param tumor_only A logical value indicating if variant calling was performed in tumor only mode. Default is FALSE, indicating it was performed in paired mode.
#' @param oncokb A logical value indicating whether OncoKB annotation was performed. Default is FALSE.
#' @param cgi A logical value indicating whether CGI annotation was performed. Default is FALSE.
#' @param cosmic A logical value indicating if the Cosmic columns should be used to filter. Default is FALSE
#' @param cosmicpref A string value indicating the version of Cosmic used to annotate the files. Default is COSMIC99
#'
#' @return A data frame containing the processed variant data in MAF format.
#' @import readxl
#' @import dplyr
#' @import stringr
#'
#' @examples
#' variants_df <- fromParse2MAF(c("/path/to/file1","/path/to/file2"), tumor_only=TRUE)
#'
#' @export
fromParse2MAF <- function(path_to_parse, tumor_only = FALSE, oncokb = FALSE, cgi = FALSE, cosmic = FALSE, cosmicpref = "COSMIC99") {

  # ---- 1. Check input files ----------------------------------------------
  missing_files <- path_to_parse[!file.exists(path_to_parse)]
  if (length(missing_files) > 0) {
    stop("These files do not exist: ", paste(missing_files, collapse = ", "), call. = FALSE)
  }

  # ---- 2. Read all Excel files --------------------------------------------
  # col_types = "text": every column is read as text, so all files get the same
  # types (readxl would otherwise guess the type from the first 1000 rows of each
  # file and the guess can differ between files). The numeric columns are
  # converted explicitly below.
  # bind_rows() (unlike rbind) matches columns by name and fills with NA the
  # columns that are missing in some file.
  variants_df <- dplyr::bind_rows(
    lapply(path_to_parse, readxl::read_excel, col_types = "text")
  )

  # ---- 3. Column maps ------------------------------------------------------
  # Named vectors where   names = NEW name (MAF)   and   values = OLD name (parseVCF)
  # Columns shared by both modes:
  common_map <- c(
    Chromosome             = "CHROM",
    Start_Position         = "POS",
    Reference_Allele       = "REF",
    Tumor_Seq_Allele2      = "ALT",
    Hugo_Symbol            = "Gene_Name",
    Variant_Classification = "Annotation",
    IMPACT                 = "Annotation_Impact",
    ExonicFunction_refGene = "Transcript_BioType",
    BNChange               = "HGVS.c",
    AAChange_refGene       = "HGVS.p"
  )

  if (tumor_only) {
    numeric_columns <- c("POS", "ADref", "ADalt", "DPtotal", "VAF")
    mode_map <- c(
      Tumor_Sample_Barcode         = "Sample",
      ReferenceReads_TumorSample   = "ADref",
      AlternativeReads_TumorSample = "ADalt",
      TotalReads                   = "DPtotal"
    )
    # Note: FILTER and VAF keep their name in this mode, no renaming is needed.
  } else {
    numeric_columns <- c(
      "POS", "Tumor_Reference_Reads", "Tumor_Alternative_Reads", "Tumor_Total_Reads", "Tumor_VAF_Alternative",
      "Control_Reference_Reads", "Control_Alternative_Reads", "Control_Total_Reads", "Control_VAF_Alternative"
    )
    mode_map <- c(
      Tumor_Sample_Barcode         = "Tumor",
      TotalReads_Normal            = "Control_Total_Reads",
      ReferenceReads_NormalSample  = "Control_Reference_Reads",
      AlternativeReads_NormalSample = "Control_Alternative_Reads",
      VAF_normal                   = "Control_VAF_Alternative",
      TotalReads_Tumor             = "Tumor_Total_Reads",
      ReferenceReads_TumorSample   = "Tumor_Reference_Reads",
      AlternativeReads_TumorSample = "Tumor_Alternative_Reads",
      VAF                          = "Tumor_VAF_Alternative"
    )
  }
  rename_map <- c(common_map, mode_map)

  # Optional annotation columns
  if (cgi) {
    rename_map <- c(rename_map,
                    "CGI-SUMMARY"    = "CGI-Oncogenicity summary",
                    "CGI-PREDICTION" = "CGI-Oncogenicity prediction")
  }
  if (oncokb) {
    rename_map <- c(rename_map, "OncoKB" = "OncoKB")
  }
  # Filter by Cosmic presence: set cosmic = TRUE and give the Cosmic version prefix
  # in cosmicpref (e.g. "COSMIC99"). Whatever the version, the columns are always
  # renamed to COSMIC_Coding / COSMIC_nonCoding so the filter does not depend on it.
  if (cosmic) {
    rename_map <- c(rename_map,
                    COSMIC_Coding    = paste0(cosmicpref, "_Coding"),
                    COSMIC_nonCoding = paste0(cosmicpref, "_nonCoding"))
  }

  # ---- 4. Check that every needed column exists (clear error message) -----
  needed  <- unique(c(unname(rename_map), numeric_columns))
  absent  <- setdiff(needed, names(variants_df))
  if (length(absent) > 0) {
    stop("These columns were not found in the input file(s): ",
         paste(absent, collapse = ", "),
         ". Check tumor_only / cgi / oncokb / cosmic / cosmicpref.", call. = FALSE)
  }

  # ---- 5. Numeric conversion and renaming ---------------------------------
  variants_df <- variants_df %>%
    mutate(across(all_of(numeric_columns), as.numeric)) %>%
    rename(all_of(rename_map))

  # ---- 6. End_Position and Variant_Type -----------------------------------
  variants_df <- variants_df %>%
    mutate(
      ref_len = str_length(Reference_Allele),
      alt_len = str_length(Tumor_Seq_Allele2),
      dif_len = alt_len - ref_len,
      # MAF specification (https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/):
      # End_Position is the highest genomic position of the variant. For a deletion
      # that is given by the REF length (before it used the ALT length, which is
      # wrong for deletions); for SNP/MNP/INS it is the same as before.
      # It is not used to draw the oncoplot, but read.maf requires it.
      End_Position = Start_Position + pmax(ref_len, alt_len) - 1,
      Variant_Type = case_when(
        dif_len == 0 & alt_len == 1 ~ "SNP",
        dif_len == 0 & alt_len == 2 ~ "DNP",
        dif_len == 0 & alt_len == 3 ~ "TNP",
        dif_len == 0 & alt_len > 3  ~ "ONP",
        dif_len > 0                 ~ "INS",
        dif_len < 0                 ~ "DEL"
      ),
      # Join the two UTR-truncation terms (they come joined by '&') before keeping
      # only the first consequence
      Variant_Classification = str_replace(Variant_Classification,
                                           "prime_UTR_truncation&exon_loss_variant",
                                           "prime_UTR_truncation+exon_loss_variant"),
      # First consequence only (SnpEff lists them from most to least severe)
      consequence = str_remove(Variant_Classification, "&.*")
    )

  # ---- 7. Variant_Classification (MAF vocabulary) --------------------------
  variants_df <- .classify_variants(variants_df)

  # Remove helper columns
  variants_df <- variants_df[, setdiff(names(variants_df), c("ref_len", "alt_len", "dif_len", "consequence"))]

  return(variants_df)
}


# Rules to translate SnpEff consequences into MAF Variant_Classification.
# ORDER MATTERS: the first rule that matches wins.
#   class   : MAF Variant_Classification assigned
#   pattern : regular expression searched in the (first) SnpEff consequence
#   vtype   : Variant_Type the variant must have (NA = any)
.classification_rules <- data.frame(
  class = c("Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
            "Nonstop_Mutation", "Translation_Start_Site", "In_Frame_Ins", "In_Frame_Del",
            "Missense_Mutation", "Intron", "Silent", "RNA", "5'UTR", "3'UTR", "IGR",
            "5'Flank", "3'Flank", "Gene_Fusion"),
  pattern = c(
    "splice_acceptor_variant|splice_donor_variant|transcript_ablation|exon_loss_variant|5_prime_UTR_truncation\\+exon_loss_variant|3_prime_UTR_truncation\\+exon_loss_variant",
    "stop_gained",
    "frameshift_variant",
    "frameshift_variant",
    "stop_lost",
    "initiator_codon_variant|start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant|coding_sequence_variant|rare_amino_acid_variant",
    "transcript_amplification|intron_variant|intragenic_variant",
    "synonymous_variant|stop_retained_variant|NMD_transcript_variant|start_retained",
    "mature_miRNA_variant|exon_variant|non_coding_exon_variant|non_coding_transcript|nc_transcript_variant",
    "5_prime_UTR_variant|5_prime_UTR_premature_start_codon_gain_variant",
    "3_prime_UTR_variant",
    "TF_binding_site_variant|regulatory_region|intergenic",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "gene_fusion"),
  vtype = c(NA, NA, "DEL", "INS", NA, NA, "INS", "DEL",
            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

# Applies .classification_rules. Needs the columns `consequence` and `Variant_Type`.
# Variants that match no rule keep their original Variant_Classification.
#' @noRd
.classify_variants <- function(df) {
  assigned <- rep(NA_character_, nrow(df))
  for (i in seq_len(nrow(.classification_rules))) {
    rule <- .classification_rules[i, ]
    hit  <- is.na(assigned) &
      grepl(rule$pattern, df$consequence) &
      (is.na(rule$vtype) | df$Variant_Type %in% rule$vtype)
    assigned[hit] <- rule$class
  }
  df$Variant_Classification <- dplyr::coalesce(assigned, df$Variant_Classification)
  df
}
