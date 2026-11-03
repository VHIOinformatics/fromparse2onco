# Builds a small parseVCF-like Excel file for the tests
make_parse_xlsx <- function(path, tumor_only = TRUE) {
  df <- data.frame(
    CHROM = c("chr1", "chr1", "chrX", "chr2", "chr3"),
    POS   = c(100, 200, 300, 400, 500),
    REF   = c("A", "ACGT", "A", "AC", "G"),
    ALT   = c("G", "A", "ATT", "GT", "T"),
    FILTER = c("PASS", "PASS", "PASS", "LowQual", "PASS"),
    Gene_Name = c("TP53", "KRAS", "EGFR", "APC", "PTEN"),
    Annotation = c("missense_variant", "frameshift_variant&splice_region_variant",
                   "frameshift_variant", "stop_gained", "5_prime_UTR_truncation&exon_loss_variant"),
    Annotation_Impact = c("MODERATE", "HIGH", "HIGH", "HIGH", "HIGH"),
    Transcript_BioType = "protein_coding", HGVS.c = "c.1A>G", HGVS.p = "p.M1V",
    stringsAsFactors = FALSE)
  if (tumor_only) {
    df$Sample <- c("S1", "S1", "S2", "S2", "S3")
    df$ADref <- c(90, 80, 50, 70, 60); df$ADalt <- c(10, 20, 50, 30, 4)
    df$DPtotal <- df$ADref + df$ADalt; df$VAF <- df$ADalt / df$DPtotal
  } else {
    df$Tumor <- c("S1", "S1", "S2", "S2", "S3"); df$Control <- c("N1", "N1", "N2", "N2", "N3")
    df$Tumor_Reference_Reads <- c(90, 80, 50, 70, 60); df$Tumor_Alternative_Reads <- c(10, 20, 50, 30, 4)
    df$Tumor_Total_Reads <- df$Tumor_Reference_Reads + df$Tumor_Alternative_Reads
    df$Tumor_VAF_Alternative <- df$Tumor_Alternative_Reads / df$Tumor_Total_Reads
    df$Control_Reference_Reads <- c(100, 100, 100, 99, NA); df$Control_Alternative_Reads <- c(0, 0, 0, 1, NA)
    df$Control_Total_Reads <- c(100, 100, 100, 100, NA); df$Control_VAF_Alternative <- c(0, 0, 0, 0.01, NA)
  }
  df$COSMIC99_Coding <- c("COSV1", NA, NA, NA, NA); df$COSMIC99_nonCoding <- NA_character_
  writexl::write_xlsx(df, path)
  path
}
