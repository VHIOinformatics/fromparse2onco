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

# Excel de demostracion para maf_cvhio (tumor-only). Datos SINTETICOS.
# Variant_Classification sale de 'Annotation' y cbp_driver_annotation de CGI + OncoKB.
make_cvhio_xlsx <- function(path) {
  df <- data.frame(
    CHROM = c("chr17", "chr12", "chr3", "chr5", "chr10", "chr1", "chr7", "chr7", "chr9", "chr2"),
    POS   = c(7675088, 25245350, 179218303, 112839942, 87933147, 26729700, 55019017, 140753336, 136496092, 29000000),
    REF   = c("C", "C", "A", "AC", "C", "G", "A", "A", "ACGT", "G"),
    ALT   = c("T", "T", "G", "A", "T", "A", "G", "T", "A", "A"),
    FILTER = "PASS",
    Gene_Name = c("TP53", "KRAS", "PIK3CA", "APC", "PTEN", "ARID1A", "EGFR", "BRAF", "NOTCH1", "FAKEGENE1"),
    Annotation = c("missense_variant", "missense_variant", "missense_variant", "frameshift_variant&splice_region_variant",
                   "stop_gained", "splice_donor_variant&intron_variant", "start_lost", "missense_variant",
                   "inframe_deletion", "synonymous_variant"),
    Annotation_Impact = c("MODERATE", "MODERATE", "MODERATE", "HIGH", "HIGH", "HIGH", "HIGH", "MODERATE", "MODERATE", "LOW"),
    Transcript_BioType = "protein_coding",
    HGVS.c = c("c.524G>A", "c.35G>A", "c.3140A>G", "c.4348delC", "c.697C>T", "c.4000+1G>A", "c.1-2A>G", "c.1799T>A", "c.4800_4802del", "c.300G>A"),
    HGVS.p = c("p.R175H", "p.G12D", "p.H1047R", "p.R1450fs", "p.R233*", NA, "p.?", "p.V600E", "p.L1601del", "p.K100="),
    Sample = c("PT01", "PT01", "PT02", "PT02", "PT03", "PT03", "PT04", "PT04", "PT05", "PT05"),
    ADref = c(80, 70, 150, 60, 90, 100, 120, 40, 85, 200), ADalt = c(70, 30, 50, 40, 60, 50, 30, 60, 15, 20),
    stringsAsFactors = FALSE)
  df$DPtotal <- df$ADref + df$ADalt; df$VAF <- round(df$ADalt / df$DPtotal, 4)
  df[["CGI-Oncogenicity summary"]] <- c("Oncogenic (predicted by BoostDM)", "Oncogenic (predicted by RulesDM)",
                                        "Passenger (predicted by BoostDM)", "Potentially Oncogenic (predicted by RulesDM)",
                                        "Passenger (predicted by RulesDM)", NA, "Unknown", NA, "Passenger (predicted by BoostDM)",
                                        "Passenger (predicted by BoostDM)")
  df[["CGI-Oncogenicity prediction"]] <- "x"
  df$OncoKB <- c("Oncogenic", "Likely Oncogenic", "Likely Neutral", "Inconclusive", "Oncogenic", "Unknown", NA, NA, NA, "Likely Neutral")
  writexl::write_xlsx(df, path)
  path
}
