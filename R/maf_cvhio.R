# =====================================================================================
#  Alternative MAF for the cVHIO portal (cBioPortal)
#
#  This file contains:
#    - .snpeff_to_mut_type()              Annotation (SnpEff, 1st consequence) -> mut_type de la funcion Python
#    - .get_variant_classification()      R port of get_variant_classification (Python)
#    - .invalid_hgvsp_short_treatment()   R port of invalid_HGVSp_Short_treatment (Python)
#    - .get_cbp_driver()                  R port of get_cbp_driver_from_cbp_driver_annotation (Python)
#    - .cvhio_driver_annotation()         cbp_driver_annotation from CGI and OncoKB
#    - .read_hugo_entrez()                reads the Hugo_Symbol -> Entrez_Gene_Id JSON
#    - .make_maf_cvhio()                  builds the final table (called by fromParse2MAF)
#
#  The three ported functions have been validated against the original Python functions
#  (see tests/testthat/test-maf_cvhio.R). Only difference: 'synonymous SNV' -> 'missense'.
# =====================================================================================

# mut_type -> Variant_Classification mappings (copied exactly from get_variant_classification),
# with ONE requested addition: 'synonymous SNV' -> 'missense'.
.cvhio_vc_table <- list(
  inframe = c("frameshift deletion", "frameshift deletion?", "frameshift insertion",
                 "frameshift block substitution", "frameshift delins", "frameshift substitution",
                 "frameshift insertion / stoploss", "frameshift bloc substitution", "frameshift delation",
                 "frameshift mutation", "frameshift mutaiion"),
  missense = c("nonsynonymous SNV", "nonsynonymous substitution", "nonsynonymous insertion",
               "nonsynonymous mutation", "synonymous SNV"),
  truncating = c("startloss","stopgain", "stoploss", "upstream", "nonframeshift deletion", "nonframeshift substitution",
                 "nonframeshift insertion", "nonframeshift block substitution", "nonframeshift delins",
                 "nonframeshift mutation"),
  splice = c("splicing", "splice", "promoter")
)

# Annotation (SnpEff) -> mut_type, using the vocabulary understood by get_variant_classification.
# Applied to the FIRST consequence (the part before the first '&'). Rules are ordered: the first matching rule wins.
# The closest equivalent was chosen; terms without an equivalent (intron, UTR, intergenic...)
# are left unchanged; .get_variant_classification() keeps them as they are and issues a warning.
# Values that are not SnpEff terms (e.g. a column already containing 'frameshift deletion') pass through unchanged.
.cvhio_snpeff_rules <- data.frame(
  pattern = c("^splice_",                                              # splice_acceptor / donor / region ...
              "frameshift",                                            # frameshift_variant
              "inframe_deletion",                                      # inframe_deletion, disruptive_/conservative_inframe_deletion
              "inframe_insertion",
              "^stop_gained$",
              "^stop_lost$",
              "^start_lost$|^initiator_codon_variant$",                # loss of the start codon ~ 'stoploss'
              "^(synonymous_variant|stop_retained_variant|start_retained_variant)$",
              "^(missense_variant|coding_sequence_variant|rare_amino_acid_variant)$",
              "^upstream_gene_variant$",
              "^downstream_gene_variant$"),                            # ~ 'upstream'
  mut_type = c("splicing",
               "frameshift mutation",
               "nonframeshift deletion",
               "nonframeshift insertion",
               "stopgain",
               "stoploss",
               "startloss",
               "synonymous SNV",
               "nonsynonymous SNV",
               "upstream",
               "downstream"),
  stringsAsFactors = FALSE
)

#' @noRd
.snpeff_to_mut_type <- function(term) {
  term <- as.character(term)
  out  <- term
  done <- is.na(term)
  for (i in seq_len(nrow(.cvhio_snpeff_rules))) {
    hit <- !done & grepl(.cvhio_snpeff_rules$pattern[i], term)
    out[hit] <- .cvhio_snpeff_rules$mut_type[i]
    done <- done | hit
  }
  out
}

# mut_type -> Variant_Classification. Vectorized.
#  - Values in the table      -> their category (inframe / missense / truncating / splice)
#  - 'protein'                -> 'protein'
#  - 'non' o ''               -> ''
#  - any other value          -> returned unchanged (as in Python's else branch), with ONE warning
#  - NA                       -> NA
#' @noRd
.get_variant_classification <- function(mut_type) {
  mut_type <- as.character(mut_type)
  out <- mut_type
  for (cl in names(.cvhio_vc_table)) out[mut_type %in% .cvhio_vc_table[[cl]]] <- cl
  out[mut_type %in% c("non", "")] <- ""
  known    <- c(unlist(.cvhio_vc_table, use.names = FALSE), "protein", "non", "")
  unmapped <- setdiff(unique(mut_type[!is.na(mut_type)]), known)
  if (length(unmapped) > 0) {
    message("maf_cvhio: ", length(unmapped), " values of mut_type without equivalence: ",
            paste(utils::head(unmapped, 10), collapse = " | "), if (length(unmapped) > 10) " | ...")
  }
  out
}

# Invalid HGVSp_Short -> corrected value. Vectorized. Arguments as in Python:
#   hgvsp_short : HGVSp_Short original (HGVS.p). An empty/NA cell is treated as "".
#   mut_type    : the ALREADY calculated Variant_Classification (compared with "splice").
#   transcript  : string with c. notation (HGVS.c), e.g. "c.375+1G>A".
#  - hgvsp == "" y mut_type == "splice"   -> "p.X<round(pos_c/3)>_splice" (NA NA if there is no c.<number>)
#  - hgvsp == "p.?" y 'c.' en transcript  -> "p.M1" if it starts with "c.1-"; otherwise, NA
#  - otherwise                             -> hgvsp unchanged
# (In Python, the "NA" cases implicitly return None; here they are kept as NA.)
#' @noRd
.invalid_hgvsp_short_treatment <- function(hgvsp_short, mut_type, transcript) {
  hgvsp      <- as.character(hgvsp_short); hgvsp[is.na(hgvsp)] <- ""
  mut_type   <- as.character(mut_type)
  transcript <- as.character(transcript)
  out <- hgvsp
  
  splice_empty <- hgvsp == "" & !is.na(mut_type) & mut_type == "splice"
  has_c        <- !is.na(transcript) & grepl("c.", transcript, fixed = TRUE)
  p_unknown    <- !splice_empty & hgvsp == "p.?" & has_c
  
  for (i in which(splice_empty)) {
    m <- regmatches(transcript[i], regexec("c\\.(\\d+)", transcript[i]))[[1]]
    out[i] <- if (length(m) == 2) sprintf("p.X%.0f_splice", round(as.numeric(m[2]) / 3)) else NA_character_
  }
  for (i in which(p_unknown)) {
    c_name <- regmatches(transcript[i], regexpr("c\\.[^:]*", transcript[i]))
    out[i] <- if (startsWith(c_name, "c.1-")) "p.M1" else NA_character_
  }
  out
}

# cbp_driver_annotation -> cbp_driver. Exact match (as in Python).
#' @noRd
.get_cbp_driver <- function(cbp_driver_annotation) {
  ifelse(as.character(cbp_driver_annotation) %in% c("pathogenic / likely pathogenic", "vus"),
         "Putative_Driver", "Unknown")
}

# OncoKB labels considered pathogenic
.cvhio_oncokb_pathogenic <- c("Likely Oncogenic", "Oncogenic")
# Original Excel columns used to derive cbp_driver_annotation
.cvhio_cgi_col    <- "CGI-Oncogenicity summary"
.cvhio_oncokb_col <- "OncoKB"

# cbp_driver_annotation from CGI and OncoKB (sin 'benign'):
#   'pathogenic / likely pathogenic'  if CGI is oncogenic (same labels as .default_cgi_list in filterMAF)
#                                     or OncoKB is Oncogenic / Likely Oncogenic
#   'vus'                             if there is any CGI/OncoKB annotation but none is oncogenic
#                                     (Passenger, Likely Neutral, Inconclusive, Unknown...)
#   NA                                if the variant has no CGI/OncoKB annotation
# The comparison is case-insensitive and ignores leading/trailing whitespace.
#' @noRd
.cvhio_driver_annotation <- function(cgi = NULL, oncokb = NULL) {
  n <- max(length(cgi), length(oncokb))
  norm <- function(x) {
    if (is.null(x)) return(rep(NA_character_, n))
    x <- tolower(trimws(as.character(x))); x[x == ""] <- NA_character_; x
  }
  cgi <- norm(cgi); oncokb <- norm(oncokb)
  pathogenic <- (!is.na(cgi) & cgi %in% tolower(.default_cgi_list)) |
    (!is.na(oncokb) & oncokb %in% tolower(.cvhio_oncokb_pathogenic))
  annotated  <- !is.na(cgi) | !is.na(oncokb)
  ifelse(pathogenic, "pathogenic / likely pathogenic", ifelse(annotated, "vus", NA_character_))
}

# JSON {Hugo_Symbol: "Entrez"} -> named vector
#' @noRd
.read_hugo_entrez <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("maf_cvhio: no encuentro el JSON Hugo -> Entrez ('", path, "'). Indica su ruta con entrez_json.", call. = FALSE)
  }
  unlist(jsonlite::fromJSON(path), use.names = TRUE)
}

# Builds the alternative MAF.
#   maf  : output from fromParse2MAF (columns already renamed to MAF format)
#   raw  : original Excel file (same rows and order as `maf`), used to read columns not present in the MAF
#' @noRd
.make_maf_cvhio <- function(maf, raw, entrez, ncbi_build, mut_type_col) {
  stopifnot(nrow(maf) == nrow(raw))
  
  # Annotation -> first consequence (before '&') -> mut_type -> Variant_Classification
  first_consequence      <- sub("&.*", "", as.character(raw[[mut_type_col]]))
  mut_type               <- .snpeff_to_mut_type(first_consequence)
  variant_classification <- .get_variant_classification(mut_type)
  
  # CGI + OncoKB -> cbp_driver_annotation -> cbp_driver
  cgi    <- if (.cvhio_cgi_col    %in% names(raw)) raw[[.cvhio_cgi_col]]    else NULL
  oncokb <- if (.cvhio_oncokb_col %in% names(raw)) raw[[.cvhio_oncokb_col]] else NULL
  if (is.null(cgi) && is.null(oncokb)) {
    message("maf_cvhio: the excel displayed 0 columns '", .cvhio_cgi_col, "' ni '", .cvhio_oncokb_col,
            "'; cbp_driver_annotation = NA y cbp_driver = 'Unknown'.")
    cbp_driver_annotation <- rep(NA_character_, nrow(maf))
  } else {
    cbp_driver_annotation <- .cvhio_driver_annotation(cgi, oncokb)
    n_path <- sum(cbp_driver_annotation == "pathogenic / likely pathogenic", na.rm = TRUE)
    n_vus  <- sum(cbp_driver_annotation == "vus", na.rm = TRUE)
    message("maf_cvhio: cbp_driver_annotation from ", paste(c(if (!is.null(cgi)) "CGI", if (!is.null(oncokb)) "OncoKB"), collapse = " + "),
            ": ", n_path, " pathogenic / likely pathogenic, ", n_vus, " vus, ", sum(is.na(cbp_driver_annotation)), " without annotation.")
  }
  
  # Entrez from Hugo_Symbol
  entrez_id <- unname(entrez[maf$Hugo_Symbol])
  n_missing <- sum(is.na(entrez_id))
  if (n_missing > 0) {
    miss <- sort(unique(maf$Hugo_Symbol[is.na(entrez_id)]))
    message("maf_cvhio: ", n_missing, " variants with Hugo_Symbol but with no Entrez in the JSON (",
            length(miss), " gen(s): ", paste(utils::head(miss, 10), collapse = ", "),
            if (length(miss) > 10) ", ...", "). Entrez_Gene_Id = NA.")
  }
  
  hgvsp_short <- .invalid_hgvsp_short_treatment(maf$AAChange_refGene, variant_classification, maf$BNChange)
  hgvsp_short[!is.na(hgvsp_short) & hgvsp_short == ""] <- NA_character_   # empty cell in the Excel file
  
  data.frame(
    Hugo_Symbol            = maf$Hugo_Symbol,
    Entrez_Gene_Id         = as.numeric(entrez_id),
    NCBI_Build             = ncbi_build,
    Chromosome             = maf$Chromosome,
    Start_Position         = maf$Start_Position,
    End_Position           = maf$End_Position,
    Variant_Classification = variant_classification,
    Reference_Allele       = maf$Reference_Allele,
    Tumor_Seq_Allele2      = maf$Tumor_Seq_Allele2,
    Tumor_Sample_Barcode   = maf$Tumor_Sample_Barcode,
    HGVSp_Short            = hgvsp_short,
    t_ref_count            = maf$ReferenceReads_TumorSample,
    t_alt_count            = maf$AlternativeReads_TumorSample,
    cbp_driver             = .get_cbp_driver(cbp_driver_annotation),
    cbp_driver_annotation  = cbp_driver_annotation,
    Variant_Annotation      = maf$Variant_Classification,
    stringsAsFactors = FALSE
  )
}

