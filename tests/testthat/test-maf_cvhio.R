# Valores esperados GENERADOS EJECUTANDO las funciones originales en Python (get_variant_classification,
# invalid_HGVSp_Short_treatment y get_cbp_driver_from_cbp_driver_annotation).

test_that("get_variant_classification coincide con Python (salvo synonymous SNV -> missense)", {
  mut      <- c("frameshift deletion", "frameshift deletion?", "frameshift insertion", "frameshift block substitution", "frameshift delins", "frameshift substitution", "frameshift insertion / stoploss", "frameshift bloc substitution", "frameshift delation", "frameshift mutation", "frameshift mutaiion", "nonsynonymous SNV", "nonsynonymous substitution", "nonsynonymous insertion", "nonsynonymous mutation", "stopgain", "stoploss", "upstream", "nonframeshift deletion", "nonframeshift substitution", "nonframeshift insertion", "nonframeshift block substitution", "nonframeshift delins", "nonframeshift mutation", "splicing", "splice", "promoter", "protein", "non", "", "synonymous SNV", "missense_variant", "Nonsynonymous SNV", "stopgain ")
  expected <- c("inframe", "inframe", "inframe", "inframe", "inframe", "inframe", "inframe", "inframe", "inframe", "inframe", "inframe", "missense", "missense", "missense", "missense", "truncating", "truncating", "truncating", "truncating", "truncating", "truncating", "truncating", "truncating", "truncating", "splice", "splice", "splice", "protein", "", "", "synonymous SNV", "missense_variant", "Nonsynonymous SNV", "stopgain ")
  expected[mut == "synonymous SNV"] <- "missense"   # unica ampliacion respecto a Python (pedida)
  expect_equal(suppressMessages(.get_variant_classification(mut)), expected)
  expect_message(.get_variant_classification(c("intron_variant", "stopgain")), "sin equivalencia")
  expect_equal(suppressMessages(.get_variant_classification(NA_character_)), NA_character_)
})

test_that("invalid_HGVSp_Short_treatment coincide con la funcion original de Python", {
  hgvsp      <- c("", "", "", "", "", "", "", "p.?", "p.?", "p.?", "p.?", "p.?", "p.?", "p.R175H", "p.R175H", "", "", "p.Met1?", "p.?")
  mut_type   <- c("splice", "splice", "splice", "splice", "splice", "splice", "splice", "missense", "missense", "missense", "missense", "missense", "splice", "missense", "splice", "missense", "truncating", "missense", "")
  transcript <- c("c.375+1G>A", "NM_000546:exon5:c.376-2A>G", "c.-14G>A", "c.1A>G", "c.100000+1G>A", "c.2G>A", "c.4G>A", "c.1-1G>A", "NM_1:c.1-15C>T", "c.5A>G", "c.-5A>G", "g.100A>G", "c.1-3G>T", "c.524G>A", "c.524G>A", "c.5A>G", "c.5A>G", "c.1A>G", "c.1-1G>A")
  expected   <- c("p.X125_splice", "p.X125_splice", NA_character_, "p.X0_splice", "p.X33333_splice", "p.X1_splice", "p.X1_splice", "p.M1", "p.M1", NA_character_, NA_character_, "p.?", "p.M1", "p.R175H", "p.R175H", "", "", "p.Met1?", "p.M1")
  expect_equal(.invalid_hgvsp_short_treatment(hgvsp, mut_type, transcript), expected)
  # Celda vacia del Excel (NA) en un splice: se trata como ""
  expect_equal(.invalid_hgvsp_short_treatment(NA, "splice", "c.375+1G>A"), "p.X125_splice")
})

test_that("get_cbp_driver coincide con la funcion original de Python", {
  ann      <- c("pathogenic / likely pathogenic", "vus", "VUS", "Pathogenic / Likely pathogenic", "benign", "likely benign", "", "other")
  expected <- c("Putative_Driver", "Putative_Driver", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown")
  expect_equal(.get_cbp_driver(ann), expected)
  expect_equal(.get_cbp_driver(NA_character_), "Unknown")
})

test_that("Annotation (SnpEff) se traduce al mut_type mas parecido y luego a la categoria", {
  snpeff <- c("missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
              "inframe_deletion", "disruptive_inframe_insertion", "splice_acceptor_variant", "splice_region_variant",
              "upstream_gene_variant", "downstream_gene_variant", "intron_variant", "3_prime_UTR_variant", "nonsynonymous SNV")
  vc <- suppressMessages(.get_variant_classification(.snpeff_to_mut_type(snpeff)))
  expect_equal(vc, c("missense", "missense", "truncating", "truncating", "truncating", "inframe",
                     "truncating", "truncating", "splice", "splice",
                     "truncating", "truncating", "intron_variant", "3_prime_UTR_variant", "missense"))
  expect_true(is.na(.snpeff_to_mut_type(NA_character_)))
})

test_that("cbp_driver_annotation sale de CGI y OncoKB: pathogenic / vus / NA, sin benign", {
  cgi    <- c("Oncogenic (predicted by BoostDM)", "Potentially Oncogenic (predicted by RulesDM)", "Passenger (predicted by BoostDM)",
              "Passenger (predicted by RulesDM)", NA, "Unknown", NA, "")
  oncokb <- c(NA, "Inconclusive", "Likely Oncogenic", "Likely Neutral", "Oncogenic", NA, NA, "")
  expect_equal(.cvhio_driver_annotation(cgi, oncokb),
               c("pathogenic / likely pathogenic", "pathogenic / likely pathogenic", "pathogenic / likely pathogenic",
                 "vus", "pathogenic / likely pathogenic", "vus", NA, NA))
  expect_equal(.cvhio_driver_annotation(cgi = "ONCOGENIC (predicted by boostdm) "), "pathogenic / likely pathogenic")
  expect_equal(.cvhio_driver_annotation(oncokb = c("Oncogenic", "Likely Neutral")), c("pathogenic / likely pathogenic", "vus"))
  expect_false(any(.cvhio_driver_annotation(cgi, oncokb) %in% c("benign", "likely benign")))
})

cvhio_cols <- c("Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome", "Start_Position", "End_Position",
                "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
                "HGVSp_Short", "t_ref_count", "t_alt_count", "cbp_driver", "cbp_driver_annotation")

test_that("maf_cvhio = TRUE escribe el Excel con las 15 columnas, en orden, y no altera el MAF devuelto", {
  f <- make_cvhio_xlsx(tempfile(fileext = ".xlsx")); out <- tempfile(fileext = ".xlsx")
  maf_off <- fromParse2MAF(f, tumor_only = TRUE)
  msgs <- testthat::capture_messages(
    maf_on <- fromParse2MAF(f, tumor_only = TRUE, maf_cvhio = TRUE, cvhio_file = out, ncbi_build = "GRCh37", cvhio_preview = 0))
  expect_true(any(grepl("10 variantes escritas", msgs)))
  expect_true(any(grepl("CGI \\+ OncoKB: 4 pathogenic / likely pathogenic, 5 vus, 1 sin anotacion", msgs)))
  expect_equal(maf_on, maf_off)                       # el oncoplot recibe exactamente lo mismo
  x <- as.data.frame(readxl::read_excel(out))
  expect_equal(names(x), cvhio_cols)
  expect_equal(nrow(x), 10)
  expect_true(all(x$NCBI_Build == "GRCh37"))
  expect_equal(x$Entrez_Gene_Id[x$Hugo_Symbol == "TP53"], 7157)
  expect_equal(x$Entrez_Gene_Id[x$Hugo_Symbol == "KRAS"], 3845)
  expect_true(is.na(x$Entrez_Gene_Id[x$Hugo_Symbol == "FAKEGENE1"]))
  # 1a consecuencia (antes del &) -> categoria. APC: frameshift&splice_region -> inframe; FAKEGENE1: synonymous -> missense
  expect_equal(x$Variant_Classification, c("missense", "missense", "missense", "inframe", "truncating", "splice", "truncating", "missense", "truncating", "missense"))
  expect_equal(x$t_ref_count, c(80, 70, 150, 60, 90, 100, 120, 40, 85, 200))
  expect_equal(x$t_alt_count, c(70, 30, 50, 40, 60, 50, 30, 60, 15, 20))
  expect_equal(x$HGVSp_Short[c(1, 6, 7)], c("p.R175H", "p.X1333_splice", "p.M1"))
  expect_equal(x$cbp_driver_annotation, c("pathogenic / likely pathogenic", "pathogenic / likely pathogenic", "vus", "pathogenic / likely pathogenic",
                                          "pathogenic / likely pathogenic", "vus", "vus", NA, "vus", "vus"))
  expect_equal(x$cbp_driver, c(rep("Putative_Driver", 7), "Unknown", "Putative_Driver", "Putative_Driver"))
  expect_false(any(c("benign", "likely benign") %in% x$cbp_driver_annotation))
})

test_that("maf_cvhio avisa de genes sin Entrez y de la falta de CGI/OncoKB", {
  f <- make_cvhio_xlsx(tempfile(fileext = ".xlsx")); out <- tempfile(fileext = ".xlsx")
  msgs <- testthat::capture_messages(
    fromParse2MAF(f, tumor_only = TRUE, maf_cvhio = TRUE, cvhio_file = out, cvhio_preview = 0))
  expect_true(any(grepl("sin Entrez.*FAKEGENE1", msgs)))
  d <- as.data.frame(readxl::read_excel(f)); d <- d[, !names(d) %in% c("CGI-Oncogenicity summary", "CGI-Oncogenicity prediction", "OncoKB")]
  f2 <- tempfile(fileext = ".xlsx"); writexl::write_xlsx(d, f2)
  msgs2 <- testthat::capture_messages(
    fromParse2MAF(f2, tumor_only = TRUE, maf_cvhio = TRUE, cvhio_file = out, cvhio_preview = 0))
  expect_true(any(grepl("no tiene columnas", msgs2)))
  x <- as.data.frame(readxl::read_excel(out))
  expect_true(all(is.na(x$cbp_driver_annotation)))
  expect_true(all(x$cbp_driver == "Unknown"))
})

test_that("maf_cvhio da errores claros", {
  f <- make_cvhio_xlsx(tempfile(fileext = ".xlsx"))
  expect_error(fromParse2MAF(f, tumor_only = TRUE, maf_cvhio = TRUE, cvhio_mut_type_col = "no_existe"), "no esta en el Excel")
  expect_error(fromParse2MAF(f, tumor_only = TRUE, maf_cvhio = TRUE, entrez_json = "no_existe.json"), "JSON")
})

test_that("sin maf_cvhio no se escribe nada", {
  f <- make_cvhio_xlsx(tempfile(fileext = ".xlsx")); out <- tempfile(fileext = ".xlsx")
  fromParse2MAF(f, tumor_only = TRUE, cvhio_file = out)
  expect_false(file.exists(out))
})
