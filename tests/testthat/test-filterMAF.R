test_that("filterMAF applies VAF / reads / FILTER thresholds (tumor-only)", {
  maf <- fromParse2MAF(make_parse_xlsx(tempfile(fileext = ".xlsx")), tumor_only = TRUE)
  out <- filterMAF(maf, tumor_only = TRUE, VAF_tumor = 0.1, alt_tumor_reads = 5)
  expect_equal(out$Hugo_Symbol, c("TP53", "KRAS", "EGFR"))   # APC is LowQual, PTEN has VAF 0.06
  out_out <- filterMAF(maf, tumor_only = TRUE, tumor_samples_out = "S1")
  expect_false(any(out_out$Tumor_Sample_Barcode == "S1"))
})

test_that("cosmic filter works in BOTH modes (it failed in paired mode)", {
  for (mode in c(TRUE, FALSE)) {
    maf <- fromParse2MAF(make_parse_xlsx(tempfile(fileext = ".xlsx"), tumor_only = mode),
                         tumor_only = mode, cosmic = TRUE)
    out <- suppressMessages(filterMAF(maf, tumor_only = mode, cosmic = TRUE))
    expect_equal(out$Hugo_Symbol, "TP53")
  }
})

test_that("paired mode removes variants with VAF in the control above the limit and reports NAs", {
  maf <- fromParse2MAF(make_parse_xlsx(tempfile(fileext = ".xlsx"), tumor_only = FALSE), tumor_only = FALSE)
  expect_message(out <- filterMAF(maf, tumor_only = FALSE, VAF_control = 0), "NA in 'VAF_normal'")
  expect_false("APC" %in% out$Hugo_Symbol)   # control VAF 0.01 > 0
  expect_false("PTEN" %in% out$Hugo_Symbol)  # control VAF is NA -> removed
})
