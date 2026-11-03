test_that("fromParse2MAF renames columns and classifies variants (tumor-only)", {
  f <- make_parse_xlsx(tempfile(fileext = ".xlsx"), tumor_only = TRUE)
  maf <- fromParse2MAF(f, tumor_only = TRUE)
  expect_true(all(c("Tumor_Sample_Barcode", "Hugo_Symbol", "VAF", "End_Position") %in% names(maf)))
  expect_false(any(c("dif_len", "consequence", "ref_len", "alt_len") %in% names(maf)))
  expect_equal(maf$Variant_Classification,
               c("Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"))
  expect_equal(maf$Variant_Type, c("SNP", "DEL", "INS", "DNP", "SNP"))
})

test_that("End_Position of a deletion uses the REF length", {
  f <- make_parse_xlsx(tempfile(fileext = ".xlsx"))
  maf <- fromParse2MAF(f, tumor_only = TRUE)
  expect_equal(maf$End_Position[2], 200 + 4 - 1)   # REF = ACGT, ALT = A
  expect_equal(maf$End_Position[1], 100)           # SNP
})

test_that("fromParse2MAF stops with a clear message when files or columns are missing", {
  expect_error(fromParse2MAF("no_such_file.xlsx"), "do not exist")
  f <- make_parse_xlsx(tempfile(fileext = ".xlsx"), tumor_only = TRUE)
  expect_error(fromParse2MAF(f, tumor_only = FALSE), "not found in the input")   # wrong mode
})

test_that("several files are combined", {
  f1 <- make_parse_xlsx(tempfile(fileext = ".xlsx")); f2 <- make_parse_xlsx(tempfile(fileext = ".xlsx"))
  expect_equal(nrow(fromParse2MAF(c(f1, f2), tumor_only = TRUE)), 10)
})
