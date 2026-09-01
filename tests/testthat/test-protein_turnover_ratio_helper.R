# Tests for calculateTurnoverRatios and its helpers.

make_feature_data <- function() {
  grid <- expand.grid(
    PROTEIN = c("P1", "P2"),
    PEPTIDE = c("PEP[+80]TIDE", "SEQVENCE", "THIRDPEP"),
    GROUP = c("0hr", "1hr", "12hrs", "168hrs"),
    RUN = c("R1", "R2"),
    LABEL = c("H", "L"),
    CHARGE = c(2, 3),
    stringsAsFactors = FALSE
  )
  grid <- grid[!(grid$PROTEIN == "P2" & grid$PEPTIDE != "THIRDPEP"), ]
  grid <- grid[!(grid$PROTEIN == "P1" & grid$PEPTIDE == "THIRDPEP"), ]
  grid$INTENSITY <- seq_len(nrow(grid)) * 100
  rownames(grid) <- NULL
  grid
}

test_that("parse_timepoint converts hours, days and weeks to hours", {
  expect_equal(parse_timepoint(c("0hr", "1hr", "12hrs", "168hrs")),
               c(0, 1, 12, 168))
  expect_equal(parse_timepoint(c("2d", "1day")), c(48, 24))
  expect_equal(parse_timepoint(c("1w", "2weeks")), c(168, 336))
  expect_true(is.na(parse_timepoint("control")))
})

test_that("calculateTurnoverRatios errors on missing required columns", {
  feat <- make_feature_data()
  expect_error(
    calculateTurnoverRatios(feat[, setdiff(names(feat), "RUN")]),
    "Missing required columns: RUN"
  )
})

test_that("calculateTurnoverRatios returns the documented columns", {
  res <- calculateTurnoverRatios(make_feature_data())
  expect_true(data.table::is.data.table(res))
  expect_equal(names(res),
               c("Protein", "BaseSequence", "TimeVal", "Run",
                 "Heavy", "Light", "Total", "H_frac", "L_frac"))
  expect_gt(nrow(res), 0)
})

test_that("modifications are stripped from the peptide sequence", {
  res <- calculateTurnoverRatios(make_feature_data())
  expect_setequal(unique(res$BaseSequence),
                  c("PEPTIDE", "SEQVENCE", "THIRDPEP"))
  expect_false(any(grepl("[", res$BaseSequence, fixed = TRUE)))
})

test_that("fractions sum to one and Total is the H + L sum", {
  res <- calculateTurnoverRatios(make_feature_data())
  expect_equal(res$Total, res$Heavy + res$Light)
  expect_equal(res$H_frac + res$L_frac, rep(1, nrow(res)))
  expect_equal(res$H_frac, res$Heavy / res$Total)
})

test_that("duplicate features are collapsed by agg_function", {
  feat <- make_feature_data()
  n_expected <- nrow(unique(feat[, c("PROTEIN", "PEPTIDE", "GROUP", "RUN")]))
  res_max <- calculateTurnoverRatios(feat)
  expect_equal(nrow(res_max), n_expected)

  res_min <- calculateTurnoverRatios(feat, agg_function = min)
  expect_true(all(res_max$Heavy >= res_min$Heavy))

  cell <- feat[feat$PROTEIN == "P1" & feat$PEPTIDE == "SEQVENCE" &
                 feat$GROUP == "1hr" & feat$RUN == "R1" & feat$LABEL == "H", ]
  target <- res_max[res_max$Protein == "P1" & res_max$BaseSequence == "SEQVENCE" &
                      res_max$TimeVal == 1 & res_max$Run == "R1", ]
  expect_equal(target$Heavy, max(cell$INTENSITY))
})

test_that("unparseable timepoints and non-H/L labels are dropped", {
  feat <- make_feature_data()
  noise <- rbind(
    transform(feat[1:4, ], GROUP = "control"),
    transform(feat[1:4, ], LABEL = "M")
  )
  res_clean <- calculateTurnoverRatios(feat)
  res_noisy <- calculateTurnoverRatios(rbind(feat, noise))
  expect_equal(as.data.frame(res_noisy), as.data.frame(res_clean))
})

test_that("rows without both channels are dropped", {
  feat <- make_feature_data()
  drop <- feat$PROTEIN == "P1" & feat$PEPTIDE == "SEQVENCE" &
    feat$GROUP == "1hr" & feat$RUN == "R2" & feat$LABEL == "L"
  res <- calculateTurnoverRatios(feat[!drop, ])
  orphan <- res[res$Protein == "P1" & res$BaseSequence == "SEQVENCE" &
                  res$TimeVal == 1 & res$Run == "R2", ]
  expect_equal(nrow(orphan), 0)
  expect_false(anyNA(res$Heavy))
  expect_false(anyNA(res$Light))
})

test_that("no matching Heavy/Light labels warns and returns no rows", {
  feat <- transform(make_feature_data(), LABEL = "X")
  expect_warning(res <- calculateTurnoverRatios(feat),
                 "No data found matching Heavy/Light labels")
  expect_equal(nrow(res), 0)
})

test_that("peptide_selector is applied within each protein", {
  feat <- make_feature_data()
  keep_first <- function(df) {
    df[df$BaseSequence == sort(unique(df$BaseSequence))[1], ]
  }
  res <- calculateTurnoverRatios(feat, peptide_selector = keep_first)
  by_protein <- tapply(res$BaseSequence, res$Protein, function(x) unique(x))
  expect_equal(unname(by_protein[["P1"]]), "PEPTIDE")
  expect_equal(unname(by_protein[["P2"]]), "THIRDPEP")
})

test_that("tracer normalization rescales H_frac and floors L_frac at zero", {
  feat <- make_feature_data()
  tc <- c("0hr" = 1.0, "1hr" = 0.95, "12hrs" = 0.85, "168hrs" = 0.75)
  plain <- calculateTurnoverRatios(feat)
  norm <- calculateTurnoverRatios(feat, normalize_tracer = TRUE,
                                  tracer_constants = tc)

  expect_true("tracer_factor" %in% names(norm))
  expect_equal(norm$tracer_factor, unname(tc[match(norm$TimeVal, c(0, 1, 12, 168))]))
  expect_equal(norm$H_frac, plain$H_frac / norm$tracer_factor)
  expect_equal(norm$L_frac, pmax(1 - norm$H_frac, 0))
  expect_true(all(norm$L_frac >= 0))
})

test_that("tracer normalization requires constants", {
  expect_error(
    calculateTurnoverRatios(make_feature_data(), normalize_tracer = TRUE),
    "tracer_constants was not provided"
  )
})

test_that("custom column names and heavy/light labels are honored", {
  feat <- make_feature_data()
  names(feat)[names(feat) == "PROTEIN"] <- "Protein"
  names(feat)[names(feat) == "GROUP"] <- "Condition"
  feat$LABEL <- ifelse(feat$LABEL == "H", "Heavy", "Light")
  res <- calculateTurnoverRatios(
    feat, time_col = "Condition", peptide_col = "Protein",
    protein_col = "Protein", heavy_label = "Heavy", light_label = "Light"
  )
  expect_true(all(c("Heavy", "Light") %in% names(res)))
  expect_setequal(unique(res$BaseSequence), c("P1", "P2"))
})

test_that("a data.table input is not modified in place", {
  feat <- data.table::as.data.table(make_feature_data())
  before <- data.table::copy(feat)
  calculateTurnoverRatios(feat)
  expect_equal(feat, before)
})

test_that("data.frame and data.table inputs give the same result", {
  feat <- make_feature_data()
  expect_equal(calculateTurnoverRatios(feat),
               calculateTurnoverRatios(data.table::as.data.table(feat)))
})

test_that("kendall_monotonicity floors at zero and tolerates sparse input", {
  expect_equal(kendall_monotonicity(c(1, 2, 3), c(1, 2, 3)), 1)
  expect_equal(kendall_monotonicity(c(1, 2, 3), c(3, 2, 1)), 0)
  expect_equal(kendall_monotonicity(1, 1), 0)
  expect_equal(kendall_monotonicity(c(1, 2), c(NA, NA)), 0)
  expect_equal(kendall_monotonicity(c(1, 2, 3), c(5, 5, 5)), 0)
})

test_that("calculatePeptideWeights consumes calculateTurnoverRatios output", {
  feat <- make_feature_data()
  res <- calculateTurnoverRatios(feat[feat$RUN == "R1", ])
  weighted <- calculatePeptideWeights(res)
  expect_true(all(c("coverage_score", "monotonicity_score", "validity_flag",
                    "weight") %in% names(weighted)))
  expect_equal(nrow(weighted), nrow(res))
  expect_true(all(weighted$weight >= 0 & weighted$weight <= 1))
})
