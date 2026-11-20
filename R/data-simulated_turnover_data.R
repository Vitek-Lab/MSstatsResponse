#' Simulated protein turnover feature-level dataset
#'
#' A simulated heavy-leucine (Leu6) labeling time course, formatted like the
#' `FeatureLevelData` element returned by `MSstats::dataProcess()`. It is
#' used in the protein turnover vignette to demonstrate the MSstatsResponse
#' turnover workflow. **These are not real measurements.** Protein
#' identifiers, peptide sequences, run names, and intensities are all
#' synthetic.
#'
#' @docType data
#' @name simulated_turnover_data
#' @keywords datasets
#'
#' @format A data frame with one row per protein, peptide, label, and
#'   timepoint, and the following columns:
#' \describe{
#'   \item{PROTEIN}{Factor. Simulated protein identifier
#'     (`Protein_001` to `Protein_200`).}
#'   \item{PEPTIDE}{Factor. Simulated peptide sequence with charge state
#'     suffix (e.g. `ACDEFLGHK_2`). Labeled peptides contain exactly one
#'     leucine; unlabeled peptides contain none.}
#'   \item{TRANSITION, FEATURE}{Factor. Feature identifiers in MSstats
#'     format.}
#'   \item{LABEL}{Factor. `"L"` (light, pre-existing), `"H"` (heavy, newly
#'     synthesized), or `NA` for leucine-free peptides, which cannot carry
#'     the label and serve as normalization standards.}
#'   \item{GROUP}{Factor. Timepoint after labeling started: `1hr`, `2hrs`,
#'     `6hrs`, `12hrs`, `24hrs`, `48hrs`, `100hrs`, `170hrs`.}
#'   \item{RUN, SUBJECT, FRACTION, originalRUN}{Run, subject, fraction, and
#'     original run identifiers. One run per timepoint.}
#'   \item{censored}{Logical. Always `FALSE`.}
#'   \item{INTENSITY}{Numeric. Feature intensity on the original scale.}
#'   \item{ABUNDANCE, newABUNDANCE}{Numeric. log2 of `INTENSITY`.}
#'   \item{predicted}{Numeric. Always `NA`.}
#' }
#'
#' @details
#' The dataset contains 200 proteins drawn from four turnover profiles: 30
#' fast-turnover, 70 medium-lived, 20 long-lived, and 80 with no heavy
#' signal. For each labeled peptide, the tracer-corrected heavy fraction
#' follows first-order synthesis kinetics,
#' \eqn{f(t) = A (1 - e^{-k t})}, where the plateau \eqn{A} and rate
#' \eqn{k} are drawn per protein according to its profile. The observed
#' heavy share is \eqn{f(t)} scaled by a timepoint-specific tracer
#' enrichment (0.11, 0.10, 0.11, 0.12, 0.12, 0.13, 0.15, 0.17), plus a small
#' heavy-channel background.
#'
#' Peptide counts per protein, intensity levels, intensity-dependent
#' detection in each channel, and the scatter of heavy/light ratios at each
#' timepoint (including occasional outliers) were chosen to resemble a real
#' in vivo labeling experiment. The tracer enrichments are passed as
#' `tracer_constants` to [calculateTurnoverRatios()].
#'
#' @source Simulated.
#'
#' @seealso [calculateTurnoverRatios()], [calculateQCScore()],
#'   [classifyTurnoverProteins()]
#'
#' @examples
#' data("simulated_turnover_data", package = "MSstatsResponse")
#' head(simulated_turnover_data)
#' table(simulated_turnover_data$GROUP, simulated_turnover_data$LABEL)
"simulated_turnover_data"
