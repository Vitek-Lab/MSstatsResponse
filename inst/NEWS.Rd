\name{MSstatsResponseNews}
\title{MSstatsResponse News}
\encoding{UTF-8}

\section{Version 1.4.0 (2026-11-25)}{
    \itemize{
        \item Added a protein turnover workflow: calculateQCScore, calculateConfidence and classifyTurnoverProteins.
        \item Added the simulated_turnover_data example dataset and a protein turnover vignette.
        \item Added tracer enrichment correction to calculateTurnoverRatios.
        \item Updated calculatePeptideWeights to score coverage relative to each protein's peptides.
        \item Observations with weight 0 no longer affect isotonic fits.
        \item Fixed isotonic fits that failed when pooled observations all had weight 0.
    }
}

\section{Changes in version 1.0.0 (2025-09-22)}{\itemize{
\item Submitted to Bioconductor
}}
