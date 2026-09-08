# PhotoSeq analysis inputs

All count matrices are comma-separated tables with gene identifiers in the first column and non-negative integer UMI-deduplicated counts in the remaining columns.

`PhotoSeq_Spatial_GeneCounts.csv` and `PhotoSeq_Temporal_GeneCounts.csv` are the outputs of the upstream Spatial and Temporal count-generation scripts. Their source headers are mapped to the standardized sample IDs in `PhotoSeq_sample_metadata.csv` by `scripts/00_prepare_analysis_data.R`.

`PhotoSeq_TechnicalReplicate_GeneCounts.csv` is the output of the upstream technical-replicate count-generation script. Its two sample columns are `Breast_early_N1` and `Breast_early_N1_techrep2`.
