# Frozen PPI reproducibility snapshot

This directory records the compact inputs and expected tables for the PhotoSeq
protein-protein interaction analysis. The STRING responses were obtained on
2 September 2026 for *Mus musculus* (NCBI taxonomy identifier 10090) with a
minimum required interaction score of 400, functional network type and no
additional interactors. The returned interaction scores are represented on a
0-1 scale; the reported network retains interactions with a score greater than
0.4.

The default workflow reads the frozen STRING response so that database
updates do not change the reported result. To submit the current candidate list
to the live STRING API, follow the optional command in the downstream README.

## Files

| File | Role |
| --- | --- |
| `candidate_genes_five_way_intersection.csv` | Fifty genes significant in all five specified DESeq2 contrasts |
| `deg_summary.csv` | Numbers of tested, upregulated and downregulated genes in the five contrasts |
| `string_mapping.tsv` | Raw STRING identifier-mapping response |
| `string_network_score400.tsv` | Raw STRING network response |
| `ppi_edges.csv` | Sixteen retained interactions among fourteen non-isolated nodes |
| `ppi_node_scores.csv` | PPI strength, Degree, MCC, MCODE membership and metric percentiles for all fifty genes |
| `ppi_mcode_modules.csv` | MCODE module membership, node count, internal-edge count and score |
| `ppi_consensus_candidate_genes.csv` | Genes shared by the MCC top five, Degree top five and MCODE module |
| `ppi_unmatched_genes.csv` | Six input genes absent from the STRING mapping response |
| `ppi_run_summary.csv` | Input, mapping, node, edge and module summary |
| `ppi_figure_edges.csv` | Edge source data used by the four-panel figure |
| `ppi_figure_nodes.csv` | Node coordinates and highlighting used by the four-panel figure |
| `ppi_figure_highlights.csv` | Highlighted and non-highlighted genes for each panel |

The frozen network contains fourteen non-isolated nodes and sixteen
edges. MCC and Degree rank the same five genes represented by the only detected
MCODE module: `Actg1`, `Itga3`, `Itga8`, `Itga9` and `Lama4`. The module contains
five nodes and nine internal edges and has an MCODE score of 4.5. In the MCODE
panel, these module genes are highlighted within the complete PPI network; the
panel is intentionally not restricted to the isolated module.

The percentile-based `combined_score` in `ppi_node_scores.csv` provides a
traceable ordering of the node table. Candidate-gene identification is based
on agreement among the MCC top five, Degree top five and MCODE module, not on a
separate combined-score threshold.
