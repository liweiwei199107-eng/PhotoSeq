"""Build and score the PhotoSeq protein-protein interaction network.

By default, this script uses the frozen STRING responses in ``reference/ppi``
so that the publication result is exactly reproducible.  Pass
``--refresh-string`` to submit the current five-way DEG intersection to the
STRING API and save a new response in the generated output directory.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path
import urllib.parse
import urllib.request


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = PROJECT_ROOT / "reference" / "ppi"
INTERSECTION_PATH = PROJECT_ROOT / "12_Venn_All" / "02.Venn_DEGs.csv"
OUTPUT_DIR = PROJECT_ROOT / "12_Venn_All" / "04.PPI" / "network"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

STRING_SPECIES = 10090
STRING_REQUIRED_SCORE = 400
INTERACTION_THRESHOLD = 0.4
MCODE_DEGREE_CUTOFF = 2
MCODE_NODE_SCORE_CUTOFF = 0.2
MCODE_K_CORE = 2
MCODE_MAX_DEPTH = 100


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh-string",
        action="store_true",
        help="Query the live STRING API instead of using the frozen response.",
    )
    return parser.parse_args()


def read_genes(allow_snapshot_mismatch: bool) -> list[str]:
    path = (
        INTERSECTION_PATH
        if INTERSECTION_PATH.exists()
        else REFERENCE_DIR / "candidate_genes_five_way_intersection.csv"
    )
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        gene_column = "gene" if "gene" in (reader.fieldnames or []) else "x"
        if gene_column not in (reader.fieldnames or []):
            raise ValueError(f"No gene column found in {path}")
        genes = [row[gene_column].strip() for row in reader if row[gene_column].strip()]
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError(f"No genes found in {path}")
    frozen_path = REFERENCE_DIR / "candidate_genes_five_way_intersection.csv"
    with frozen_path.open(encoding="utf-8-sig", newline="") as handle:
        frozen = sorted(row["gene"].strip() for row in csv.DictReader(handle))
    if sorted(genes) != frozen and not allow_snapshot_mismatch:
        raise ValueError(
            "Generated candidate genes differ from the frozen STRING snapshot; "
            "use --refresh-string to obtain STRING responses for this gene set"
        )
    return genes


def string_post(endpoint: str, params: dict[str, str]) -> str:
    payload = urllib.parse.urlencode(params).encode()
    request = urllib.request.Request(
        f"https://string-db.org/api/tsv/{endpoint}",
        data=payload,
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8")


def obtain_string_responses(genes: list[str], refresh: bool) -> tuple[str, str, str]:
    if refresh:
        parameters = {"identifiers": "\r".join(genes), "species": str(STRING_SPECIES)}
        mapping_raw = string_post("get_string_ids", parameters)
        network_raw = string_post(
            "network",
            {
                **parameters,
                "required_score": str(STRING_REQUIRED_SCORE),
                "network_type": "functional",
                "add_nodes": "0",
            },
        )
        (OUTPUT_DIR / "string_mapping.tsv").write_text(mapping_raw, encoding="utf-8")
        (OUTPUT_DIR / "string_network_score400.tsv").write_text(network_raw, encoding="utf-8")
        return mapping_raw, network_raw, "live STRING API"

    mapping_raw = (REFERENCE_DIR / "string_mapping.tsv").read_text(encoding="utf-8")
    network_raw = (REFERENCE_DIR / "string_network_score400.tsv").read_text(encoding="utf-8")
    return mapping_raw, network_raw, "frozen publication snapshot"


def parse_edges(
    genes: list[str], mapping_raw: str, network_raw: str
) -> tuple[list[tuple[str, str, float]], set[str]]:
    mapped_rows = list(csv.DictReader(io.StringIO(mapping_raw), delimiter="\t"))
    gene_by_lower = {gene.lower(): gene for gene in genes}
    string_id_to_gene: dict[str, str] = {}
    for row in mapped_rows:
        preferred_name = row.get("preferredName", "")
        if preferred_name.lower() in gene_by_lower:
            string_id_to_gene[row["stringId"]] = gene_by_lower[preferred_name.lower()]

    network_rows = list(csv.DictReader(io.StringIO(network_raw), delimiter="\t"))
    edges: list[tuple[str, str, float]] = []
    seen_pairs: set[tuple[str, str]] = set()
    for row in network_rows:
        gene_a = string_id_to_gene.get(row.get("stringId_A", ""))
        gene_b = string_id_to_gene.get(row.get("stringId_B", ""))
        if not gene_a or not gene_b or gene_a == gene_b:
            continue
        score = float(row["score"])
        if score <= INTERACTION_THRESHOLD:
            continue
        pair = tuple(sorted((gene_a, gene_b)))
        if pair in seen_pairs:
            raise ValueError(f"Duplicate undirected STRING edge: {pair}")
        seen_pairs.add(pair)
        edges.append((gene_a, gene_b, score))

    mapped_genes = set(string_id_to_gene.values())
    return edges, mapped_genes


def core_numbers(nodes: set[str], adjacency: dict[str, set[str]]) -> dict[str, int]:
    remaining = set(nodes)
    degree = {node: len(adjacency[node] & remaining) for node in remaining}
    result: dict[str, int] = {}
    current = 0
    while remaining:
        node = min(remaining, key=lambda item: (degree[item], item))
        current = max(current, degree[node])
        result[node] = current
        remaining.remove(node)
        for neighbor in adjacency[node] & remaining:
            degree[neighbor] -= 1
    return result


def maximal_cliques(nodes: list[str], adjacency: dict[str, set[str]]) -> list[list[str]]:
    cliques: list[list[str]] = []

    def bron_kerbosch(r: list[str], p: set[str], x: set[str]) -> None:
        if not p and not x:
            cliques.append(list(r))
            return
        union = p | x
        pivot = max(union, key=lambda node: len(adjacency[node] & p)) if union else None
        candidates = p - (adjacency[pivot] if pivot else set())
        for vertex in list(candidates):
            bron_kerbosch(
                r + [vertex],
                p & adjacency[vertex],
                x & adjacency[vertex],
            )
            p.remove(vertex)
            x.add(vertex)

    bron_kerbosch([], set(nodes), set())
    return cliques


def density(nodes: set[str], adjacency: dict[str, set[str]]) -> float:
    count = len(nodes)
    if count < 2:
        return 0.0
    directed_edge_count = sum(
        1 for gene_a in nodes for gene_b in adjacency[gene_a] if gene_b in nodes
    )
    return directed_edge_count / (count * (count - 1))


def percentile_average(values: dict[str, float], nodes: list[str]) -> dict[str, float]:
    groups: dict[float, list[str]] = {}
    for node in nodes:
        groups.setdefault(values[node], []).append(node)
    position = 1
    count = len(nodes)
    result: dict[str, float] = {}
    for value in sorted(groups, reverse=True):
        tied = len(groups[value])
        average_rank = (position + position + tied - 1) / 2.0
        for node in groups[value]:
            result[node] = (count - average_rank) / max(count - 1, 1)
        position += tied
    return result


def mcode_node_weights(
    genes: list[str], adjacency: dict[str, set[str]]
) -> dict[str, float]:
    weights: dict[str, float] = {}
    for gene in genes:
        if len(adjacency[gene]) < MCODE_DEGREE_CUTOFF:
            weights[gene] = 0.0
            continue
        neighborhood = adjacency[gene] | {gene}
        subgraph = {node: adjacency[node] & neighborhood for node in neighborhood}
        numbers = core_numbers(neighborhood, subgraph)
        highest_core = max(numbers.values()) if numbers else 0
        core = {node for node, value in numbers.items() if value >= highest_core}
        core_adjacency = {node: subgraph[node] & core for node in core}
        weights[gene] = highest_core * density(core, core_adjacency) if core else 0.0
    return weights


def find_mcode_modules(
    genes: list[str], adjacency: dict[str, set[str]]
) -> list[tuple[float, list[str], int]]:
    node_weights = mcode_node_weights(genes, adjacency)
    assigned: set[str] = set()
    modules: list[tuple[float, list[str], int]] = []
    seed_order = sorted(genes, key=lambda node: (node_weights[node], node), reverse=True)

    for seed in seed_order:
        if seed in assigned or node_weights[seed] <= 0:
            continue
        minimum_weight = node_weights[seed] * (1.0 - MCODE_NODE_SCORE_CUTOFF)
        module = {seed}
        stack = [(seed, 0)]
        while stack:
            node, depth = stack.pop()
            if depth >= MCODE_MAX_DEPTH:
                continue
            for neighbor in adjacency[node]:
                if (
                    neighbor not in assigned
                    and neighbor not in module
                    and node_weights[neighbor] >= minimum_weight
                ):
                    module.add(neighbor)
                    stack.append((neighbor, depth + 1))

        if len(module) < MCODE_K_CORE + 1:
            continue
        numbers = core_numbers(module, {node: adjacency[node] & module for node in module})
        # Haircut: retain the maximal 2-core of the detected cluster.
        module = {node for node, value in numbers.items() if value >= MCODE_K_CORE}
        if len(module) < MCODE_K_CORE + 1:
            continue
        edge_count = sum(len(adjacency[node] & module) for node in module) // 2
        score = len(module) * density(module, adjacency)
        modules.append((score, sorted(module), edge_count))
        assigned |= module

    modules.sort(key=lambda item: (-item[0], item[1]))
    return modules


def score_network(
    genes: list[str], edges: list[tuple[str, str, float]]
) -> tuple[list[dict[str, object]], list[tuple[float, list[str], int]], dict[str, set[str]]]:
    adjacency = {gene: set() for gene in genes}
    weights: dict[frozenset[str], float] = {}
    for gene_a, gene_b, score in edges:
        adjacency[gene_a].add(gene_b)
        adjacency[gene_b].add(gene_a)
        weights[frozenset((gene_a, gene_b))] = score

    degree = {gene: len(adjacency[gene]) for gene in genes}
    strength = {
        gene: sum(weights[frozenset((gene, neighbor))] for neighbor in adjacency[gene])
        for gene in genes
    }
    mcc = {gene: 0.0 for gene in genes}
    for clique in maximal_cliques(genes, adjacency):
        if len(clique) >= 2:
            value = math.factorial(len(clique) - 1)
            for gene in clique:
                mcc[gene] += value

    modules = find_mcode_modules(genes, adjacency)
    module_score = {gene: 0.0 for gene in genes}
    module_id = {gene: "" for gene in genes}
    for index, (score, module_genes, _) in enumerate(modules, 1):
        for gene in module_genes:
            if score > module_score[gene]:
                module_score[gene] = score
                module_id[gene] = f"MCODE_{index}"

    strength_percentile = percentile_average(strength, genes)
    degree_percentile = percentile_average(degree, genes)
    mcc_percentile = percentile_average(mcc, genes)
    mcode_percentile = percentile_average(module_score, genes)

    records: list[dict[str, object]] = []
    for gene in genes:
        record: dict[str, object] = {
            "gene": gene,
            "ppi_strength": strength[gene],
            "degree": degree[gene],
            "mcc": mcc[gene],
            "mcode_module": module_id[gene],
            "mcode_module_score": module_score[gene],
            "ppi_percentile": strength_percentile[gene],
            "degree_percentile": degree_percentile[gene],
            "mcc_percentile": mcc_percentile[gene],
            "mcode_percentile": mcode_percentile[gene],
        }
        record["combined_score"] = sum(
            float(record[key])
            for key in (
                "ppi_percentile",
                "degree_percentile",
                "mcc_percentile",
                "mcode_percentile",
            )
        ) / 4.0
        records.append(record)
    records.sort(key=lambda record: (-float(record["combined_score"]), str(record["gene"])))
    return records, modules, adjacency


def write_csv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, object]],
    *,
    lineterminator: str = "\r\n",
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator=lineterminator,
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    arguments = parse_arguments()
    genes = read_genes(arguments.refresh_string)
    mapping_raw, network_raw, query_source = obtain_string_responses(
        genes, arguments.refresh_string
    )
    edges, mapped_genes = parse_edges(genes, mapping_raw, network_raw)
    records, modules, adjacency = score_network(genes, edges)

    with (OUTPUT_DIR / "ppi_edges.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_a", "gene_b", "combined_score"])
        writer.writerows(edges)

    fields = list(records[0])
    write_csv(OUTPUT_DIR / "ppi_node_scores.csv", fields, records)

    with (OUTPUT_DIR / "ppi_mcode_modules.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["module", "module_score", "gene_count", "edge_count", "genes"])
        for index, (score, module_genes, edge_count) in enumerate(modules, 1):
            writer.writerow(
                [f"MCODE_{index}", score, len(module_genes), edge_count, ";".join(module_genes)]
            )

    unmatched = sorted(set(genes) - mapped_genes)
    write_csv(
        OUTPUT_DIR / "ppi_unmatched_genes.csv",
        ["gene", "status"],
        [
            {"gene": gene, "status": "not_found_in_string_mapping"}
            for gene in unmatched
        ],
    )

    network_nodes = sorted(gene for gene in genes if adjacency[gene])
    summary_rows = [
        {"metric": "query_source", "value": query_source},
        {"metric": "input_genes", "value": len(genes)},
        {"metric": "string_mapping_rows", "value": len(mapped_genes)},
        {"metric": "unmatched_input_genes", "value": len(unmatched)},
        {"metric": "network_nodes_with_edge", "value": len(network_nodes)},
        {"metric": "retained_edges_score_gt_0.4", "value": len(edges)},
        {
            "metric": "isolated_input_genes",
            "value": sum(not adjacency[gene] for gene in genes),
        },
        {"metric": "mcode_modules", "value": len(modules)},
    ]
    write_csv(OUTPUT_DIR / "ppi_run_summary.csv", ["metric", "value"], summary_rows)

    top_mcc = {record["gene"] for record in sorted(records, key=lambda x: -float(x["mcc"]))[:5]}
    top_degree = {
        record["gene"] for record in sorted(records, key=lambda x: -int(x["degree"]))[:5]
    }
    first_module = set(modules[0][1]) if modules else set()
    consensus = sorted(top_mcc & top_degree & first_module)
    write_csv(
        OUTPUT_DIR / "ppi_consensus_candidate_genes.csv",
        ["gene", "top5_mcc", "top5_degree", "mcode_1"],
        [
            {
                "gene": gene,
                "top5_mcc": "TRUE" if gene in top_mcc else "FALSE",
                "top5_degree": "TRUE" if gene in top_degree else "FALSE",
                "mcode_1": "TRUE" if gene in first_module else "FALSE",
            }
            for gene in consensus
        ],
        lineterminator="\n",
    )

    print(
        json.dumps(
            {
                "query_source": query_source,
                "input_genes": len(genes),
                "mapped_genes": len(mapped_genes),
                "network_nodes": len(network_nodes),
                "network_edges": len(edges),
                "mcode_modules": len(modules),
                "consensus_candidate_genes": consensus,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
