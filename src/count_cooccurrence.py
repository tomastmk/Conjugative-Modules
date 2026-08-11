"""A module for creating coocurence graph."""

import argparse
import sys
import time
from array import array
from collections import defaultdict
from functools import partial
from math import exp, sqrt
from multiprocessing import Pool

import polars as pl


# Usage: python count_cooccurrence.py -i <input_cluster> -o <output_dir> -mg <max_distance>
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="Input tsv file with clusters")
    parser.add_argument("-o", type=str, help="Output dir")
    parser.add_argument("-mg", type=str, default=1, help="Max distance between neighbor genes")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Pré-processamento dos dados
    print(f"\n{bcolors.BOLD}Preparando...{bcolors.ENDC}")
    n_clusters: int = tsv_cluster_reader(args.i, args.o)

    # Contagem de Frequência
    print(f"\n{bcolors.BOLD}Calculando frequência de cada cluster:{bcolors.ENDC}")
    frequency = calculate_frequency(args.o + "_cluster_sequence.tsv", n_clusters)

    # Cálculo das métricas
    print(f"\n{bcolors.BOLD}Calculando métricas de coocorrência{bcolors.ENDC}")

    print("serializado")
    matrix = serializado(args.o + "_cluster_sequence.tsv", frequency, int(args.mg))

    df = pairmetrics_to_polars(matrix, frequency)
    df.write_csv(args.o + "_graph.tsv")



class bcolors:
    """Classe para armazenar alguns valores de cores para print"""

    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


# Função wrapper para contabilizar tempo que uma função leva
def time_func(func):
    """time_func.

    Parameters
    ----------
    func :
        func
    """

    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        delta_time = end_time - start_time
        print(f"{func.__name__} demorou {delta_time:.2f} segundos.")
        return result

    return wrapper


class PairMetrics:
    """
    Classe criada representar uma sequência de valores para um determinado par
    """

    def __init__(self):
        # Contagem de quantas genomas possuem o par de genes
        self.cooccurrence_count: int = 0
        # Contagem de quantos genomas possuem o par de genes a no máximo uma distância n
        self.gap_count: int = 0
        # Menor distância entre os genes no genoma
        self.min_dist: int = None
        # Soma dos resultados da métrica de distância exponencial
        self.sum_exp: float = 0

    def add_dist(self, dist):
        """Salva a distância fornecida somente se ele for menor que a anterior.

        Parameters
        ----------
        dist :
            distânica de genes
        """
        if self.min_dist is None or self.min_dist > dist:
            self.min_dist = dist

    def merge(self, other):
        """Une dois objetos do tipo PairMetrics, somando seus atributos.

        Parameters
        ----------
        other :
            outro objeto do tipo PairMetrics
        """
        self.cooccurrence_count += other.cooccurrence_count
        self.gap_count += other.gap_count
        self.sum_exp += other.sum_exp


def calculate_genome_cooccurrence(
    genome: list[int],
    max_dist: int,
    freq: list[int],
    circular: bool = True,
    MIN_FREQ: int = 0,
):
    """Calcula a coocorrência de ints em uma lista de ints.

    Parameters
    ----------
    genome : list[int]
        genome
    max_dist : int
        Distância máxima para considerar vizinhos
    freq : list[int]
        Lista com frequência total de cada int
    circular : bool
        Se a lista é circular ou não
    MIN_FREQ : int
        Valor mínimo de frequência para considerar coocorrência
    """

    local_matrix: dict[PairMetrics] = defaultdict(PairMetrics)
    length = len(genome)

    for pos1, gene1 in enumerate(genome):
        if freq[gene1] < MIN_FREQ:
            continue  # Evita calcular pares com pouca frequência

        lin_dist = 0  # Distância sem considerar circularidade do genoma

        for pos2, gene2 in enumerate(genome[pos1 + 1 :]):

            if freq[gene2] < MIN_FREQ:
                continue

            lin_dist += 1

            if gene1 == gene2:
                continue

            if circular:
                dist = min(lin_dist, length - lin_dist)

            pair_index = frozenset((gene1, gene2))
            local_matrix[pair_index].add_dist(dist)

    for pair in local_matrix:
        dist = local_matrix[pair].min_dist

        if local_matrix[pair].min_dist is None:
            continue  # Não coocorreu

        local_matrix[pair].cooccurrence_count += 1
        local_matrix[pair].sum_exp += exp(1 - dist)

        if dist <= max_dist:
            local_matrix[pair].gap_count += 1

        local_matrix[pair].min_dist = None

    return local_matrix


def pairmetrics_to_polars(
    matrix: dict[frozenset, PairMetrics], frequency: list[int]
) -> pl.DataFrame:
    """Tranforma uma matriz contendo os pares de genes e suas métricas em um Polars DataFrame.

    Parameters
    ----------
    matrix : dict[frozenset, PairMetrics]
        matrix
    frequency : list[int]
        frequency

    Returns
    -------
    pl.DataFrame

    """
    rows = []

    for pair, metrics in matrix.items():

        gene1, gene2 = tuple(pair)
        f1, f2 = frequency[gene1], frequency[gene2]

        rows.append(
            {
                "gene1": gene1,
                "gene2": gene2,
                "freq1": f1,
                "freq2": f2,
                "cooccurrence_count": metrics.cooccurrence_count,
                "gap_count": metrics.gap_count,
                "jaccard": metrics.gap_count / (f1 + f2 - metrics.gap_count),
                "ochiai": metrics.gap_count / sqrt(f1 * f2),
                "exp": metrics.sum_exp / metrics.cooccurrence_count,
            }
        )

    return pl.DataFrame(rows)




@time_func
def tsv_cluster_reader(clustering_path: str, output_dir: str = None):
    """
    Recebe como entrada um arquivo tsv de clusterização do Pyrodigal. Escreve um arquivo tsv com
    replicon na primeira coluna e na segunda a sequência de ids de cluster dos genes, em ordem e
    outro com uma tabela contendo o nome do gene na primeira coluna e na segunda seu id de cluster


    Parameters
    ----------
    clustering_path : str
        clustering_path
    output_path : str
        output_path

    Returns
    -------
    pl.DataFrame

    """

    centroidId_df: pl.LazyFrame = pl.scan_csv(
        clustering_path, separator="\t", has_header=True
    ).with_columns(
        (pl.col("centroid") != pl.col("centroid").shift())
        .fill_null(True)
        .cum_sum()
        .sub(1)
        .alias("cluster_n")
    )

    max_val: int = centroidId_df.select(pl.col("cluster_n").max()).collect().item()

    cluster_df: pl.LazyFrame = (
        centroidId_df.with_columns(
            gene_n=(
                pl.col("member")
                .cast(pl.Utf8)
                .str.replace(
                    r"^.*_", ""
                )  # Indo pela esquerda, tudo antes do primeiro underline
                .cast(pl.Int64)
            ),
            plasmid_name=(
                pl.col("member")
                .cast(pl.Utf8)
                .str.replace(
                    r"_[^_]*$", ""
                )  # Indo pela direita, tudo depois do primeiro underline
                .cast(pl.Categorical)
            ),
        )
        .sort("gene_n")
        .group_by("plasmid_name")
        .agg(pl.col("cluster_n").alias("genome"))
        .with_columns(
            pl.col("genome").list.eval(pl.element().cast(pl.Utf8)).list.join(",")
        )
    )
    centroidId_df = centroidId_df.drop("centroid")
    if output_dir is not None:

        cluster_df.sink_csv(output_dir + "_cluster_sequence.tsv", separator="\t")
        centroidId_df.sink_csv(output_dir + "_member_id.tsv", separator="\t")
        return max_val

    return cluster_df.collect()



@time_func
def calculate_frequency(cluster_path: str, max_cluster: int):
    """Calcula frequência de cada gene.

    Parameters
    ----------
    cluster_path : str
        cluster_path
    max_cluster : int
        max_cluster
    """
    counts = array("I", [0]) * (max_cluster + 1)
    with open(cluster_path, encoding="utf-8") as file:
        next(file)
        for line in file:
            values = line.rstrip().split("\t", 1)[1]
            for v in values.split(","):
                counts[int(v)] += 1

    return counts


def process_line(line, frequency, MAX_GAP):
    genome = list(map(int, line.split("\t")[1].strip("\n").split(",")))
    temp_matrix = calculate_genome_cooccurrence(genome, MAX_GAP, frequency, MIN_FREQ=30)
    return temp_matrix


@time_func
def paralelized(OUTPUT_PATH, frequency, MAX_GAP):
    matrix = defaultdict(PairMetrics)
    worker = partial(process_line, MAX_GAP=MAX_GAP, frequency=frequency)
    with open(OUTPUT_PATH, "r", encoding="utf-8") as file:
        next(file)
        with Pool() as pool:
            for temp_matrix in pool.imap_unordered(worker, file, chunksize=200):
                for pair in temp_matrix:
                    matrix[pair].merge(temp_matrix[pair])


@time_func
def serializado(OUTPUT_PATH, frequency, MAX_GAP):
    matrix = defaultdict(PairMetrics)
    with open(OUTPUT_PATH, "r", encoding="utf-8") as file:
        next(file)
        for line in file:
            genome = list(map(int, line.split("\t")[1].strip("\n").split(",")))
            temp_matrix = calculate_genome_cooccurrence(
                genome, MAX_GAP, frequency, MIN_FREQ=30
            )
            for pair in temp_matrix:
                matrix[pair].merge(temp_matrix[pair])
    return matrix


if __name__ == "__main__":
    main()
