import argparse
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import polars as pl
from tqdm import tqdm


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="Input tsv file with your modules")
    parser.add_argument("-c", type=str, help="Input tsv file with CONJScan modules")
    parser.add_argument("-w", type=str, help="Number of workers for parallel processing", default=8)
    parser.add_argument("-o", type=str, help="Output file path without extension")
    args = parser.parse_args()

    n_workers = int(args.w)
    output_path = args.o
    conjscan_modules = read_conjscan(args.c)


    # Get all files in the input directory
    files = sorted(
        os.path.join(args.i, entry.name)
        for entry in os.scandir(args.i)
        if entry.is_file()
    )


    with open(output_path, "w", encoding="utf-8") as outfile:

        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=init_pool,
            initargs=(conjscan_modules,)
        ) as executor:

            futures = [
                executor.submit(process_file, file)
                for file in files
            ]

            for future in tqdm(as_completed(futures), total=len(futures)):

                name, values_list = future.result()

                for value in values_list:
                    outfile.write(
                        f"{name},{','.join(map(str, value))}\n"
                    )



def read_my_modules(file_path: str) -> defaultdict[str, list[set[int]]]:

    modules = defaultdict(list)

    df = pl.read_csv(
        file_path,
        separator="\t",
        has_header=True,
    )

    result = (
        df
        .group_by(["community_id", "replicon"])
        .agg(pl.col("gene"))
    )

    for row in result.iter_rows(named=True):

        modules[row["replicon"]].append(
            set(row["gene"])
        )

    return modules


def read_conjscan(file_path) -> defaultdict[str, set[int]]:

    modules = defaultdict(set)

    with open(file_path, "r", encoding="utf-8") as file:
        next(file)  # Pula header
        for line in file:
            _, plasmid, gene_n, _ = line.strip().split("\t")
            modules[plasmid].add(int(gene_n))

    return modules


def init_pool(conjscan_modules):

    global GLOBAL_CONJSCAN
    GLOBAL_CONJSCAN = conjscan_modules



def process_file(filepath) -> tuple[str,list[list[float]]]:

    global GLOBAL_CONJSCAN
    my_modules: defaultdict[str, list[set[int]]] = read_my_modules(filepath)

    metrics_values = []
    for genome, conjscan_module in GLOBAL_CONJSCAN.items():

        if not conjscan_module:
            continue
        modules = my_modules.get(genome)

        if not modules:
            metrics_values.append([0,0,len(conjscan_module)])
        else: 
            metrics_values.append(metrics(conjscan_module, modules))

    return os.path.basename(filepath), metrics_values

def metrics(conjscan_module: set[int], my_modules: list[set[int]], beta=1) -> list[int]:

    metrics_max = [0,0,len(conjscan_module)]
    f1_max = 0
    for module in my_modules:

        TP = len(conjscan_module & module)
        FP = len(module - conjscan_module)
        FN = len(conjscan_module) - TP

        f1 = 2*TP/(2*TP+FP+FN) if 2*TP+FP+FN > 0 else 0

        if f1>f1_max:
            f1_max = f1
            metrics_max = [TP,FP,FN]

        '''
        f1 = 2*TP/(2*TP+FP+FN) if 2*TP+FP+FN > 0 else 0
        recall = TP/(TP+FN) if TP+FN>0 else 0
        precision = TP/(TP+FP) if TP+FP >0 else 0 
        delta_len = len(conjscan_module)-len(module)

        metrics = [f1, recall, precision, delta_len]

        metrics_max = metrics_max if metrics[0] < metrics_max[0] else metrics
        '''

    return metrics_max


if __name__ == "__main__":
    main()