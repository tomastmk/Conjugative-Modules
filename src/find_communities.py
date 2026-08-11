import argparse
import sys

import igraph as ig
import leidenalg as leiden
import polars as pl
from tqdm import tqdm


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="Input tsv file with graph")
    parser.add_argument("-o", type=str, help="Output dir")
    parser.add_argument("-mi", type=str, help="Member Id path")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    graph_df = pl.read_csv(args.i, separator=",")

    exp_df = (
        graph_df
        .select(["gene1", "gene2", "exp"])
    )
    jaccard_df = (
        graph_df
        .select(["gene1", "gene2", "jaccard"])
    )
    ochiai_df = (
        graph_df
        .select(["gene1", "gene2", "ochiai"])
    )




    res_list = [i/10 for i in range(1,81)] 
    cutoff_list = [i/100 for i in range(50)]

    exp_graph: ig.Graph = df_to_graph(exp_df)
    jaccard_graph: ig.Graph = df_to_graph(jaccard_df)
    ochiai_graph: ig.Graph = df_to_graph(ochiai_df)

    memberId: pl.DataFrame = pl.read_csv(args.mi, separator="\t")


    member_map = dict(
        zip(
            memberId["cluster_n"],
            memberId["gene"]
        )
    )

    total = len(cutoff_list) * len(res_list) * 3
    with tqdm(total=total) as pbar:

        for cutoff in cutoff_list:

            cut_edges(exp_graph, cutoff)
            cut_edges(jaccard_graph, cutoff)
            cut_edges(ochiai_graph, cutoff)

            for res in res_list:

                find_communities(
                    exp_graph,
                    res,
                    f"{args.o}/exp_com_r{res}_c{cutoff}",
                    member_map
                )
                pbar.update(1)

                find_communities(
                    jaccard_graph,
                    res,
                    f"{args.o}/jaccard_com_r{res}_c{cutoff}",
                    member_map
                )
                pbar.update(1)

                find_communities(
                    ochiai_graph,
                    res,
                    f"{args.o}/ochiai_com_r{res}_c{cutoff}",
                    member_map
                )
                pbar.update(1)

def cut_edges(g: ig.Graph, cutoff: float) -> None:
    edges_remove = g.es.select(weight_lt=cutoff)
    g.delete_edges(edges_remove)



def df_to_graph(df: pl.DataFrame) -> ig.Graph:

    colunas = df.columns

    if len(colunas) != 3:
        print("Erro, formato errado")
        return

    edges = df.select(colunas[0], colunas[1]).rows()
    weights = df.select(colunas[2]).to_series().to_list()

    # Cria o grafo não direcionado
    g = ig.Graph.TupleList(edges, weights=True, directed=False)
    g.es["weight"] = weights

    return g
    

def find_communities(g: ig.Graph, resolution: float, output_path: str, memberId: dict) -> None:


    part = leiden.find_partition(
        g,
        leiden.RBConfigurationVertexPartition,
        weights=g.es["weight"] if "weight" in g.es.attributes() else None,
        resolution_parameter=resolution,
    )

    # Export
    with open(f"{output_path}.tsv", "w", encoding="utf-8") as file:
        file.write("id\tgene_n\tcommunity_n\n")

        for i, community in enumerate(part):
            ids = [g.vs[idx]["name"] for idx in community]
            

            for number in ids:

                protein = memberId[number]

                replicon, gene_n = protein.split("_")
                file.write(f"{replicon}\t{gene_n}\t{i}\n")


if __name__ == "__main__":
    main()
