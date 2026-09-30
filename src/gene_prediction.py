"""
Script simples para predição de genes de procariotas de dados metagenômicos
"""

# Dependecias
import argparse
import sys
from multiprocessing.pool import ThreadPool

from Bio import SeqIO, SeqRecord
from pyrodigal import GeneFinder

global GENE_FINDER
GENE_FINDER: GeneFinder = GeneFinder(meta=True)

# Usage: python gene_prediction.py -i <input_fasta> -o <output_path_without_extension>
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="Input fasta file with genomes")
    parser.add_argument("-o", type=str, help="Output file path. Do not add extension.")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    input_path: str = args.i
    output_dir: str = args.o

    data_iterator: SeqIO.FastaIO.FastaIterator = SeqIO.parse(input_path, "fasta")

    with (
        ThreadPool() as pool,
        open(output_dir + "/pred_genes.faa", "w", encoding="utf-8") as transcription_file,
        open(output_dir + "/pred_proteins.faa", "w", encoding="utf-8") as translation_file,
        open(output_dir + "/gene_count.tsv", "w", encoding="utf8") as cont_file,
    ):
        sequences = pool.imap(predict_genes, data_iterator)
        for i, sequence in enumerate(sequences):
            
            sequence[1].write_gff(
                transcription_file, sequence[0], include_translation_table=True 
            )
            #sequence[1].write_genes(
            #    transcription_file, sequence[0] 
            #)
            #sequence[1].write_translations(
            #    translation_file, sequence[0], include_stop=False
            #)
            cont_file.write(f"{sequence[0]} {len(sequence[1])}\n")
            print(f"Sequências lidas: {i:02d}\r")

def predict_genes(genome_sequence: SeqRecord.SeqRecord) -> tuple[str, str]:
    """Função que recebe uma sequência e retorna uma tupla, contendo o id e uma lista com todos genes.

    Parameters
    ----------
    genome_sequence : SeqRecord.SeqRecord
        genome_sequence

    Returns
    -------
    tuple[str,str]

    """
    return genome_sequence.id, GENE_FINDER.find_genes(
        genome_sequence.seq._data.decode()
    )




if __name__ == "__main__":
    main()
