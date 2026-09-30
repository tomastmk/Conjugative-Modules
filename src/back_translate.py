import argparse
import os


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-p", type=str, help="Input dir with fasta files with protein alignment")
    parser.add_argument("-n", type=str, default=1, help="Input fasta dir with files with nucleotide sequence")
    parser.add_argument("-o", type=str, help="Output file")

    args = parser.parse_args()

    list_of_files = os.listdir(args.p)

    for file_path in list_of_files:

        protein_alignment = f'{args.p}/{file_path}'
        nucleotide_sequence = f'{args.n}/{file_path[:-3]}fna'
        output = f'{args.o}/{file_path[:-3]}fna'

        back_translate(protein_alignment, nucleotide_sequence, output)



def back_translate(protein_alignment_path, fasta_path, output_path):

    with open(protein_alignment_path, 'r') as protein_alignment_file:
     
        protein_alignment = dict()
        for line in protein_alignment_file:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
                protein_alignment[seq_id] = ''
                continue
            protein_alignment[seq_id] = protein_alignment[seq_id] + line.strip("\n")
    


    with open(fasta_path, 'r') as fasta_file:
            
        fasta = dict()
        for line in fasta_file:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
                fasta[seq_id] = ''
                continue
            fasta[seq_id] = fasta[seq_id] + line.strip("\n")

    
    result = dict()


    for protein_id in protein_alignment:

        codon_seq = ''
        protein_seq = protein_alignment[protein_id]
        print(protein_alignment_path)
        nucleotide_seq = fasta[protein_id]
        index = 0

        for amino_acid in protein_seq:
            if amino_acid == '-':
                codon_seq += 3*'-'

            else:
                codon_seq += nucleotide_seq[index:index+3]
                index+=3

        if index != len(nucleotide_seq):
            print(f'Warning: {protein_id}')
            print(f'Used {index} nucleotides, but sequence has {len(nucleotide_seq)}')

        result[protein_id] = codon_seq




    with open(output_path, 'w') as output_file:
        for protein in result:
            output_file.write('>' + protein + '\n')
            output_file.write(result[protein] + '\n')

if __name__ == '__main__':
    main()
