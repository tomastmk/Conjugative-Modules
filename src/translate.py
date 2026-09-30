

def main():


def translate(nucleotide_dir, protein_dir):
    nucleotide_file = open(nucleotide_dir, 'r')
    for line in nucleotide_file:

        save = {}
        if line.startswith('>'):
            separated = line.split()[1:]
            save[s

if __name__ == '__main__':
    main()
