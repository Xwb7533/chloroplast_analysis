import argparse
from Bio import SeqIO
from collections import defaultdict


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Compare gene counts in two GenBank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com
        """

        return f"{description}\n\n{help_text}"


def compare_genbank_files(ref_file, test_file):
    ref_dict, test_dict = defaultdict(int), defaultdict(int)

    # ref file
    for rec in SeqIO.parse(ref_file, 'genbank'):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA':
                if feature.qualifiers.get('gene'):
                    gene_name = feature.qualifiers.get('gene')[0]
                    ref_dict[gene_name] += 1
                else:
                    gene_loc = feature.location
                    gene_start = int(gene_loc.start) + 1
                    gene_end = int(gene_loc.end)
                    print(f"[{gene_start, gene_end}]  in reference gb file has no gene name!")

    # test file
    for rec in SeqIO.parse(test_file, 'genbank'):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA':
                if feature.qualifiers.get('gene'):
                    gene_name = feature.qualifiers.get('gene')[0]
                    test_dict[gene_name] += 1
                else:
                    gene_loc = feature.location
                    gene_start = int(gene_loc.start) + 1
                    gene_end = int(gene_loc.end)
                    print(f"[{gene_start, gene_end}]  in test gb file has no gene name!")

    print(f"Reference genbank file has {len(ref_dict)} genes!")
    print(f"Test genbank file has {len(test_dict)} genes!")

    # Comparative
    for gene in ref_dict.keys():
        if ref_dict[gene] != test_dict[gene]:
            print(f"{gene} in reference file number is {ref_dict[gene]} \
                but in test file is {test_dict[gene]}.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-r', '--ref_file', help='reference GenBank file', required=True)
    parser.add_argument('-i', '--test_file', help='testing GenBank file', required=True)
    args = parser.parse_args()

    compare_genbank_files(args.ref_file, args.test_file)