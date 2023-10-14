import os
from Bio import SeqIO
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Extract common gene sequences from genbank files
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def common_gene_extract(input_file):
    work_dir = os.path.dirname(input_file)
    all_gene = []
    for rec in SeqIO.parse(input_file, format='genbank'):
        for feature in rec.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                if feature.qualifiers['gene'][0] not in all_gene:
                    all_gene.append(feature.qualifiers['gene'][0])
    print(len(all_gene))
    for files in os.listdir(work_dir):
        if files.endswith('gb') or files.endswith('gbk'):
            single_gene = []
            gb_file = os.path.join(work_dir, files)
            print(gb_file)
            for rec in SeqIO.parse(gb_file, format='genbank'):
                for feature in rec.features:
                    if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                        if feature.qualifiers['gene'][0] not in single_gene:
                            single_gene.append(feature.qualifiers['gene'][0])
                print(len(single_gene))
            # delete unigue gene
            for gene_index in range(len(all_gene)-1, -1, -1):
                if all_gene[gene_index].lower() not in [y.lower() for y in single_gene]:
                    all_gene.remove(all_gene[gene_index])
    print("Total gene number is : ", end='\t')
    print(len(all_gene))
    gene_name_file = os.path.join(os.path.dirname(work_dir), 'gene_cp_sort.txt')
    with open(gene_name_file, 'w') as ff:
        for i in all_gene:
            ff.write(f'{i}\n')
    save_dir = os.path.join(os.path.dirname(work_dir), 'common_gene')
    os.mkdir(save_dir)
    for gene_name in all_gene:
        file_name = str(gene_name) + '.fasta'
        file_path = os.path.join(save_dir, file_name)
        with open(file_path, 'w') as fasta_file:
            for gb_file in os.listdir(work_dir):
                if gb_file.endswith('gb') or gb_file.endswith('gbk'):
                    gb_file_path = os.path.join(work_dir, gb_file)
                    fasta_file.write(f">{gb_file.split('.')[0]}\n")
                    for rec in SeqIO.parse(gb_file_path, format='genbank'):
                        my_seqs = []
                        for feature in rec.features:
                            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                                if feature.qualifiers['gene'][0].lower() == gene_name.lower():
                                    my_seqs.append(feature.extract(rec.seq))
                        if len(my_seqs) == 1:
                            fasta_file.write(f"{my_seqs[0]}\n")
                        if len(my_seqs) == 2:
                            my_seqs.remove(my_seqs[0]) if len(my_seqs[0]) <= len(my_seqs[1]) else my_seqs.remove(my_seqs[1])
                            fasta_file.write(f"{my_seqs[0]}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=CustomHelpFormatter
    )
    parser.add_argument('-i', '--input', type=str, required=True, help='Input genbank file path')
    args = parser.parse_args()

    common_gene_extract(args.input)
    print("All Done!")

# next step
# cd common_gene/
# mkdir align_gene
# for i in ./*.fasta; do mafft --thread 10 --auto --thread 10 $i > ./align_gene/$i ;done
