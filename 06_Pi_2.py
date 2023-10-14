import argparse
import os, sys
from Bio import SeqIO


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Common Intergenic Sequences Extraction
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"



def common_IGS(input_file):
    all_common = []
    for rec in SeqIO.parse(input_file, format='fasta'):
        all_common.append(rec.id)
    print(len(all_common))
    work_dir = os.path.dirname(input_file)
    for fasta_file in os.listdir(work_dir):
        single_IGS = []
        if fasta_file.split('.')[1] == 'fasta':
            fasta_path = os.path.join(work_dir, fasta_file)
            print(fasta_path)
            for rec2 in SeqIO.parse(fasta_path, format='fasta'):
                single_IGS.append(rec2.id)
            for name1_index in range(len(all_common)-1, -1, -1):
                if all_common[name1_index].lower() not in [single_name.lower() for single_name in single_IGS]:
                    all_common.remove(all_common[name1_index])
    print(all_common)
    cp_sort_IGS_file = open(os.path.join(work_dir, 'cp_sort_IGS.txt'), 'w')
    for i in all_common:
        print(i)
        cp_sort_IGS_file.write(f"{i}\n")
    cp_sort_IGS_file.close()
    print(len(all_common))
    save_dir = os.path.join(os.path.dirname(input_file), 'unalign_common_IGS')
    os.mkdir(save_dir)
    for common_name in all_common:
        save_file_path = os.path.join(save_dir, common_name + '.fasta')
        save_file = open(save_file_path, 'w')
        for fasta_file in os.listdir(work_dir):
            fasta_path = os.path.join(work_dir, fasta_file)
            if os.path.isfile(fasta_path):
                if fasta_path.split('.')[1] == 'fasta':
                    for rec in SeqIO.parse(fasta_path, format='fasta'):
                        if rec.id == common_name:
                            save_file.write(f">{fasta_file.split('_IGS')[0]}\n{rec.seq}\n")
        save_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--ref-fasta", help="Input reference of fasta files")
    args = parser.parse_args()
    common_IGS(args.ref_fasta)
# next step
# cd unalign_common_IGS/
# mkdir align_IGS
# for i in ./*.fasta; do mafft --thread 10 --auto  $i > ./align_IGS/$i ;done