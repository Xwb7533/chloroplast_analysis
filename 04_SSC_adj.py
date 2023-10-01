import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def adjust_SSC_forward(work_dir, save_dir, info_text):
    info_cont = open(info_text, 'r')
    info_line = info_cont.readline()

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    while info_line:
        file_name = info_line.split('\t')[0]
        SSC_loc_start = int(info_line.split('\t')[3].split('-')[0].split(':')[1]) - 1
        SSC_loc_end = int(info_line.split('\t')[3].split('-')[1].strip())
        print(file_name)
        save_name = file_name
        work_file = os.path.join(work_dir, file_name)
        save_file = os.path.join(save_dir, save_name)

        with open(save_file, 'w') as ff:
            for rec in SeqIO.parse(work_file, 'fasta'):
                SSC_seq = Seq(str(rec.seq[SSC_loc_start:SSC_loc_end])).reverse_complement()
                rev_seq = str(rec.seq[:SSC_loc_start]) + str(SSC_seq) + str(rec.seq[SSC_loc_end:])
                ff.write(f'>{rec.id}\n{rev_seq}\n')
            info_line = info_cont.readline()

    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to adjust SSC forward")
    parser.add_argument("-i", "--work_dir", help="Input directory of fasta files")
    parser.add_argument("-o", "--save_dir", help="Output directory for save files")
    parser.add_argument("-f", "--info_text", help="Path to the info text file")
    args = parser.parse_args()

    adjust_SSC_forward(args.work_dir, args.save_dir, args.info_text)