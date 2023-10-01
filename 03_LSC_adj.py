import os
import argparse
from Bio import SeqIO


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """Adjust the LSC start in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def adjust_start_to_lsc(work_dir, save_dir, info_text):
    info_cont = open(info_text, 'r')
    info_line = info_cont.readline()

    while info_line:
        file_name = info_line.split('\t')[0]
        LSC_loc = info_line.split('\t')[1].split('-')[0].split(':')[1]

        print(file_name)
        save_name = file_name

        work_file = os.path.join(work_dir, file_name)
        save_file = os.path.join(save_dir, save_name)

        # 检查并创建 save_dir 目录
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        with open(save_file, 'w') as ff:
            for rec in SeqIO.parse(work_file, 'fasta'):
                if LSC_loc == '1':
                    ff.write(f'>{file_name}\n{rec.seq}')
                    info_line = info_cont.readline()
                else:
                    adj_seq = str(rec.seq[int(LSC_loc) - 1:]) + str(rec.seq[:int(LSC_loc) - 1])
                    ff.write(f'>{file_name}\n{adj_seq}')
                    info_line = info_cont.readline()

    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", help="Input directory of fasta files")
    parser.add_argument("-o", "--save_dir", help="Output directory for save files")
    parser.add_argument("-f", "--info_text", help="Path to the info text file")
    args = parser.parse_args()

    adjust_start_to_lsc(args.work_dir, args.save_dir, args.info_text)
