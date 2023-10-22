from Bio import SeqIO
import os
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        To get RSCU values from genbank files.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def RSCU_filter(input_file):
    if os.path.basename(input_file).split('.')[-1] != 'gb' and os.path.basename(input_file).split('.')[-1] != 'gbk':
        # print(os.path.basename(input_file).split('.')[-1])
        print("Please input Genbank format file")
    else:
        aa = []
        for rec in SeqIO.parse(input_file, 'genbank'):
            for feature in rec.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['gene'][0] not in aa:
                        aa.append(feature.qualifiers['gene'][0])
        print(f"Total gene number: {len(aa)}")
        # print(aa)

    # delete duplicated regions
    save_results_dir = os.path.abspath(input_file).split('.')[0]
    os.mkdir(save_results_dir)
    remove_duplicated = os.path.join(save_results_dir, 'remove_duplicated.fasta')
    # print(remove_duplicated)
    with open(remove_duplicated, 'w') as ff:
        for i in aa:
            for rec in SeqIO.parse(input_file, format='genbank'):
                c = []
                for feature in rec.features:
                    if feature.type == 'CDS':
                        if feature.qualifiers['gene'][0] == i:
                            c.append(feature.extract(rec.seq))

                if len(c) == 1:
                    ff.write('>{}\n{}\n'.format(i, str(c[0])))
                if len(c) == 2:
                    c.remove(c[0]) if len(c[0]) <= len(c[1]) else c.remove(c[1])
                    ff.write('>{}\n{}\n'.format(i, str(c[0])))
    ff.close()
    # delete cds < 300 bp
    filter_300 = os.path.join(save_results_dir, 'after_filter_300.fasta')
    filter_sequences = open(os.path.join(save_results_dir, 'filter_sequences.txt'), 'w')
    save_sequences = open(os.path.join(save_results_dir, 'save_sequences_name.txt'), 'w')
    with open(filter_300, 'w') as ff:
        number = 0
        for rec in SeqIO.parse(remove_duplicated, format='fasta'):
            if len(rec.seq) >= 300 and rec.seq[:3] == 'ATG':
                save_sequences.write(f"{rec.id}\n")
                number += 1
                ff.write('>{}\n{}\n'.format(rec.id, str(rec.seq)))
            else:
                if len(rec.seq) < 300:
                    filter_sequences.write(f"{rec.id}\t seqeunce length samller than 300 bp\n")
                if rec.seq[:3] != 'ATG':
                    filter_sequences.write(f"{rec.id}\t seqeunce not start with ATG\n")
        print(f"After fillter, only {number} sequences saved!")
        save_sequences.write(f"{number}")
        ff.close()
        save_sequences.close()
        filter_sequences.close()

    # merge sequences
    merge_fasta = os.path.join(save_results_dir, 'merge.fasta')
    merge_fasta_w = open(merge_fasta, 'w')
    merge_fasta_w.write(">all_merge_cds\n")
    for rec in SeqIO.parse(filter_300, format='fasta'):
        merge_fasta_w.write(str(rec.seq))
    merge_fasta_w.close()

    # final step: calculate rscu values
    codon_2_aa = {
        'TGG': 'Trp',   'GGT': 'Gly',   'AGG': 'Arg',   'CGA': 'Arg',
        'TGT': 'Cys',   'GGG': 'Gly',   'AGA': 'Arg',   'CGC': 'Arg',
        'TGC': 'Cys',   'GGA': 'Gly',   'AGC': 'Ser',   'CGG': 'Arg',
        'TAC': 'Tyr',   'GGC': 'Gly',   'AGT': 'Ser',   'CGT': 'Arg',
        'TAT': 'Tyr',   'GAT': 'Asp',   'AAG': 'Lys',   'CAA': 'Gin',
        'TCT': 'Ser',   'GAC': 'Asp',   'AAA': 'Lys',   'CAG': 'Gin',
        'TCG': 'Ser',   'GAA': 'Glu',   'AAC': 'Asn',   'CAT': 'His',
        'TCA': 'Ser',   'GAG': 'Glu',   'AAT': 'Asn',   'CAC': 'His',
        'TCC': 'Ser',   'GCA': 'Ala',   'ACC': 'Thr',   'CCT': 'Pro',
        'TTT': 'Phe',   'GCC': 'Ala',   'ACA': 'Thr',   'CCG': 'Pro',
        'TTC': 'Phe',   'GCG': 'Ala',   'ACG': 'Thr',   'CCA': 'Pro',
        'TTG': 'Leu',   'GCT': 'Ala',   'ACT': 'Thr',   'CCC': 'Pro',
        'TTA': 'Leu',   'GTA': 'Val',   'ATC': 'Ile',   'CTA': 'Leu',
        'TAG': 'Ter',   'GTC': 'Val',   'ATA': 'Ile',   'CTT': 'Leu',
        'TGA': 'Ter',   'GTG': 'Val',   'ATT': 'Ile',   'CTC': 'Leu',
        'TAA': 'Ter',   'GTT': 'Val',   'ATG': 'Met',   'CTG': 'Leu',

        }
    seqs = [rec.seq for rec in SeqIO.parse(merge_fasta, format='fasta')][0]
    print(len(seqs))
    aa_list = []
    seqs = str(seqs)
    # stat number of codons
    for i in range(0, len(seqs), 3):
        amino_codon = seqs[i:i + 3]
        aa_list.append(amino_codon)
    number_codon = []
    for key, value in codon_2_aa.items():
        # print(key, aa_list.count(key), value)
        number_codon.append([key, aa_list.count(key), value])
    # print(number_codon)
    # calculate RSCU values
    save_results_file = os.path.join(os.path.dirname(merge_fasta), 'RSCU_results.txt')
    calculate_results = open(save_results_file, 'w')
    all_single_amino_number = {}
    all_same_number = {}
    for i in range(len(number_codon)):
        amino_number = number_codon[i][1]
        amino_name = number_codon[i][2]
        if amino_name not in all_single_amino_number.keys():
            all_single_amino_number[amino_name] = amino_number
            all_same_number[amino_name] = 1
        else:
            all_single_amino_number[amino_name] += amino_number
            all_same_number[amino_name] += 1
    for i in range(len(number_codon)):
        codon_name = number_codon[i][0]
        codon_number = number_codon[i][1]
        amino_name = number_codon[i][2]
        # rfsc_value = round(codon_number / (all_single_amino_number[amino_name]), 2)
        rscu_value = round(codon_number/(all_single_amino_number[amino_name]/all_same_number[amino_name]), 2)
        calculate_results.write(f"{amino_name}\t{codon_name}\t{codon_number}\t{rscu_value}\n")
        # print(amino_name, codon_name, codon_number, rfsc_value)
    calculate_results.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument("-i", "--work_dir", help="Input directory of genbank files")
    args = parser.parse_args()
    for i in os.listdir(args.work_dir):
        if i.endswith('gb') or i.endswith('gbk'):
            file_path = os.path.join(args.work_dir, i)
            print(file_path)
            RSCU_filter(file_path)
        else:
            print(f"please input genbank format files and endwith 'gb' or 'gbk'! ")

print("All Done !!!")
