import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """Find four regions in chloroplast genomes.
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""

        return f"{description}\n\n{help_text}"


def find_repeat_regions(input_file):
    replace_dit = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c"
    }
    print(os.path.basename(input_file))
    my_seq = ''
    for rec in SeqIO.parse(input_file, 'fasta'):
        my_seq = str(rec.seq).lower()
    rev_seq = str(Seq(my_seq).reverse_complement())
    # seed test
    all_list = []
    for i in range(0, len(my_seq), 500):
        if i + 500 <= len(my_seq):
            test_seq = my_seq[i:i + 500]
            if test_seq in rev_seq:
                all_list.append([i, i + 500])
    # print(all_list)
    if not all_list:
        print("No repeated sequenced longer than 1,000 bp was detected!")
    else:
        merged_sublists = []
        start = all_list[0][0]
        end = all_list[0][1]

        for sublist in all_list[1:]:
            sublist_start = sublist[0]
            sublist_end = sublist[1]

            if sublist_start == end:
                end = sublist_end
            else:
                merged_sublists.append([start, end])
                start = sublist_start
                end = sublist_end

        merged_sublists.append([start, end])

        # 找出最长的重复序列
        longest_sublist = max(merged_sublists, key=lambda x: x[1] - x[0])

        # 找到最长序列的反向互补序列的位置
        start_pos, start_end = longest_sublist[0], longest_sublist[1]
        max_seq = str(Seq(my_seq[start_pos:start_end]).reverse_complement())
        start_loc = my_seq.find(max_seq)
        seq_end = start_loc + len(max_seq)
        start_list = [start_loc, start_loc + len(max_seq)]

        # 判断前后关系

        def sort_lists(list1, list2):
            if list1[0] < list2[0]:
                return [list1, list2]
            else:
                return [list2, list1]

        _list = sort_lists(start_list, longest_sublist)
        sorted_list = [item for sublist in _list for item in sublist]

        # 延伸并确定出重复区域的序列

        start1, end1, start2, end2 = sorted_list[0], sorted_list[1], sorted_list[2], sorted_list[3]

        # 内侧延伸
        while my_seq[end1 - 1] == replace_dit[my_seq[start2]]:
            end1 += 1
            start2 -= 1

        # print(start1, end1, start2, end2)
        # IR 不在1处
        if start1 != 0:
            if end2 != len(my_seq):
                # 继续延伸
                if my_seq[start1 - 1] == replace_dit[my_seq[end2]]:
                    new_start = start1 - 1
                    new_end = end2
                    while new_end <= len(my_seq) - 1 and my_seq[new_start] == replace_dit[my_seq[new_end]]:
                        new_start -= 1
                        new_end += 1
                    # print("******")
                    # print(new_start, new_end)
                    if new_end == len(my_seq):
                        # 跨区域延伸
                        if my_seq[new_start] == replace_dit[my_seq[0]]:
                            new_end = 1
                            new_start -= 1
                            while my_seq[new_start] == replace_dit[my_seq[new_end]]:
                                new_start -= 1
                                new_end += 1
                            # LSC在前
                            if (start2 - end1) <= (new_start - new_end):
                                result = (
                                    f"LSC:{new_end + 1}-{new_start + 1}\nIRb:{new_start + 2}-{end1 - 1}\n"
                                    f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{len(my_seq)},1-{new_end}"
                                )
                            # SSC在前
                            else:
                                result = (
                                    f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{len(my_seq)},1-{new_end}\n"
                                    f"SSC:{new_end + 1}-{new_start + 1}\nIRa:{start1 + 2}-{end1 - 1}"
                                )
                        # 不跨区域延伸
                        else:
                            # LSC在前
                            if (start2 - end1) <= new_start:
                                result = (
                                    f"LSC:{1}-{new_start + 1}\nIRb:{new_start + 2}-{end1 - 1}\n"
                                    f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{new_end}"
                                )
                            # SSC在前
                            else:
                                result = (
                                    f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{new_end}\n"
                                    f"SSC:1-{new_start + 1}\nIRa:{new_start + 2}-{end1 - 1}"
                                )
                    # 不跨区域
                    else:
                        # LSC在前
                        if (start2 - end1) <= (len(my_seq) - new_end + new_start):
                            result = (
                                f"LSC:{new_end + 1}-{len(my_seq)},1-{new_start + 1}\nIRb:{new_start + 2}-{end1 - 1}\n"
                                f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{new_end}"
                            )
                        # SSC在前
                        else:
                            result = (
                                f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{new_end}\n"
                                f"SSC:{new_end + 1}-{len(my_seq)},1-{new_start + 1}\nIRa:{start1 + 2}-{end1 - 1}"
                            )
                # 停止延伸
                else:
                    # LSC在前
                    if (start2 - end1) <= (len(my_seq) - end2 + start1):
                        result = (
                            f"LSC:{end2 + 1}-{len(my_seq)},1-{start1 + 1}\nIRb:{start1 + 2}-{end1 - 1}\n"
                            f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{end2}"
                        )
                    # SSC在前
                    else:
                        result = (
                            f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{end2}\n"
                            f"SSC:{end2 + 1}-{len(my_seq)},1-{start1 + 1}\nIRa:{start1 + 2}-{end1 - 1}"
                        )
            # 直接跨区域
            else:
                # 延伸
                if my_seq[start1 - 1] == replace_dit[my_seq[0]]:
                    new_end = 0
                    new_start = start1 - 1
                    while my_seq[new_start] == replace_dit[my_seq[new_end]]:
                        new_start -= 1
                        new_end += 1
                    # LSC在前
                    if (start2 - end1) <= (new_start - new_end):

                        result = (
                            f"LSC:{new_end + 1}-{new_start + 1}\nIRb:{new_start + 2}-{end1 - 1}\n"
                            f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{len(my_seq)},1-{new_end}"
                        )
                    # SSC在前
                    else:
                        result = (
                            f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{len(my_seq)},1-{new_end}\n"
                            f"SSC:{new_end + 1}-{new_start}\nIRa:{start1 + 1}-{end1 - 1}"
                        )
                # 不延伸
                else:
                    if (start2 - end1) <= start1:
                        result = (
                            f"LSC:1-{start1}\nIRb:{start1 + 1}-{end1 - 1}\n"
                            f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{end2}"
                        )
                    else:
                        result = (
                            f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{end2}\n"
                            f"SSC:{end2 + 1}-{len(my_seq)},1-{start1}\nIRa:{start1 + 1}-{end1 - 1}"
                        )

        # IR 区域在1处
        else:
            # 跨区域延伸
            if my_seq[-1] == replace_dit[my_seq[end2]]:
                new_start = len(my_seq) - 1
                new_end = end2
                while my_seq[new_start] == replace_dit[my_seq[new_end]]:
                    new_start -= 1
                    new_end += 1
                # 判断LSC 和 SSC
                # LSC在前
                if (start2 - end1) <= (new_start - new_end):
                    result = (
                        f"LSC:{new_end + 1}-{new_start + 1}\nIRb:{new_start + 2}-{len(my_seq)},1-{end1 - 1}\n"
                        f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{end2}"
                    )
                # SSC在前
                else:
                    result = (
                        f"LSC:{end1}-{start2 + 1}\nIRb{start2 + 2}-{new_end}\nSSC{new_end + 1}-{new_start + 1}\n"
                        f"IRa:{new_start + 2}-{len(my_seq)},1-{end1 - 1}"
                    )
            # 不跨区域延伸
            else:
                # LSC 在前
                if (start2 - end1) <= (len(my_seq) - end2 + start1):
                    result = (
                        f"LSC:{end2 + 1}-{len(my_seq)}\nIRb:{start1 + 1}-{end1 - 1}\n"
                        f"SSC:{end1}-{start2 + 1}\nIRa:{start2 + 2}-{end2}"
                    )
                else:
                    result = (
                        f"LSC:{end1}-{start2 + 1}\nIRb:{start2 + 2}-{end2}\n"
                        f"SSC:{end2 + 1}-{len(my_seq)}\nIRa:{start1 + 1}-{end1 - 1}"
                    )
        return result


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--fasta_file', help='fasta format file', required=True)
    args = parser.parse_args()

    input_file = args.fasta_file
    result = find_repeat_regions(input_file)
    if result:
        print(result)


if __name__ == '__main__':
    main()
