import os
import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    def format_help(self):
        help_text = super().format_help()
        description = """
        Sorted IGS as cp genome order
        Author: Xu wenbo
        Org:    China Pharmaceutical University
        Email:  xwb7533@163.com"""
        return f"{description}\n\n{help_text}"


def sort_as_cp_order(input_file1, input_file2):
    pi_results = open(input_file1, 'r')
    cp_order_results = open(input_file2, 'r')
    results_file_path = os.path.join(os.path.dirname(input_file1), 'IGS_sort_as_cp_order.txt')
    reuslts_file = open(results_file_path, 'w')
    file1_line_list = pi_results.readlines()
    file2_line_list = cp_order_results.readlines()
    for IGS2 in file2_line_list:
        for IGS1 in file1_line_list:
            if IGS1.split('\t')[0] == IGS2.strip():
                print(IGS1, end='')
                reuslts_file.write(IGS1)
    reuslts_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=CustomHelpFormatter)
    parser.add_argument('-i', '--input', type=str, help='Input Pi results file path')
    parser.add_argument('-r', '--reference', type=str, help='Input cp order reference file path')
    args = parser.parse_args()

    if args.input and args.reference:
        sort_as_cp_order(args.input, args.reference)
    else:
        parser.print_help()

print("All Done !!!")

# input_file1 -> Pi_results
# input_file2 -> IGS order in cp genome
# results save in the same directory with input_file1 and named "IGS_sort_as_cp_order.txt"

