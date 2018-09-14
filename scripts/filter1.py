#/usr/bin/python

import os
import sys
from collections import OrderedDict

#filter the files
def filter_files(folder_path):
    file_list = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".final.out") and filename != "Log.final.out":
            file_list.append(filename)
        else:
            continue
    return file_list

#filter the rows
def filter_rows_in_file(file_path):
    line_list = []
    with open(file_path, 'r') as fr:
        for line in fr.readlines():
            if ("Number of input reads" in line or
                "Uniquely mapped reads number" in line or
                "Number of reads mapped to multiple loci" in line or
                "Number of chimeric reads" in line
               ):
                line_list.append(line)
    return line_list

def split_to_key_value (file_name, lines):
    d_col_value = OrderedDict()
    d_col_value["filename                               "] = file_name
    for line in lines:
        col = line.split('|')[0].strip()
        val = line.split('|')[1].strip()
        col = col.replace(" ", "_")
        d_col_value[col] = val
    return d_col_value


def save_dicts_to_file(output_file, dicts):
    file_lines = []
    cols = dicts[0].keys()

    # save width for each column
    d_width = {}
    for c in cols:
        d_width[c] = len(c)


    header = " ".join(cols)
    header += "\n"

    file_lines.append(header)

    for d_item in dicts:
        # construct each line
        dict_line = ""
        for k, v in d_item.items():
            wth  = d_width[k]
            final_str = "{value:{width}}".format(value=v, width=wth)
            dict_line += final_str + " "
        dict_line += "\n"
        file_lines.append(dict_line)

    # write all to fils
    with open(output_file, "w") as file_out:
        # write table header
        file_out.writelines(file_lines)

    return file_lines

def merge_dicts_to_line(dicts):
    # merge all dicts to lines
    file_lines = []
    cols = dicts[0].keys()
    header = " ".join(cols) + '\n'
    file_lines.append(header)
    for d_item in dicts:
        # construct each line
        dict_line = ' '.join(d_item.values()) + "\n"
        file_lines.append(dict_line)

    return file_lines

def save_lines_to_file(output_file, file_lines):
    # write all lines to file
    with open(output_file, "w") as file_out:
        # write table header
        file_out.writelines(file_lines)


def save_dicts_to_file2(output_file, dicts):
    lines = merge_dicts_to_line(dicts)
    save_lines_to_file(output_file, lines)

if __name__ == "__main__":
    all_file_dict = []
    file_list = filter_files(".")
    for file in file_list:
        lines = filter_rows_in_file(file)
        file_dict = split_to_key_value(file, lines)
        all_file_dict.append(file_dict)

    save_dicts_to_file("filter_output.txt", all_file_dict)


    with open("filter_output.txt", "r") as f:
        for line in f.readlines():
            sys.stdout.write(line)
