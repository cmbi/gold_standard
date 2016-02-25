#!/usr/bin/python
import argparse
import logging
import os
import re

from var_parser import convert_var_to_aln

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def get_var_files(input_dir, ref_seq):
    """
    Get all var files from the input_dir that contain the reference sequence
    """
    if not os.path.exists(input_dir):
        _log.error("Directory not found: {}".format(input_dir))

    filenames_list = os.listdir(input_dir)
    # filter out files not containing the reference sequence and without a .Var
    # extension
    filenames_list = filter(lambda x: ref_seq in x and x.endswith('.Var'),
                            filenames_list)
    file_list = []
    for i in filenames_list:
        with open(os.path.join(input_dir, i)) as a:
            f = a.read().splitlines()
            file_list.append(f)
    return file_list


def remove_ref_gaps(aln, ref_seq):
    ids = aln.keys()
    ids.remove(ref_seq)
    assert len(ids) == 1
    non_ref_id = ids[0]
    aln_seq = ""
    for i in range(len(aln[ref_seq])):
        if aln[ref_seq][i] == "-":
            continue
        aln_seq += aln[non_ref_id][i]
    return {non_ref_id: aln_seq}


def merge_alignments(aln_list, ref_seq):
    ref_ungapped = re.sub('-', '', aln_list[0][ref_seq])
    # first sequence without gaps
    final_aln = {ref_seq: ref_ungapped}
    for aln in aln_list:
        new = remove_ref_gaps(aln, ref_seq)
        final_aln.update(new)
    return final_aln


def merge_goldies(input_dir, ref_seq, output):
    file_list = get_var_files(input_dir, ref_seq)
    aln_list = [convert_var_to_aln(i) for i in file_list]
    print aln_list[0]
    final_aln = merge_alignments(aln_list, ref_seq)
    final_aln_text = [">{}\n{}\n".format(ref_seq, aln_list[0][ref_seq])]
    final_aln_text.extend([">{}\n{}\n".format(seq_id, seq)
                           for seq_id, seq in final_aln.iteritems()])
    with open(output, 'w') as out:
        out.write('\n'.join(final_aln_text))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("input_dir")
    parser.add_argument("ref_seq")
    parser.add_argument("output")
    args = parser.parse_args()
    merge_goldies(args.input_dir, args.ref_seq, args.output)
