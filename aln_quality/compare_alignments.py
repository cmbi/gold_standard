import argparse
import logging

from aln_analyzer import compare_alignments
from html_handler import write_html_comparison
from num_seq import core_aln_to_num
from parsers.fasta import parse_fasta


FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def run_comparison(aln1_path, aln2_path, outpath, full_seq_path):
    aln_dict1 = parse_fasta(aln1_path)
    aln_dict2 = parse_fasta(aln2_path)
    full_seq = parse_fasta(full_seq_path)
    num_aln_dict1 = core_aln_to_num(aln_dict1, full_seq, final_core=None)
    num_aln_dict2 = core_aln_to_num(aln_dict2, full_seq, final_core=None)
    comp_result = compare_alignments(num_aln_dict1, num_aln_dict2)
    write_html_comparison(comp_result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compare two alignments - "
                                                 "write output to html with "
                                                 "highlighted differences")
    parser.add_argument("aln1")
    parser.add_argument("aln2")
    parser.add_argument("full_seq")
    parser.add_argument("out")
    args = parser.parse_args()

    run_comparison(args.aln1, args.aln2, args.out, args.full_seq)
