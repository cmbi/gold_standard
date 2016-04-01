import argparse
import logging

from src.gold_standard.aln_analyzer import compare_alignments
from src.gold_standard.html_handler import HtmlHandler
from src.gold_standard.num_seq import core_aln_to_num
from src.gold_standard.parsers.fasta import parse_fasta


FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def run_comparison(aln1_path, aln2_path, outprefix, full_seq_path):
    aln_dict1 = parse_fasta(aln1_path)
    aln_dict2 = parse_fasta(aln2_path)
    full_seq = parse_fasta(full_seq_path)
    num_aln_dict1 = core_aln_to_num(aln_dict1, full_seq, final_core=None)[0]
    num_aln_dict2 = core_aln_to_num(aln_dict2, full_seq, final_core=None)[0]
    comp_result = compare_alignments(num_aln_dict1, num_aln_dict2)
    hh = HtmlHandler()
    hh.write_html(aln_dict1, comp_result["diff_cols1"], outprefix + '1')
    hh.write_html(aln_dict2, comp_result["diff_cols2"], outprefix + '2')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compare two alignments - "
                                                 "write output to html with "
                                                 "highlighted differences")
    parser.add_argument("aln1")
    parser.add_argument("aln2")
    parser.add_argument("full_seq")
    parser.add_argument("outprefix")
    args = parser.parse_args()

    run_comparison(args.aln1, args.aln2, args.outprefix, args.full_seq)
