#!/usr/bin/python
import argparse
import logging

from gold_standard_src.gold_standard.html_handler import HtmlHandler
from gold_standard_src.gold_standard.parsers.aln3SSP import parse_3SSP
from gold_standard_src.gold_standard.parsers.gold import (parse_gold_pairwise,
                                                          parse_gold_multi)
from gold_standard_src.gold_standard.parsers.fasta import parse_fasta
from gold_standard_src.gold_standard.num_seq import (
    aln_seq_to_num, core_aln_to_num, aln_3SSP_to_num, get_core_indexes)
from gold_standard_src.gold_standard.aln_analyzer import (calc_scores,
                                                          calc_scores_3dm)
from gold_standard_src.result_processor import process_results


fs = frozenset

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def calculate_aln_quality(paths, output, in_format, multi):
    # read the final_core file if provided
    if paths['final_core']:
        core_indexes = get_core_indexes(paths['final_core'])
    else:
        core_indexes = None
    # read the gold standard alignments
    if multi:
        gold_in = parse_gold_multi(paths['gold_path'])
    else:
        gold_in = parse_gold_pairwise(paths['gold_dir'])
    # parse and assess test alignments
    if in_format == '3dm':
        _log.info("Calculating alignment quality in 3DM mode")
        aln_dict = parse_fasta(paths['aln_path'], gold_in['ids'])
        num_aln_dict = core_aln_to_num(
            aln_dict, gold_in['full_seq'], core_indexes,
            golden_ids=gold_in['ids'])
        scores = calc_scores_3dm(gold_in['alns'], num_aln_dict)
    elif in_format == '3SSP':
        aln_dict = parse_3SSP(paths['aln_path'])
        num_aln_dict = aln_3SSP_to_num(aln_dict, gold_in['full_seq'])
        scores = calc_scores_3dm(gold_in['alns'], num_aln_dict)
    elif in_format == 'fasta':
        _log.info("Calculating alignment quality")
        aln_dict = parse_fasta(paths['aln_path'], gold_in['ids'])
        num_aln_dict = {seq_id: aln_seq_to_num(seq)
                        for seq_id, seq in aln_dict.iteritems()}
        scores = calc_scores(gold_in['alns'], num_aln_dict)
    else:
        raise Exception("Invalid input format: {}".format(in_format))
    process_results(scores['pairwise'], scores['full'],
                    scores['sp_score'], output)
    return {
        'wrong_cols': scores["wrong_cols"],
        'aa_aln': aln_dict,
        'num_aln': num_aln_dict,
        'full': gold_in['full_seq'],
        'core_indexes': core_indexes
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates quality of the"
                                                 " multiple sequence alignment"
                                                 " based on the provided golden"
                                                 " standard pairwise or"
                                                 " multiple alignments")
    parser.add_argument("test_aln_path")
    parser.add_argument("output")
    parser.add_argument("--gold_dir")
    parser.add_argument("--html", action="store_true")
    parser.add_argument("--html_var", action="store_true")
    parser.add_argument("--html_var_short", action="store_true")
    parser.add_argument("--in3dm", default=False, action="store_true")
    parser.add_argument("--in3SSP", default=False, action="store_true")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--final_core", help="final core file")
    parser.add_argument("--multi", action="store_true")
    parser.add_argument("--gold_path")

    args = parser.parse_args()
    # check args
    if args.multi and not args.gold_path:
        raise parser.error("In the 'multi' mode you must provide the gold_path "
                           "argument")
    elif not args.multi and not args.gold_dir:
        raise parser.error("In the pairwise (default) mode you must provide the"
                           " gold_dir argument")
    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)
    # check input format
    input_format = "fasta"
    if args.in3dm:
        input_format = "3dm"
    elif args.in3SSP:
        input_format = "3SSP"

    input_paths = {
        'gold_dir': args.gold_dir,
        'gold_path': args.gold_path,
        'aln_path': args.test_aln_path,
        'final_core': args.final_core
    }

    quality_data = calculate_aln_quality(input_paths, args.output, input_format,
                                         args.multi)
    if args.html or args.html_var or args.html_var_short:
        # create html output
        hh = HtmlHandler(var=args.html_var, var_short=args.html_var_short)
        hh.write_html(quality_data, args.output)
