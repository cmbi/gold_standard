#!/usr/bin/python
import argparse
import json
import logging
import sys

from gold_standard.html_handler import HtmlHandler
from gold_standard.parsers.aln3SSP import parse_3SSP
from gold_standard.parsers.csv_parser import (
    parse_csv_alignment)
from gold_standard.parsers.gold import (parse_gold_pairwise,
                                                          parse_gold_multi)
from gold_standard.parsers.fasta import parse_fasta
from gold_standard.parsers.fatcat import parse_fatcat
from gold_standard.num_seq import (core_aln_to_num,
                                                     get_core_indexes)
from gold_standard.aln_analyzer import calc_scores_3dm
from gold_standard.result_processor import process_results

# use frozensets of sequence ids as keys in dictionaries
# (regular sets cannot be used because they are mutable)
fs = frozenset

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def calculate_aln_quality(paths, output, in_format, multi, write_json):
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
    if not gold_in['ids']:
        raise RuntimeError("No gold standard alignments were found")
    _log.debug("Sequences in the gold alignment: %s", gold_in['ids'])

    # parse and assess test alignments
    strcts_order = []
    if in_format != 'csv':
        if in_format == 'fatcat' or in_format == '3dm':
            aln_dict = parse_fatcat(paths['aln_path'], gold_in['ids'])
        elif in_format == 'fasta':
            aln_dict = parse_fasta(paths['aln_path'], gold_in['ids'])
        elif in_format == '3SSP':
            aln_dict, strcts_order = parse_3SSP(paths['aln_path'])
        else:
            raise Exception("Invalid input format: {}".format(in_format))
        # create alignment of grounded sequences
        num_aln_dict, core_indexes = core_aln_to_num(
            aln_dict, gold_in['full_seq'], golden_ids=gold_in['ids'])
    else:
        # input format is 'csv'
        aln_dict, num_aln_dict, core_indexes = parse_csv_alignment(
            paths['aln_path'], gold_in['ids'])
    tmpl_no = len(num_aln_dict['cores'])
    if tmpl_no < 2:
        raise Exception("Test alignment has fewer than 2 sequences from the "
                        "gold standard alignment. Did you provide any "
                        "sequences from the gold standard alignment?")
    _log.debug("Sequences in the test alignment: %s",
               str(num_aln_dict['cores'].keys()))
    scores = calc_scores_3dm(gold_in['alns'], num_aln_dict, multi)
    stats = process_results(scores['pairwise'], scores['full'], scores['sp_scores'],
                    output, tmpl_no)

    if write_json:
        # write scores to a json file
        with open(output + ".json", 'w') as o:
            json.dump(stats, o)

    return {
        'wrong_cols': scores["wrong_cols"],
        'aa_aln': aln_dict,
        'num_aln': num_aln_dict,
        'full': gold_in['full_seq'],
        'core_indexes': sorted(core_indexes),
        'order': strcts_order
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
    parser.add_argument("--html", help="HTML output (without variable regions)",
                        action="store_true")
    parser.add_argument("--html_var", help="HTML output with variable regions", action="store_true")
    parser.add_argument("--html_var_short", help="HTML output with shortened "
                        "variable regions", action="store_true")
    parser.add_argument("--input_format", default="fasta")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--json", default=False, action="store_true")
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
    else:
        # no exception traceback when not in debug mode
        sys.tracebacklimit = 0

    # check input format
    allowed_formats = ["fasta", "3dm", "3SSP", "fatcat", "csv"]
    # 3dm - fasta format but variable regions are not in the alignment
    # fatcat - 'final_core'-like format
    # 3SSP - sequence id (no whitespaces) and sequence (corvar) on one line
    #   separated by a comma
    if args.input_format not in allowed_formats:
        parser.error("{} is not an allowed formats. Input format needs to be "
                     "one of the following: {}".format(args.input_format,
                                                       allowed_formats))

    input_paths = {
        'gold_dir': args.gold_dir,
        'gold_path': args.gold_path,
        'aln_path': args.test_aln_path,
        'final_core': args.final_core
    }

    quality_data = calculate_aln_quality(input_paths, args.output,
                                         args.input_format, args.multi, args.json)
    if args.html or args.html_var or args.html_var_short:
        # create html output
        hh = HtmlHandler(var=False, var_short=False)
        hh.write_html(quality_data, args.output)

        if args.html_var or args.html_var_short:
            hh = HtmlHandler(var=args.html_var, var_short=args.html_var_short)
            hh.write_html(quality_data, args.output + "_varshort")
