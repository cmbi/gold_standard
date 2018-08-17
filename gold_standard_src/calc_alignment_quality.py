#!/usr/bin/python
import argparse
import json
import logging
import sys

from copy import deepcopy

from gold_standard.aln_processor import make_master_seq_full
from gold_standard.html_handler import HtmlHandler
from gold_standard.parsers.aln3SSP import parse_3SSP
from gold_standard.parsers.csv_parser import (
    parse_csv_alignment)
from gold_standard.parsers.gold import (parse_gold_pairwise,
                                        parse_gold_multi, parse_gold_json)
from gold_standard.parsers.fasta import parse_fasta
from gold_standard.parsers.fatcat import parse_fatcat
from gold_standard.num_seq import (core_aln_to_num,
                                   get_core_indexes)

from gold_standard.aln_analyzer import calc_scores_3dm, calc_scores_3dm_complex
from gold_standard.result_processor import process_results

# use frozensets of sequence ids as keys in dictionaries
# (regular sets cannot be used because they are mutable)
fs = frozenset

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger("__main__")


def detect_input_format(aln_path):
    """
    Check if input format is 3dm or 3SSP
    """
    with open(aln_path) as a:
        first_line = a.read().splitlines()[0]

    if first_line.isupper() and len(first_line.split()) > 2:
        return "3SSP"
    else:
        return "3dm"


def parse_input_alignment(aln_path, full_seq, gold_ids, in_format, final_core_path, master_id):
    """
    parse and assess test alignments
    """
    # read the final_core file if provided
    if final_core_path:
        core_indexes = get_core_indexes(final_core_path)
    else:
        core_indexes = None

    strcts_order = []
    write_pairwise_html = True
    if in_format != 'csv':
        if not in_format:
            # detect input format (only for fatcat, 3dm, 3SSP)
            in_format = detect_input_format(aln_path)

        if in_format == 'fatcat' or in_format == '3dm':
            aln_dict, strcts_order = parse_fatcat(aln_path, gold_ids)
        elif in_format == 'fasta':
            aln_dict = parse_fasta(aln_path, gold_ids)
        elif in_format == '3SSP':
            aln_dict, strcts_order = parse_3SSP(aln_path)
        else:
            raise Exception("Invalid input format: {}".format(in_format))

        # fill in the alignment with gaps so that the full master sequence is in
        # the test alignment
        old_aln_dict = deepcopy(aln_dict)
        if master_id in aln_dict:
            aln_dict = make_master_seq_full(aln_dict, full_seq, gold_ids, master_id)
        else:
            _log.warning("Will not be able to create a pairwise comparison because "
                         "target sequence is not present in the test alignment")
            write_pairwise_html = False

        # create alignment of grounded sequences
        num_aln_dict, core_indexes = core_aln_to_num(aln_dict, full_seq, golden_ids=gold_ids)
    else:
        # input format is 'csv'
        aln_dict, num_aln_dict, core_indexes = parse_csv_alignment(aln_path, gold_ids)

    # check if enough sequences
    tmpl_no = len(num_aln_dict['cores'])
    if tmpl_no < 2:
        raise Exception("Test alignment has fewer than 2 sequences from the "
                        "gold standard alignment. Did you provide any "
                        "sequences from the gold standard alignment?")
    _log.debug("Sequences in the test alignment: %s",
               str(num_aln_dict['cores'].keys()))
    return aln_dict, strcts_order, num_aln_dict, core_indexes, write_pairwise_html


def process_per_residue_data(per_residue_scores, target_id, target_seq, max_scores):
    """
    Convert it to a simple dict, {seq_id: res_index: score}}
    """
    wrong_cols = {seq_id: {} for seq_id in per_residue_scores}
    for seq_id, residue_scores in per_residue_scores.iteritems():
        for res_index, score in residue_scores.iteritems():
            try:
                max_score_on_pos = max_scores[seq_id][str(res_index)]
            except:
                max_score_on_pos = 0

            if max_score_on_pos != 0 and score[0]:
                normalized_score = score[1] / max_score_on_pos
            elif score[0]:
                normalized_score = 1
            else:
                normalized_score = -1
            if max_score_on_pos < score[1]:
                msg = "Score on position {} is higher than the max score. {} vs {}".format(
                    res_index, score, max_score_on_pos)
                _log.error(msg)
                raise RuntimeError(msg)

            wrong_cols[seq_id][res_index] = (score[0], normalized_score)

    # add target sequence with a full score (it will not be in
    # the per_residue_scores if it wasn't in the final_core.json)
    if target_id not in wrong_cols:
        wrong_cols[target_id] = {i + 1: (True, 1.0) for i in range(len(target_seq))}

    return wrong_cols


def calculate_aln_quality_complex(paths, output, in_format, write_json):
    # read the gold standard alignments

    gold_path = paths['gold_path']
    corvar_path = gold_path.replace('.json', '.txt.Var')
    gold_in = parse_gold_json(gold_path, corvar_path)

    if not gold_in['ids']:
        raise RuntimeError("No gold standard alignments were found")
    _log.debug("Sequences in the gold alignment: %s", gold_in['ids'])

    # that's the test alignment in a final core format, not necessary
    final_core = paths.get("final_core")

    aln_dict, strcts_order, num_aln_dict, core_indexes, write_pairwise_html = \
        parse_input_alignment(
            paths['aln_path'], gold_in['full_seq'], gold_in['ids'], in_format, final_core,
            gold_in["target"])

    # calculate scores
    scores = calc_scores_3dm_complex(gold_in, num_aln_dict)
    if write_json:
        # write scores to a json file
        with open(output + ".json", 'w') as o:
            json.dump(scores, o, indent=4)

    target_id = gold_in['target']

    wrong_cols = process_per_residue_data(scores['per_residue_scores'], target_id, gold_in['full_seq'][target_id], scores["max_scores"])

    return {
        'target_id': target_id,
        'write_pairwise_html': write_pairwise_html,
        'overall_score': scores['overall_score'],
        'per_residue_scores': scores['per_residue_scores'],
        'wrong_cols': wrong_cols,
        'aa_aln': aln_dict,
        'gold_aln': gold_in['alns'],
        'num_aln': num_aln_dict,
        'full': gold_in['full_seq'],
        'core_indexes': sorted(core_indexes),
        'order': strcts_order
    }

    # stats = process_results(scores['pairwise'], scores['full'], scores['sp_scores'],
    #                         output, len(strcts_order))

    # if write_json:
    #     # write scores to a json file
    #     with open(output + ".json", 'w') as o:
    #         json.dump(stats, o)

    # return {
    #     'wrong_cols': scores["wrong_cols"],
    #     'aa_aln': aln_dict,
    #     'gold_aln': gold_in['alns'],
    #     'num_aln': num_aln_dict,
    #     'full': gold_in['full_seq'],
    #     'core_indexes': sorted(core_indexes),
    #     'order': strcts_order
    # }


def calculate_aln_quality_simple(paths, output, in_format, multi, write_json, gold_json):
    # read the gold standard alignments
    if multi:
        gold_in = parse_gold_multi(paths['gold_path'])
    else:
        gold_in = parse_gold_pairwise(paths['gold_dir'])

    if not gold_in['ids']:
        raise RuntimeError("No gold standard alignments were found")
    _log.info("'SIMPLE' score calculation")
    _log.debug("Sequences in the gold alignment: %s", gold_in['ids'])

    aln_dict, strcts_order, num_aln_dict, core_indexes, write_pairwise_html = parse_input_alignment(
        paths['aln_path'], gold_in['full_seq'], gold_in['ids'], in_format, paths['final_core'], gold_in["ids"])

    # calculate scores
    scores = calc_scores_3dm(gold_in['alns'], num_aln_dict, multi)
    stats = process_results(scores['pairwise'], scores['full'], scores['sp_scores'],
                            output, len(strcts_order))

    if write_json:
        # write scores to a json file
        with open(output + ".json", 'w') as o:
            json.dump(stats, o)

    return {
        "write_pairwise_html": write_pairwise_html,
        'wrong_cols': scores["wrong_cols"],
        'aa_aln': aln_dict,
        'gold_aln': gold_in['alns'],
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
    parser.add_argument("--html_pair", help="HTML output with pairwise comparisons", action="store_true")
    parser.add_argument("--html_var_short", help="HTML output with shortened "
                        "variable regions", action="store_true")
    parser.add_argument("--input_format")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--json", default=False, action="store_true")
    parser.add_argument("--final_core", help="final core file")
    parser.add_argument("--multi", action="store_true")
    parser.add_argument("--gold_path")
    parser.add_argument("--gold_json", default=False, action='store_true')

    args = parser.parse_args()
    # check args
    if args.multi and not args.gold_path:
        raise parser.error("In the 'multi' mode you must provide the gold_path "
                           "argument")
    elif not args.multi and not args.gold_dir and not args.gold_json:
        raise parser.error("In the pairwise (default) mode you must provide the"
                           " gold_dir argument")
    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)
    else:
        # no exception traceback when not in debug mode
        sys.tracebacklimit = 0

    # check input format
    allowed_formats = ["fasta", "3dm", "3SSP", "fatcat", "csv", "json"]
    # 3dm - fasta format but variable regions are not in the alignment
    # fatcat - 'final_core'-like format
    # 3SSP - sequence id (no whitespaces) and sequence (corvar) on one line
    #   separated by a comma
    if args.input_format and args.input_format not in allowed_formats:
        parser.error("{} is not an allowed formats. Input format needs to be "
                     "one of the following: {}".format(args.input_format,
                                                       allowed_formats))

    input_paths = {
        'gold_dir': args.gold_dir,
        'gold_path': args.gold_path,
        'aln_path': args.test_aln_path,
        'final_core': args.final_core
    }
    if args.gold_json:
        quality_data = calculate_aln_quality_complex(input_paths, args.output,
                                                     args.input_format, args.json)
    else:
        quality_data = calculate_aln_quality_simple(input_paths, args.output,
                                                    args.input_format, args.multi, args.json,
                                                    args.gold_json)
    hh = HtmlHandler()
    if not args.gold_json:
        if args.html_pair and quality_data["write_pairwise_html"]:
            # write pairwise html output
            hh.write_html(quality_data, args.output + "_pairwise", mode="pairwise")

        if args.html or args.html_var or args.html_var_short:
            # create html output
            hh.write_html(quality_data, args.output, mode="cores")

        if args.html_var:
            # create html output with variable regions (full or trimmed)
            hh.write_html(quality_data, args.output + "_var", mode="var")

        if args.html_var_short:
            # create html output with variable regions (full or trimmed)
            hh.write_html(quality_data, args.output + "_varshort", mode="var_short")
    else:
        if args.html_pair and quality_data["write_pairwise_html"]:
            # write pairwise html output
            hh.write_html(quality_data, args.output + "_pairwise", mode="pairwise_complex")

        if args.html or args.html_var or args.html_var_short:
            # create html output
            hh.write_html(quality_data, args.output, mode="cores_complex")
