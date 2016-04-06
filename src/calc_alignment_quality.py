#!/usr/bin/python
import argparse
import logging

from src.gold_standard.html_handler import HtmlHandler
from src.gold_standard.parsers.aln3SSP import parse_3SSP
from src.gold_standard.parsers.gold import parse_gold_pairwise
from src.gold_standard.parsers.var_file import parse_var_file
from src.gold_standard.parsers.fasta import parse_fasta
from src.gold_standard.parsers.fasta import parse_fasta
from src.gold_standard.num_seq import (aln_seq_to_num, core_aln_to_num,
                                       aln_3SSP_to_num)
from src.gold_standard.aln_analyzer import calc_scores, calc_scores_3dm


fs = frozenset

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def calc_sums(aln, id1, id2):
    sum_pos = 0
    sum_neg = 0
    for i in range(len(aln[id1])):
        if ((aln[id1][i] == "-" and aln[id2][i] != "-") or
                (aln[id1][i] != "-" and aln[id2][i] == "-")):
            sum_neg += 1
        elif aln[id1][i] != "-" and aln[id2][i] != "-":
            sum_pos += 2
    return {"pos": sum_pos, "neg": sum_neg}


def calc_stats(confusion_matrices):
    _log.info("Calculating stats")
    stats = {}
    for m_id, m in confusion_matrices.iteritems():
        stats[m_id] = {'specificity': None, 'sensitivity': None,
                       'ppv': None, 'npv': None}
        if m['TN'] + m['FP'] != 0:
            stats[m_id]['specificity'] = float(m['TN']) / (m['TN'] + m['FP'])
        if m['TP'] + m['FN'] != 0:
            stats[m_id]['sensitivity'] = float(m['TP']) / (m['TP'] + m['FN'])
        if m['TP'] + m['FP'] != 0:
            stats[m_id]['ppv'] = float(m['TP']) / (m['TP'] + m['FP'])
        if m['TN'] + m['FN'] != 0:
            stats[m_id]['npv'] = float(m['TN']) / (m['TN'] + m['FN'])
    return stats


def process_results(matrices, full_matrix, sp_scores, output):
    _log.info("Processing the results")
    out_txt = ""

    # FULL MATRIX #
    out_txt += "#### RESULTS ####\n"
    # sensitivity, specificity, ppv, npv
    out_txt += ' '.join(["{}: {}".format(k, v)
                         for k, v in full_matrix.iteritems()]) + '\n'
    # FP, TP, FN, TN values
    full_stats = calc_stats({'full': full_matrix})['full']
    out_txt += ''.join(["{}: {}\n".format(k, v)
                        for k, v in full_stats.iteritems()]) + '\n'
    # average SP score
    out_txt += "SP score: {}\n".format(sum(sp_scores.values()) / len(sp_scores))

    # PAIRWISE stats
    stats = calc_stats(matrices)
    for s_id, s in stats.iteritems():
        header = s_id if (s_id is str) else ' '.join(list(s_id))
        out_txt += "# {}\n".format(header)
        # sensitivity, specificity, ppv, npv
        out_txt += ' '.join(["{}: {}".format(k, v)
                             for k, v in matrices[s_id].iteritems()]) + '\n'
        # FP, TP, FN, TN values
        out_txt += ''.join(["{}: {}\n".format(k, v)
                            for k, v in s.iteritems()]) + '\n'
        # SP score
        out_txt += "SP score: {}\n".format(sp_scores[s_id])

    with open(output, 'w') as out:
        out.write(out_txt)
    _log.info("Created the output file: %s", output)


def calculate_aln_quality(paths, output, in_format, html_out, multi):
    if multi:
        gold_in = parse_var_file(paths['gold_path'])
    else:
        gold_in = parse_gold_pairwise(
            paths['gold_dir'])
    if in_format == '3dm':
        _log.info("Calculating alignment quality in 3DM mode")
        aln_dict = parse_fasta(paths['aln_path'], golden_ids)
        num_aln_dict, core_indexes = core_aln_to_num(
            aln_dict, full_seq, paths['final_core'], golden_ids=golden_ids)
        scores = calc_scores_3dm(gold_alns, num_aln_dict)
    elif in_format == '3SSP':
        aln_dict = parse_3SSP(paths['aln_path'])
        num_aln_dict = aln_3SSP_to_num(aln_dict, full_seq)
        scores = calc_scores_3dm(gold_alns, num_aln_dict)
    elif in_format == 'fasta':
        _log.info("Calculating alignment quality")
        aln_dict = parse_fasta(paths['aln_path'], golden_ids)
        num_aln_dict = {seq_id: aln_seq_to_num(seq)
                        for seq_id, seq in aln_dict.iteritems()}
        scores = calc_scores(gold_alns, num_aln_dict)
    else:
        raise Exception("Invalid input format: {}".format(in_format))
    process_results(scores['pairwise'], scores['full'], scores['SP'], output)
    if html_out:
        return {
            'wrong_cols': scores["wrong_cols"],
            'aa_aln': aln_dict,
            'num_aln': num_aln_dict,
            'full': full_seq,
            'core_indexes': core_indexes
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates quality of the"
                                                 " multiple sequence alignment"
                                                 " based on the provided golden"
                                                 " standard pairwise "
                                                 " alignments")
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
    elif not args.multi and not args.gold_path:
        raise parser.error("In the pairwise (default) mode you must provide the"
                           " gold_dir argument")

    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)
    input_format = "fasta"
    if args.in3dm:
        input_format = "3dm"
    elif args.in3SSP:
        input_format = "3SSP"
    html = (args.html or args.html_var or args.html_var_short)

    input_paths = {
        'gold_dir': args.gold_dir,
        'gold_path': args.gold_path,
        'aln_path': args.test_aln_path,
        'final_core': args.final_core
    }
    quality_data = calculate_aln_quality(input_paths, args.output, input_format,
                                         html, args.multi)
    if html:
        hh = HtmlHandler(var=args.html_var, var_short=args.html_var_short)
        hh.write_html(quality_data, args.output)
