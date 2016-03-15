#!/usr/bin/python
import argparse
import logging

from custom_exceptions import CustomException
from html_handler import write_html
from parsers.aln3SSP import parse_3SSP
from parsers.golden import parse_golden_alns
from parsers.fasta import parse_fasta
from num_seq import aln_seq_to_num, core_aln_to_num, aln_3SSP_to_num
from aln_analyzer import calc_scores, calc_scores_3dm


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
    _log.info("Created the output file: {}".format(output))


def calculate_aln_quality(golden_dir, test_aln_path, output, in3dm, in3SSP,
                          html, final_core=None):
    golden_alns, golden_ids, full_seq = parse_golden_alns(golden_dir)

    if in3dm:
        _log.info("Calculating alignment quality in 3DM mode")
        aln_dict = parse_fasta(test_aln_path, golden_ids)
        num_aln_dict = core_aln_to_num(aln_dict, full_seq, final_core,
                                       golden_ids=golden_ids)
        scores = calc_scores_3dm(golden_alns, num_aln_dict)
    elif in3SSP:
        aln_dict = parse_3SSP(test_aln_path)
        num_aln_dict = aln_3SSP_to_num(aln_dict, full_seq)
        scores = calc_scores_3dm(golden_alns, num_aln_dict)
    else:
        _log.info("Calculating alignment quality")
        aln_dict = parse_fasta(test_aln_path, golden_ids)
        num_aln_dict = {seq_id: aln_seq_to_num(seq)
                        for seq_id, seq in aln_dict.iteritems()}
        scores = calc_scores(golden_alns, num_aln_dict)
    process_results(scores['pairwise'], scores['full'], scores['SP'], output)
    print aln_dict
    if html:
        return {
            'wrong_cols': scores["wrong_cols"],
            'aln': aln_dict,
            'num_aln': num_aln_dict,
            'full_seq': full_seq
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates quality of the"
                                                 " multiple sequence alignment"
                                                 " based on the provided golden"
                                                 " standard pairwise "
                                                 " alignments")
    parser.add_argument("golden_dir")
    parser.add_argument("test_aln_path")
    parser.add_argument("output")
    parser.add_argument("--html", action="store_true")
    parser.add_argument("--htmlvar", action="store_true")
    parser.add_argument("--in3dm", default=False, action="store_true")
    parser.add_argument("--in3SSP", default=False, action="store_true")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--final_core", help="final core file")

    args = parser.parse_args()

    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)

    try:
        quality_data = calculate_aln_quality(
            args.golden_dir, args.test_aln_path, args.output, args.in3dm,
            args.in3SSP, args.html, args.final_core)
        if args.html or args.htmlvar:
            write_html(quality_data["aln"], quality_data["wrong_cols"],
                       args.output, var=args.htmlvar)
    except CustomException as e:
        _log.error("{}".format(e.message))
        exit(1)
