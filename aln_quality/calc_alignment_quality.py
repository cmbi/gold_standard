#!/usr/bin/python
import argparse
import logging
import os

from custom_exceptions import CustomException
from html_handler import aln_to_html
from paths import CSS, TEMPLATE
from parsers.aln3SSP import parse_3SSP
from parsers.golden import parse_golden_alns
from parsers.fasta import parse_fasta
from num_seq import aln_seq_to_num, aln_3dm_to_num, aln_3SSP_to_num


fs = frozenset

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def merge_dicts(dict1, dict2):
    """
    Sum up values with the same keys
    """
    return {k: dict1.get(k, 0) + dict2.get(k, 0)
            for k in set(dict1) | set(dict2)}


def get_aligned_res(res_num, query_id, id2, golden_aln):
    aln_pos = golden_aln[query_id].index(res_num)
    return golden_aln[id2][aln_pos]


def calc_pairwise_score(golden_aln, id1, seq1, id2, seq2):
    if len(seq1) != len(seq2):
        raise CustomException("Aligned sequences {} and {} are not of the same "
                              "length")

    matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    sp_score = 0
    wrong_cols = {
        id1: {},
        id2: {}
    }

    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                sp_score += 2
                matrix["TP"] += 2
            else:
                sp_score -= 2
                matrix["FP"] += 2
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq1[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if res2_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq2[i] != '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        # if both are gaps do nothing

    # check output sanity
    res_only = filter(lambda x: x != "-", seq1 + seq2)
    if len(res_only) != sum(matrix.values()):
        raise CustomException("Sum of values in the confusion matrix({}) should"
                              " be equal to the total number of"
                              " residues({})".format(len(res_only),
                                                     sum(matrix.values())))
    sp_max = get_max_sp_score(golden_aln)
    sp_score = float(sp_score) / sp_max

    return {"matrix": matrix, "SP": sp_score, "wrong_cols": wrong_cols}


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


def score_var_regions(golden_aln, id1, id2, var1):
    matrix = {"FN": 0, "TN": 0}
    for p in var1:
        aln_res = get_aligned_res(p, id1, id2, golden_aln)
        if aln_res == '-':
            matrix["TN"] += 1
        else:
            matrix["FN"] += 1
    return matrix


def get_max_sp_score(golden_aln):
    seq1, seq2 = golden_aln.values()
    id1, id2 = golden_aln.keys()
    score = 0
    for i in range(len(golden_aln[id1])):
        if seq1[i] != '-' and seq2[i] != '-':
            score += 2
        elif seq1[i] != '-' or seq2[i] != '-':
            score += 1
    return score


def calc_pairwise_score_3dm(golden_aln, id1, seq1, var1, id2, seq2, var2):
    if len(seq1) != len(seq2):
        raise CustomException("Aligned sequences {} and {} are not of the same "
                              "length: {} and {}".format(seq1, seq2, len(seq1),
                                                         len(seq2)))
    matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    sp_score = 0
    wrong_cols = {
        id1: {},
        id2: {}
    }
    # score core regions
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                sp_score += 2
                matrix["TP"] += 2
            else:
                sp_score -= 2
                matrix["FP"] += 2
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq1[i] != '-' and seq2[i] == '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if res2_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq2[i] != '-' and seq1[i] == '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        # if both are gaps do nothing
    # score variable regions in seq1
    var_matrix = score_var_regions(golden_aln, id1, id2, var1)
    matrix = merge_dicts(matrix, var_matrix)

    # score variable regions in seq2
    var_matrix = score_var_regions(golden_aln, id2, id1, var2)
    matrix = merge_dicts(matrix, var_matrix)

    # check output sanity
    res_only = filter(lambda x: x != "-", seq1 + seq2 + var1 + var2)
    if len(res_only) != sum(matrix.values()):
        raise CustomException("Sum of values in the confusion matrix({}) should"
                              " be equal to the total number of"
                              " residues({})".format(len(res_only),
                                                     sum(matrix.values())))

    sp_max = get_max_sp_score(golden_aln)
    sp_score = float(sp_score) / sp_max

    return {"matrix": matrix, "SP": sp_score, "wrong_cols": wrong_cols}


def calc_scores_3dm(golden_alns, test_aln):
    _log.info("Calculating confusion matrices [3DM mode]")
    matrices_dict = {}
    full_matrix = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    sp_scores = {}
    wrong_cols = {seq_id: {} for seq_id in test_aln["cores"].keys()}
    for id1, seq1 in test_aln["cores"].iteritems():
        for id2, seq2 in test_aln["cores"].iteritems():
            id_set = fs([id1, id2])
            if (id1 != id2 and id_set not in matrices_dict.keys() and
                    id_set in golden_alns.keys()):
                _log.debug("Calculating confusion matrix for sequences {} "
                           "and {}".format(id1, id2))
                scores = calc_pairwise_score_3dm(
                    golden_alns[id_set], id1, seq1, test_aln["var"][id1],
                    id2, seq2, test_aln["var"][id2])
                wrong_cols[id1] = merge_dicts(wrong_cols[id1],
                                              scores["wrong_cols"][id1])
                wrong_cols[id2] = merge_dicts(wrong_cols[id2],
                                              scores["wrong_cols"][id2])
                matrices_dict[id_set] = scores['matrix']
                # add the values to the matrix with overall scores
                full_matrix = merge_dicts(full_matrix, scores['matrix'])
                sp_scores[id_set] = scores['SP']
    return {'pairwise': matrices_dict, 'full': full_matrix, 'SP': sp_scores,
            'wrong_cols': wrong_cols}


def calc_scores(golden_alns, test_aln):
    _log.info("Calculating confusion matrices")
    matrices_dict = {}
    # full matrix holds scores summed up from all sequence pairs
    full_matrix = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    sp_scores = {}
    wrong_cols = set()
    for id1, seq1 in test_aln.iteritems():
        for id2, seq2 in test_aln.iteritems():
            id_set = fs([id1, id2])
            # check if not the same sequence and if not already calculated
            if (id1 != id2 and id_set not in matrices_dict.keys() and
                    id_set in golden_alns.keys()):
                scores = calc_pairwise_score(golden_alns[id_set], id1, seq1,
                                             id2, seq2)
                wrong_cols[id1] = merge_dicts(wrong_cols[id1],
                                              scores["wrong_cols"][id1])
                wrong_cols[id2] = merge_dicts(wrong_cols[id2],
                                              scores["wrong_cols"][id2])
                matrices_dict[id_set] = scores['matrix']
                # add the values to the matrix with overall scores
                full_matrix = merge_dicts(full_matrix, scores['matrix'])
                sp_scores[id_set] = scores['SP']
    return {'pairwise': matrices_dict, 'full': full_matrix, 'SP': sp_scores,
            'wrong_cols': wrong_cols}


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
        num_aln_dict = aln_3dm_to_num(aln_dict, full_seq, golden_ids,
                                      final_core)
        scores = calc_scores_3dm(golden_alns, num_aln_dict)
    elif in3SSP:
        aln_dict = parse_3SSP(test_aln_path)
        num_aln_dict = aln_3SSP_to_num(aln_dict, full_seq, golden_ids,
                                       final_core)
        scores = calc_scores_3dm(golden_alns, num_aln_dict)
    else:
        _log.info("Calculating alignment quality")
        aln_dict = parse_fasta(test_aln_path, golden_ids)
        num_aln_dict = {seq_id: aln_seq_to_num(seq)
                        for seq_id, seq in aln_dict.iteritems()}
        scores = calc_scores(golden_alns, num_aln_dict)
    process_results(scores['pairwise'], scores['full'], scores['SP'], output)
    if html:
        return {
            'wrong_cols': scores["wrong_cols"],
            'aln': aln_dict,
        }


def write_html(quality_data, outname):
    outtxt = aln_to_html(quality_data)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    css_full_path = os.path.join(script_dir, CSS)
    with open(os.path.join(script_dir, TEMPLATE)) as a:
        template_fmt = a.read()
    with open(outname, 'w') as out:
        out.write(template_fmt.format(css_full_path, outtxt))


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
    parser.add_argument("--in3dm", default=False, action="store_true")
    parser.add_argument("--in3SSP", default=False, action="store_true")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--final_core", help="final core file")

    args = parser.parse_args()

    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)

    try:
        quality_data = calculate_aln_quality(args.golden_dir,
                                             args.test_aln_path, args.output,
                                             args.in3dm, args.in3SSP, args.html,
                                             args.final_core)
        if args.html:
            write_html(quality_data, args.output + ".html")
    except CustomException as e:
        _log.error("{}".format(e.message))
        exit(1)
