import logging

from gold_standard_src.gold_standard.dict_utils import merge_dicts
from gold_standard_src.gold_standard.sanity_checker import (
    check_pairwise_score, check_pairwise_score_3dm)


fs = frozenset

_log = logging.getLogger(__name__)


def compare_pairwise(cores):
    diff_cols1 = {cores['1']['id']: {}, cores['2']['id']: {}}
    diff_cols2 = {cores['1']['id']: {}, cores['2']['id']: {}}
    # res_i1 seq1 aln1
    # res_j1 seq2 aln1
    # res_i2 seq1 aln2
    # res_j2 seq2 aln2
    for i, res_i1 in enumerate(cores['1']['aln1']):
        res_j1 = cores['2']['aln1'][i]
        if res_i1 == '-':
            continue
        if res_i1 in cores['1']['aln2'] and res_j1 in cores['2']['aln2']:
            aln2_index = cores['1']['aln2'].index(res_i1)
            res_j2 = cores['2']['aln2'][aln2_index]
            if res_j2 != res_j1:
                diff_cols1[cores['1']['id']][i] = 1
                diff_cols2[cores['1']['id']][aln2_index] = 1
    return {"diff_cols1": diff_cols1, "diff_cols2": diff_cols2}


def compare_alignments(aln1, aln2):
    diff_cols1 = {seq_id: {} for seq_id in aln1["cores"].keys()}
    diff_cols2 = {seq_id: {} for seq_id in aln2["cores"].keys()}
    for id1, cores1_aln1 in aln1["cores"].iteritems():
        for id2, cores2_aln1 in aln1["cores"].iteritems():
            if (id1 != id2 and id1 in aln2["cores"].keys() and
                    id2 in aln1["cores"].keys()):
                cores1_aln2 = aln2["cores"][id1]
                cores2_aln2 = aln2["cores"][id2]
                cores = {'1':
                         {'id': id1, 'aln1': cores1_aln1, 'aln2': cores1_aln2},
                         '2':
                         {'id': id2, 'aln1': cores2_aln1, 'aln2': cores2_aln2}}
                scores = compare_pairwise(cores)
                diff_cols1[id1] = merge_dicts(diff_cols1[id1],
                                              scores["diff_cols1"][id1])
                diff_cols1[id2] = merge_dicts(diff_cols1[id2],
                                              scores["diff_cols1"][id2])
                diff_cols2[id1] = merge_dicts(diff_cols2[id1],
                                              scores["diff_cols2"][id1])
                diff_cols2[id2] = merge_dicts(diff_cols2[id2],
                                              scores["diff_cols2"][id2])
    return {'diff_cols1': diff_cols1, "diff_cols2": diff_cols2}


def get_max_sp_score(golden_aln):
    seq1, seq2 = golden_aln.values()
    id1 = golden_aln.keys()[0]
    score = 0
    for i in range(len(golden_aln[id1])):
        if seq1[i] == '-' and seq2[i] == '-':
            # ignore this position if both are gaps
            continue
        if seq1[i] == '-' or seq2[i] == '-':
            # insertion / deletion
            score += 1
        else:
            # match / mismatch
            score += 2
    return score


def score_var_regions(golden_aln, id1, id2, var1):
    matrix = {"FN": 0, "TN": 0}
    for p in var1:
        aln_res = get_aligned_res(p, id1, id2, golden_aln)
        if aln_res == '-':
            matrix["TN"] += 1
        else:
            matrix["FN"] += 1
    return matrix


def score_core_regions_3dm(sequences, golden_aln, id1, id2):
    result = {"matrix": {"TP": 0, "FP": 0, "FN": 0, "TN": 0},
              "sp_score": 0,
              "wrong_cols": {i: {} for i in sequences.keys()}}
    # score core regions
    for i, res_i in enumerate(sequences[id1]):
        if res_i != '-' and sequences[id2][i] != '-':
            res2_gold = get_aligned_res(res_i, id1, id2, golden_aln)
            if sequences[id2][i] == res2_gold:
                result['sp_score'] += 2
                result['matrix']['TP'] += 2
            else:
                result['sp_score'] -= 2
                result['matrix']['FP'] += 2
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        elif res_i != '-' and sequences[id2][i] == '-':
            res2_gold = get_aligned_res(res_i, id1, id2, golden_aln)
            if res2_gold == '-':
                result['sp_score'] += 1
                result['matrix']['TN'] += 1
            else:
                result['matrix']['FN'] += 1
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        elif sequences[id2][i] != '-' and sequences[id1][i] == '-':
            res1_gold = get_aligned_res(sequences[id2][i], id2, id1, golden_aln)
            if res1_gold == '-':
                result['sp_score'] += 1
                result['matrix']["TN"] += 1
            else:
                result['matrix']["FN"] += 1
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        # if both are gaps do nothing
    return result


def calc_pairwise_score_3dm(golden_aln, sequences, var_regs):
    id1, id2 = sequences.keys()
    if len(sequences[id1]) != len(sequences[id2]):
        raise Exception("Aligned sequences {} and {} are not of the same "
                        "length: {} and {}".format(
                            sequences[id1], sequences[id2], len(sequences[id1]),
                            len(sequences[id2])))
    result = score_core_regions_3dm(sequences, golden_aln, id1, id2)
    # score variable regions in seq1
    var_matrix = score_var_regions(golden_aln, id1, id2, var_regs[id1])
    result['matrix'] = merge_dicts(result['matrix'], var_matrix)

    # score variable regions in seq2
    var_matrix = score_var_regions(golden_aln, id2, id1, var_regs[id2])
    result['matrix'] = merge_dicts(result['matrix'], var_matrix)

    # check output sanity
    check_pairwise_score_3dm(sequences, var_regs, result, id1, id2)

    sp_max = get_max_sp_score(golden_aln)
    result['sp_score'] = float(result['sp_score']) / sp_max
    return result


def get_aligned_res(res_num, query_id, id2, golden_aln):
    aln_pos = golden_aln[query_id].index(res_num)
    return golden_aln[id2][aln_pos]


def calc_pairwise_score(golden_aln, id1, seq1, id2, seq2):
    if len(seq1) != len(seq2):
        raise Exception("Aligned sequences {} and {} are not of the same "
                        "length".format(seq1, seq2))

    result = {
        "matrix": {"TP": 0, "FP": 0, "FN": 0, "TN": 0},
        "sp_score": 0,
        "wrong_cols": {
            id1: {},
            id2: {}
            }
    }

    for i, res_i in enumerate(seq1):
        if res_i != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(res_i, id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                result['sp_score'] += 2
                result['matrix']["TP"] += 2
            else:
                result['sp_score'] -= 2
                result['matrix']["FP"] += 2
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        elif seq1[i] != '-':
            res2_gold = get_aligned_res(res_i, id1, id2, golden_aln)
            if res2_gold == '-':
                result['sp_score'] += 1
                result['matrix']["TN"] += 1
            else:
                result['matrix']["FN"] += 1
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        elif seq2[i] != '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                result['sp_score'] += 1
                result['matrix']["TN"] += 1
            else:
                result['matrix']["FN"] += 1
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        # if both are gaps do nothing

    # check output sanity
    check_pairwise_score(seq1, seq2, result['matrix'])
    # normalize SP score
    sp_max = get_max_sp_score(golden_aln)
    result['sp_score'] = float(result['sp_score']) / sp_max
    return result


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
                _log.debug("Calculating confusion matrix for sequences %s "
                           "and %s", id1, id2)
                sequences = {
                    id1: seq1,
                    id2: seq2
                }
                var_regs = {
                    id1: test_aln['var'][id1],
                    id2: test_aln['var'][id2]
                }
                scores = calc_pairwise_score_3dm(golden_alns[id_set], sequences,
                                                 var_regs)
                wrong_cols[id1] = merge_dicts(wrong_cols[id1],
                                              scores["wrong_cols"][id1])
                wrong_cols[id2] = merge_dicts(wrong_cols[id2],
                                              scores["wrong_cols"][id2])
                matrices_dict[id_set] = scores['matrix']
                # add the values to the matrix with overall scores
                full_matrix = merge_dicts(full_matrix, scores['matrix'])
                sp_scores[id_set] = scores['sp_score']
    return {'pairwise': matrices_dict, 'full': full_matrix,
            'sp_score': sp_scores, 'wrong_cols': wrong_cols}


def calc_scores(golden_alns, test_aln):
    _log.info("Calculating confusion matrices")
    matrices_dict = {}
    # full matrix holds scores summed up from all sequence pairs
    full_matrix = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    sp_scores = {}
    wrong_cols = dict()
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
                sp_scores[id_set] = scores['sp_score']
    return {'pairwise': matrices_dict, 'full': full_matrix,
            'sp_score': sp_scores, 'wrong_cols': wrong_cols}
