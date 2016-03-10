import logging

from custom_exceptions import CustomException
from dict_utils import merge_dicts


fs = frozenset

_log = logging.getLogger(__name__)


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


def score_var_regions(golden_aln, id1, id2, var1):
    matrix = {"FN": 0, "TN": 0}
    for p in var1:
        aln_res = get_aligned_res(p, id1, id2, golden_aln)
        if aln_res == '-':
            matrix["TN"] += 1
        else:
            matrix["FN"] += 1
    return matrix


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
