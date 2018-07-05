import logging

from .dict_utils import merge_dicts
from .sanity_checker import (
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


def get_max_sp_score(golden_aln, id1, id2, multi=False):
    if multi:
        # each residue in the var region gives one point to the overall score
        score = len(golden_aln['var'][id1]) + len(golden_aln['var'][id2])
        seq1 = golden_aln['cores'][id1]
        seq2 = golden_aln['cores'][id2]
    else:
        score = 0
        seq1 = golden_aln[id1]
        seq2 = golden_aln[id2]
    for i, res_i in enumerate(seq1):
        if res_i == '-' and seq2[i] == '-':
            # ignore this position if both are gaps
            continue
        if res_i == '-' or seq2[i] == '-':
            # insertion / deletion
            score += 1
        else:
            # match / mismatch
            score += 2
    return score


def score_var_regions(golden_aln, id1, id2, var1, multi):
    matrix = {"FN": 0, "TN": 0}
    for p in var1:
        aln_res = get_aligned_res(p, id1, id2, golden_aln, multi)
        if aln_res == '-':
            matrix["TN"] += 1
        else:
            matrix["FN"] += 1
    return matrix


def score_core_regions_3dm(sequences, golden_aln, id1, id2, multi):
    result = {"matrix": {"TP": 0, "FP": 0, "FN": 0, "TN": 0},
              "sp_score": 0,
              "wrong_cols": {i: {} for i in sequences.keys()}}
    # score core regions (regions that are part of the core in the TEST
    # alignments - not necessarily core regions in the golden alignments)
    for i, res_i in enumerate(sequences[id1]):
        if res_i != '-' and sequences[id2][i] != '-':
            res2_gold = get_aligned_res(res_i, id1, id2, golden_aln, multi)
            if sequences[id2][i] == res2_gold:
                result['sp_score'] += 2
                result['matrix']['TP'] += 2
            else:
                result['sp_score'] -= 2
                result['matrix']['FP'] += 2
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        elif res_i != '-' and sequences[id2][i] == '-':
            res2_gold = get_aligned_res(res_i, id1, id2, golden_aln, multi)
            if res2_gold == '-':
                result['sp_score'] += 1
                result['matrix']['TN'] += 1
            else:
                result['matrix']['FN'] += 1
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        elif sequences[id2][i] != '-' and sequences[id1][i] == '-':
            res1_gold = get_aligned_res(sequences[id2][i], id2, id1, golden_aln,
                                        multi)
            if res1_gold == '-':
                result['sp_score'] += 1
                result['matrix']["TN"] += 1
            else:
                result['matrix']["FN"] += 1
                result['wrong_cols'][id1][i] = 1
                result['wrong_cols'][id2][i] = 1
        # if both are gaps do nothing
    return result


def calc_pairwise_score_3dm(golden_aln, sequences, var_regs, multi):
    id1, id2 = sequences.keys()
    if len(sequences[id1]) != len(sequences[id2]):
        raise Exception("Aligned sequences {} and {} ({} and {}) are not of "
                        "the same length: {} and {}".format(
                            id1, id2, sequences[id1], sequences[id2],
                            len(sequences[id1]), len(sequences[id2])))
    result = score_core_regions_3dm(sequences, golden_aln, id1, id2, multi)
    # score variable regions in seq1
    var_matrix = score_var_regions(golden_aln, id1, id2, var_regs[id1], multi)
    result['matrix'] = merge_dicts(result['matrix'], var_matrix)

    # score variable regions in seq2
    var_matrix = score_var_regions(golden_aln, id2, id1, var_regs[id2], multi)
    result['matrix'] = merge_dicts(result['matrix'], var_matrix)

    # check output sanity
    check_pairwise_score_3dm(sequences, var_regs, result, id1, id2)

    sp_max = get_max_sp_score(golden_aln, id1, id2, multi)
    result['sp_score'] = float(result['sp_score']) / sp_max
    return result


def get_aligned_res(res_num, query_id, id2, golden_aln, multi=False):
    if not multi:
        aln_pos = golden_aln[query_id].index(res_num)
        aligned_res = golden_aln[id2][aln_pos]
    else:
        if res_num in golden_aln['cores'][query_id]:
            aln_pos = golden_aln['cores'][query_id].index(res_num)
            aligned_res = golden_aln['cores'][id2][aln_pos]
        elif res_num in golden_aln['var'][query_id]:
            aligned_res = '-'
        else:
            raise Exception("Residue {} from sequence {} is not present in "
                            "the gold alignment".format(res_num, query_id))
    return aligned_res


def calc_pairwise_score(golden_aln, id1, seq1, id2, seq2):
    if len(seq1) != len(seq2):
        raise Exception("Aligned sequences {} and {} ({} and {}) are not of the same "
                        "length".format(id1, id2, seq1, seq2))

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
    sp_max = get_max_sp_score(golden_aln, id1, id2)
    result['sp_score'] = float(result['sp_score']) / sp_max
    return result


def calc_scores_3dm(golden_alns, test_aln, multi):
    _log.info("Calculating confusion matrices [3DM mode]")
    result = {
        'pairwise': {},
        'full': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0},
        'sp_scores': {},
        'wrong_cols': {seq_id: {} for seq_id in test_aln["cores"].keys()}
    }
    for id1, seq1 in test_aln["cores"].iteritems():
        for id2, seq2 in test_aln["cores"].iteritems():
            id_set = fs([id1, id2])
            if id1 != id2 and id_set not in result['pairwise'].keys():
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
                if not multi:
                    golden_aln = golden_alns[id_set]
                else:
                    golden_aln = golden_alns

                scores = calc_pairwise_score_3dm(golden_aln, sequences,
                                                 var_regs, multi)
                result['wrong_cols'][id1] = merge_dicts(
                    result['wrong_cols'][id1], scores["wrong_cols"][id1])

                result['wrong_cols'][id2] = merge_dicts(
                    result['wrong_cols'][id2], scores["wrong_cols"][id2])

                result['pairwise'][id_set] = scores['matrix']
                # add the values to the matrix with overall scores
                result['full'] = merge_dicts(result['full'], scores['matrix'])
                result['sp_scores'][id_set] = scores['sp_score']
    return result


def calc_scores_3dm_complex(gold_aln_data, test_aln, mode="strict"):
    """
    """
    target_id = gold_aln_data["target"]
    gold_alns = gold_aln_data["alns"]
    score_mods = gold_aln_data["score_modifiers"]

    result_cores = compare_cores_complex(gold_alns, score_mods, target_id, test_aln)
    n = result_cores["n"]
    overall_score = result_cores["overall_score"]

    if mode == "strict":
        # this is to penalize false negatives, only done in the strict mode
        result_vars = compare_vars_complex(gold_alns, score_mods, target_id, test_aln)
        n += result_vars["n"]
        overall_score += result_vars["overall_score"]

    overall_score /= n
    return overall_score


def compare_cores_complex(gold_alns, score_mods, target_id, test_aln):
    # number of comparisons for score normalization
    n = 0
    overall_score = 0
    # check the aligned ones in test alignment (FPs and TPs)
    test_target_aln = test_aln["cores"][target_id]

    for seq_id, aln in test_aln["cores"].iteritems():
        if seq_id == target_id:
            continue

        gold_aln = gold_alns[seq_id]

        pairwise_score = 0

        for i, res_number in enumerate(aln):
            if res_number == '-':
                continue
            # alignment scores on this position
            gold_residue_scores = gold_aln[str(res_number)]

            # target residue aligned with residue 'res_number' from sequence 'seq_id'
            test_target_res = str(test_target_aln[i])

            found_aln = False
            for res_number_gold, score_category in gold_residue_scores.iteritems():
                if res_number_gold == '*':
                    # means this residue should not be aligned with anything, give a penalty
                    break

                if res_number_gold == test_target_res:
                    n += 1
                    if score_category.startswith('m'):
                        pairwise_score += get_m_score(score_category, score_mods["a"])
                        found_aln = True
                    else:
                        pairwise_score += score_mods[score_category]
                        found_aln = True
                    break

            if not found_aln and test_target_res != '-':
                # means this residue is aligned while it should not be aligned with anything
                pairwise_score -= 1
                n += 1

        overall_score += pairwise_score
    return {"n": n, "overall_score": overall_score}


def get_m_score(m, full_score):
    """
    Get a score for a position with multiple solutions
    :param m: m modifier name, m + int
    :param full_score: full score
    :return: score (float)
    """

    part = int(m.lstrip('m'))

    return (1. / part) * full_score


def compare_vars_complex(gold_alns, score_mods, target_id, test_aln):
    # number of comparisons for score normalization
    n = 0
    overall_score = 0
    # check the aligned ones in test alignment (FPs and TPs)
    # test_target_aln = test_aln["vars"][target_id]

    penalty = score_mods["p"]
    alrights = 0
    for seq_id, gold_aln in gold_alns.iteritems():
        test_seq_nonaligned = test_aln["var"][seq_id]

        for test_res_number in test_seq_nonaligned:
            test_res_number = str(test_res_number)
            if test_res_number not in gold_aln:
                alrights += 1
                # residue not aligned in gold one, OK
                continue

            scores = gold_aln[test_res_number]
            if scores.keys()[0] == '*':
                alrights += 1
                # residue not aligned in gold one, OK
                # means residue should not be aligned with anything, also OK
                continue

            # if we're here that means the residue is not aligned in the
            # test alignment and aligned in gold
            overall_score -= penalty
            n += 1
    return {"n": n, "overall_score": overall_score}
