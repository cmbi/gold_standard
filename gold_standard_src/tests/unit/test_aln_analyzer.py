from nose.tools import eq_

import gold_standard_src.gold_standard.aln_analyzer as aa

from gold_standard_src.gold_standard.parsers.aln3SSP import parse_3SSP
from gold_standard_src.gold_standard.num_seq import core_aln_to_num


def test_score_var_regions():
    golden = {"1": [1, 2, 3, '-', '-', 4, 5, 6, '-', '-'],
              "2": ['-', '-', '-', 1, 2, 3, 4, 5, '-', '-']}
    var = [1, 2, 3, 4]
    matrix = aa.score_var_regions(golden, '1', '2', var, multi=False)
    eq_({"TN": 3, "FN": 1}, matrix)


def test_compare_pairwise():
    id1 = '1'
    id2 = '2'
    cores1_aln1 = "-1234--"
    cores2_aln1 = "-1234--"
    cores1_aln2 = "-123-4-"
    cores2_aln2 = "1-23-4-"
    cores = {'1':
             {'id': id1, 'aln1': cores1_aln1, 'aln2': cores1_aln2},
             '2':
             {'id': id2, 'aln1': cores2_aln1, 'aln2': cores2_aln2}}
    res = aa.compare_pairwise(cores)

    expected_result = {
        "diff_cols1": {
            id1: {1: 1},
            id2: {}
        },
        "diff_cols2": {
            id1: {1: 1},
            id2: {}
        }
    }

    eq_(res, expected_result)


def get_max_sp_score():
    golden_aln = {
        '1': '-BCDEF-',
        '2': 'ABCD---'
    }

    expected = 9
    score = aa.get_max_sp_score(golden_aln)
    eq_(score, expected)


def test_calc_pairwise_score():
    golden_aln = {
        '1': '-123--4',
        '2': '12-34-5'
    }

    test_aln = {
        '1': '1-23--4-',
        '2': '12-34--5'
    }

    expected = {
        'matrix': {
            'TP': 2,
            'FP': 2,
            'TN': 2,
            'FN': 3
        },
        'wrong_cols': {
            '1': {
                0: 1,
                1: 1,
                6: 1,
                7: 1
            },
            '2': {
                0: 1,
                1: 1,
                6: 1,
                7: 1
            }
        },
        'sp_score': float(2) / 9
    }

    result = aa.calc_pairwise_score(golden_aln, '1', test_aln['1'], '2',
                                    test_aln['2'])
    eq_(result, expected)


def test_score_core_regions_3dm():
    golden_aln = {
        '1': '-123--45-',
        '2': '12-34-567'
    }

    test_aln = {
        '1': '13--4-',
        '2': '134--5'
    }

    id1 = '1'
    id2 = '2'
    result = aa.score_core_regions_3dm(test_aln, golden_aln, id1, id2,
                                       multi=False)
    expected = {
        'matrix': {
            'TP': 2,
            'FP': 2,
            'TN': 1,
            'FN': 2
        },
        'wrong_cols': {
            '1': {
                0: 1,
                4: 1,
                5: 1
            },
            '2': {
                0: 1,
                4: 1,
                5: 1
            }
        },
        'sp_score': 1
    }

    eq_(result, expected)


def test_calc_scores_3dm_complex():
    # -- prepare input data --
    with open("gold_standard_src/tests/testdata/gold_alns_datadict.txt") as a:
        g = a.read()

    # retrieve gold aln
    gold_alns = eval(g)

    # retrieve test aln
    aln_path = "gold_standard_src/tests/testdata/tautomerase_final_core.txt"
    aln_dict, strcts_order = parse_3SSP(aln_path)

    # convert to grounded
    full_seq = gold_alns["full_seq"]
    num_aln_dict, core_indexes = core_aln_to_num(
            aln_dict, full_seq, golden_ids=gold_alns["ids"])

    # -- run test --
    result = aa.calc_scores_3dm_complex(gold_alns, num_aln_dict)

