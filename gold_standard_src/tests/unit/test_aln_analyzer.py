from nose.tools import eq_

import gold_standard_src.gold_standard.aln_analyzer as aa


def test_score_var_regions():
    golden = {"1": [1, 2, 3, '-', '-', 4, 5, 6, '-', '-'],
              "2": ['-', '-', '-', 1, 2, 3, 4, 5, '-', '-']}
    var = [1, 2, 3, 4]
    matrix = aa.score_var_regions(golden, '1', '2', var)
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
        '1': '-BCDEF',
        '2': 'ABCD--'
    }

    expected = 9
    score = aa.get_max_sp_score(golden_aln)
    eq_(score, expected)
