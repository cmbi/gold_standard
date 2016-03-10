from nose.tools import eq_

from aln_quality.aln_analyzer import score_var_regions, compare_pairwise


def test_score_var_regions():
    golden = {"1": [1, 2, 3, '-', '-', 4, 5, 6, '-', '-'],
              "2": ['-', '-', '-', 1, 2, 3, 4, 5, '-', '-']
              }

    var = [1, 2, 3, 4]
    matrix = score_var_regions(golden, '1', '2', var)
    eq_({"TN": 3, "FN": 1}, matrix)


def test_compare_pairwise():
    id1 = '1'
    id2 = '2'
    cores1_aln1 = "-1234--"
    cores2_aln1 = "-1234--"
    cores1_aln2 = "-123-4-"
    cores2_aln2 = "1-23-4-"
    res = compare_pairwise(id1, cores1_aln1, cores2_aln1, id2, cores1_aln2,
                           cores2_aln2)

    expected_result = {
        "diff_cols1": {
            id1: {1: 1},
            id2: {1: 1}
        },
        "diff_cols2": {
            id1: {1: 1},
            id2: {0: 1}
        }
    }

    eq_(res, expected_result)
