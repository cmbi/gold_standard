from nose.tools import eq_

from gold_standard_src.result_processor import calc_stats


def test_calc_stats():
    matrices = {'m1': {"TP": 2, "FN": 0, "TN": 4, "FP": 6}}

    expected = {'m1': {'specificity': 0.4, 'sensitivity': 1, 'ppv': 0.25,
                       'npv': 1}}
    eq_(expected, calc_stats(matrices))
