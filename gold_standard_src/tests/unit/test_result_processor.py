from nose.tools import assert_almost_equals, eq_

from gold_standard_src.gold_standard.result_processor import calc_stats


def test_calc_stats():
    matrices = {'m1': {"TP": 2, "FN": 0, "TN": 4, "FP": 6}}

    expected = {'m1': {'specificity': 0.4, 'sensitivity': 1, 'ppv': 0.25,
                       'npv': 1, 'mcc': 0.316227766}}
    result = calc_stats(matrices)
    eq_(len(result['m1']), len(expected['m1']))
    for key, value in expected['m1'].iteritems():
        assert_almost_equals(result['m1'][key], value)
