from nose.tools import eq_

from gold_standard_src.gold_standard.dict_utils import merge_dicts


def test_merge_dicts():
    dict1 = {
        'A': 1,
        'B': 2,
        'C': 3
    }
    dict2 = {
        'A': 3,
        'B': 6,
        'D': 11
    }

    expected_dict = {
        'A': 4,
        'B': 8,
        'C': 3,
        'D': 11
    }

    merged = merge_dicts(dict1, dict2)
    eq_(merged, expected_dict)
