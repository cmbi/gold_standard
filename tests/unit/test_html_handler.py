from nose.tools import eq_

from aln_quality.html_handler import (make_corvar, aln_to_html_var, write_html,
                                      split_vars)


def test_make_corvar():
    aa_aln = {'1': "ABCDEFGHIJKLMNOP"}
    num_aln = {
        '1': {
            'cores': [
                [1, 2],
                [6, 7, 8],
                [9, 10, 11]
            ],
            'vars': [[3, 4, 5],
                     [12, 13, 14, 15, 16]
                     ]}}
    expected = {
        '1': {
            'cores': [
                "AB",
                "FGH",
                "IJK"
            ],
            'vars': [
                "CDE",
                "LMNOP"
            ]}}
    res = make_corvar(aa_aln, num_aln)
    eq_(res, expected)


def test_aln_to_html_var():
    num_aln = {
        "1": {
            "cores": [[1], [3, 4, 5, 6, 7], [11, 12, 13, 14, 15, 16]],
            "var": [[], [2], [8, 9, 10], []]
        },
        '2': {
            "cores": [[1], [2, 3, '-', '-', '-'], [4, 5, 6, '-', '-', '-']],
            "var": [[], [], [], []]
        }
    }
    aa_aln = {'1': 'ACDEFGKLMNOP',
              '2': 'ABC---DEF---'}
    wrong = {'1': {0: 1, 10: 1, 11: 1}, '2': {}}
    full_seq = {'1': 'ABCDEFGHIJKLMNOP',
                '2': 'ABCDEF'}
    expected = ""
    res = aln_to_html_var(num_aln, aa_aln, wrong, full_seq)
    # eq_(res, expected)
    print res, expected

    write_html(aa_aln, wrong, "tests/testdata/test", var=True,
               num_aln=num_aln, full_seq=full_seq)


def test_split_vars():
    num_aln = {
        'var': {'1': [1, 4, 9, 10]},
        'cores': {'1': [[2, 3], ['-'], ['-', 5, 6], [7, 8]]}
    }
    expected = {
        'var': {'1': [[1], [], [4], [], [9, 10]]},
        'cores': {'1': [[2, 3], ['-'], ['-', 5, 6], [7, 8]]}
    }
    res = split_vars(num_aln)
    eq_(expected, res)

    num_aln = {
        'var': {'1': [4]},
        'cores': {'1': [[1, 2, 3], [5, 6], [7, 8]]}
    }
    expected = {
        'var': {'1': [[], [4], [], []]},
        'cores': {'1': [[1, 2, 3], [5, 6], [7, 8]]}
    }
    res = split_vars(num_aln)
    eq_(expected, res)
