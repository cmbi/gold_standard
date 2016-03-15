import os

from nose.tools import eq_, ok_

from aln_quality.html_handler import (make_corvar, aln_to_html_var, write_html,
                                      split_vars)


def test_make_corvar():
    aa_aln = {'1': "ABCDEFGHIJKLMNOP"}
    num_aln = {
        'cores': {
            '1': [
                [1, 2],
                [6, 7, 8],
                [9, 10, 11]
            ],
        },
        'var': {
            '1': [[3, 4, 5],
                  [12, 13, 14, 15, 16]
                  ]}}
    expected = {
        '1': {
            'cores':
            [
                "AB",
                "FGH",
                "IJK"
            ],
            'var': [
                "CDE",
                "LMNOP"
            ]}}
    res = make_corvar(aa_aln, num_aln)

    eq_(res, expected)


def test_aln_to_html_var():
    num_aln = {
        'cores': {
            '1': [1, 3, 4, 5, 6, 7, 11, 12, 13, 14, 15, 16],
            '2': [1, 2, 3, '-', '-', '-', 4, 5, 6, '-', '-', '-']
        },
        'var': {
            '1': [2, 9, 10],
            '2': []
        }
    }
    aa_aln = {'1': 'ACDEFGKLMNOP',
              '2': 'ABC---DEF---'}
    wrong = {'1': {0: 1, 10: 1, 11: 1}, '2': {}}
    full_seq = {'1': 'ABCDEFGHIJKLMNOP',
                '2': 'ABCDEF'}
    expected = "1      <span class=featWRONG3>A</span><span class=featOK>C</" \
               "span><span class=featOK>D</span>  <span class=featWRONG3>E</" \
               "span><span class=featOK>F</span><span class=featOK>G</span> " \
               " \n" \
               "2      <span class=featOK>A</span><span class=featOK>B</span" \
               "><span class=featOK>C</span>  ---  \n"
    core_indexes = [0, 3, 6]
    res = aln_to_html_var(num_aln, aa_aln, wrong, full_seq, core_indexes)
    eq_(res, expected)

    expected_path = "tests/testdata/expected.html"
    with open(expected_path) as a:
        expected = a.read()
    res_path = "tests/testdata/test.html"
    write_html(aa_aln, wrong, "tests/testdata/test", var=True, num_aln=num_aln,
               full_seq=full_seq)

    ok_(os.path.exists(res_path))
    with open(res_path) as a:
        res = a.read()
    os.remove(res_path)
    eq_(expected, res)


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
