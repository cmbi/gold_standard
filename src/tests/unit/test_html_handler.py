import os

from nose.tools import eq_, ok_
from mock import patch

from src.gold_standard.html_handler import (make_corvar, aln_to_html_var,
                                            make_short_var, write_html,
                                            split_vars)
from src.gold_standard.num_seq import core_aln_to_num


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
                  [12, 13, 14, 15, 16]]
        }
    }
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


@patch('src.gold_standard.num_seq.get_core_indexes')
def test_aln_to_html_var(mock_get_indexes):
    aa_aln = {'1': 'ACDEFGKLMNOP',
              '2': 'ABC---DEF---'}
    wrong = {'1': {0: 1, 10: 1, 11: 1}, '2': {}}
    full_seq = {'1': 'ABCDEFGHIJKLMNOP',
                '2': 'ABCDEF'}
    core_indexes = [0, 1, 6]
    mock_get_indexes.return_value = core_indexes
    final_core = 'sth'
    num_aln = core_aln_to_num(aa_aln, full_seq, final_core)[0]
    expected = "1      <span class=featWRONG3>A</span> b <span class=featOK>C" \
               "</span><span class=featOK>D</span><span class=featOK>E</span>" \
               "<span class=featOK>F</span><span class=featOK>G</span> hij <s" \
               "pan class=featOK>K</span><span class=featOK>L</span><span cla" \
               "ss=featOK>M</span><span class=featOK>N</span><span class=feat" \
               "WRONG3>O</span><span class=featWRONG3>P</span> \n2      <span" \
               " class=featOK>A</span>   <span class=featOK>B</span><span cla" \
               "ss=featOK>C</span>---     <span class=featOK>D</span><span cl" \
               "ass=featOK>E</span><span class=featOK>F</span>--- \n"
    res = aln_to_html_var(num_aln, aa_aln, wrong, full_seq, core_indexes)
    eq_(res, expected)

    expected_path = "tests/testdata/expected.html"
    with open(expected_path) as a:
        expected = a.read()
    res_path = "tests/testdata/test.html"
    write_html(aa_aln, wrong, "tests/testdata/test", var=True, var_short=False,
               num_aln=num_aln, full_seq=full_seq, core_indexes=core_indexes)

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


def test_make_short_var():
    l = 6
    var = "ABCDEFGH"
    res = make_short_var(var, l)
    expected = "ABC 2 FGH"
    eq_(expected, res)
    var = "ABCFGH"
    res = make_short_var(var, l)
    expected = "ABCFGH"
    eq_(expected, res)
