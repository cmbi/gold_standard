import os

from nose.tools import eq_, ok_

from gold_standard_src.gold_standard.html_handler import HtmlHandler
from gold_standard_src.gold_standard.num_seq import core_aln_to_num


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
    hh = HtmlHandler()
    res = hh.make_corvar(aa_aln, num_aln)

    eq_(res, expected)


def test_aln_to_html_var():
    aa_aln = {'1': 'ACDEFGKLMNOP',
              '2': 'ABC---DEF---'}
    wrong = {'1': {0: 1, 10: 1, 11: 1}, '2': {}}
    full_seq = {'1': 'ABCDEFGHIJKLMNOP',
                '2': 'ABCDEF'}
    num_aln, core_indexes = core_aln_to_num(aa_aln, full_seq)
    expected = "1      <span class=featWRONG3>A</span> b <span class=featOK>C" \
               "</span><span class=featOK>D</span><span class=featOK>E</span>" \
               "<span class=featOK>F</span><span class=featOK>G</span> hij <s" \
               "pan class=featOK>K</span><span class=featOK>L</span><span cla" \
               "ss=featOK>M</span><span class=featOK>N</span><span class=feat" \
               "WRONG3>O</span><span class=featWRONG3>P</span> \n2      <span" \
               " class=featOK>A</span>   <span class=featOK>B</span><span cla" \
               "ss=featOK>C</span>---     <span class=featOK>D</span><span cl" \
               "ass=featOK>E</span><span class=featOK>F</span>--- \n"
    hh = HtmlHandler()
    quality_data = {'num_aln': num_aln,
                    'aa_aln': aa_aln,
                    'wrong_cols': wrong,
                    'full': full_seq,
                    'core_indexes': core_indexes}
    res = hh.aln_to_html_var(quality_data)
    eq_(res, expected)

    expected_path = "gold_standard_src/tests/testdata/expected.html"
    with open(expected_path) as a:
        expected = a.read()
    expected_excerpt = expected.splitlines()[3:]
    res_path = "gold_standard_src/tests/testdata/test.html"
    hh.var = True
    hh.write_html(quality_data, "gold_standard_src/tests/testdata/test")
    ok_(os.path.exists(res_path))
    with open(res_path) as a:
        res = a.read()
    result_excerpt = res.splitlines()[3:]
    os.remove(res_path)
    eq_(expected_excerpt, result_excerpt)


def test_split_vars():
    num_aln = {
        'var': {'1': [1, 4, 9, 10]},
        'cores': {'1': [[2, 3], ['-'], ['-', 5, 6], [7, 8]]}
    }
    expected = {
        'var': {'1': [[1], [], [4], [], [9, 10]]},
        'cores': {'1': [[2, 3], ['-'], ['-', 5, 6], [7, 8]]}
    }
    hh = HtmlHandler()
    res = hh.split_vars(num_aln)
    eq_(expected, res)

    num_aln = {
        'var': {'1': [4]},
        'cores': {'1': [[1, 2, 3], [5, 6], [7, 8]]}
    }
    expected = {
        'var': {'1': [[], [4], [], []]},
        'cores': {'1': [[1, 2, 3], [5, 6], [7, 8]]}
    }
    res = hh.split_vars(num_aln)
    eq_(expected, res)


def test_make_short_var():
    l = 6
    var = "ABCDEFGH"
    hh = HtmlHandler(long_len=l)
    res = hh.make_short_var(var)
    expected = "ABC 2 FGH"
    eq_(expected, res)
    var = "ABCFGH"
    res = hh.make_short_var(var)
    expected = "ABCFGH"
    eq_(expected, res)
