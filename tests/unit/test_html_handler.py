from nose.tools import eq_

from aln_quality.html_handler import make_corvar


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


def test_make_html_var_seq():
    pass
