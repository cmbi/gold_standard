import os

from nose.tools import eq_, ok_

from gold_standard_src.gold_standard.html_handler import HtmlHandler
from gold_standard_src.gold_standard.num_seq import core_aln_to_num


def test_make_corvar():
    aa_aln = {'1': "ABCDEFGHIJKLMNOP"}
    num_aln = {
        'cores': {
            '1': [
                [2],
                [6, 7, 8],
                [9, 10, 11]
            ],
        },
        'var': {
            '1': [[1], [3, 4, 5],
                  [12, 13, 14, 15, 16]]
        }
    }
    expected = {
        '1': {
            'cores':
            [
                "B",
                "FGH",
                "IJK"
            ],
            'var': [
                "A",
                "CDE",
                "LMNOP"
            ]}}
    hh = HtmlHandler()
    res = hh.make_corvar(aa_aln, num_aln)

    eq_(res, expected)
    #num_aln = {'var': {'1NDDA': [], '3UF8A': [1, 2, 3, 4, 34, 35, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192]}, 'cores': {'1NDDA': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74], '3UF8A': [5, 6, 7, 8, 9, 10, 11, 12, 13, '-', 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, '-', '-', 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, '-', '-', '-']}}
    full_seq = {
        '1NDDA':
            'MLIKVKTLTGKEIEIDIEPTDKVERIKERVEEKEGIPPQQQRLIYSGKQMNDEKTAADYKILGGSVLHLVLALR',
        '3UF8A':
            'PETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEMDSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHREQIGGSTVVTTESGLKYEDLTEGSGAEARAGQTVSVHYTGWLTDGQKFDSSKDRNDPFAFVLGGGMVIKGWDEGVQGMKVGGVRRLTIPPQLGYGARGAAGVIPPNATLVFEVELLDV'
    }
    res = hh.make_corvar(aa_aln, num_aln)
    print res

    #eq_(res, expected)


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


def test_find_residues_neighbouring_insertions():
    gold_aln_cores = {
        "1ABCA": ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10'],
        "1ABCB": ['2', '3', '4', '5', '-', '7', '8', '9', '10', '11'],
        "1ABCC": ['-', '3', '4', '-', '-', '7', '8', '9', '10', '11'],
    }

    expected = {
        "1ABCA": [],
        "1ABCB": ['5', '7'],
        "1ABCC": ['3', '4', '7'],
    }

    hh = HtmlHandler()
    result = hh.find_residues_neighbouring_insertions(gold_aln_cores)

    eq_(result, expected)


def test_get_full_seq_pos():
    hh = HtmlHandler()

    merged_corvar_seq = "vhINLKVKGQDGNEVFFRIKRSTQMRKLMNAY----SVDMNSIAFLFDGRRLRAEQTPDELEMEEGDEIDAMLH--"

    corvar_index = 0
    expected_pos = 2
    pos = hh.get_full_seq_pos(merged_corvar_seq, corvar_index)
    eq_(pos, expected_pos)

    corvar_index = 67
    expected_pos = 65
    pos = hh.get_full_seq_pos(merged_corvar_seq, corvar_index)
    eq_(pos, expected_pos)

    merged_corvar_seq = "pethINLKVSDGS-SEI"
    corvar_index = 10
    expected_pos = 13
    pos = hh.get_full_seq_pos(merged_corvar_seq, corvar_index)
    eq_(pos, expected_pos)

    merged_corvar_seq = "pethINLKVSDGS-SEIFFKIKKTTPLRRLMEAF--RQGKEMDSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHR---"
    corvar_index = 32
    expected_pos = 35
    pos = hh.get_full_seq_pos(merged_corvar_seq, corvar_index)
    eq_(pos, expected_pos)

