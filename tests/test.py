from nose.tools import eq_
from mock import patch, mock_open

from aln_quality.aln_quality import calc_alignment_quality as ca
from aln_quality.aln_quality import num_seq
from aln_quality.aln_quality.parsers.fasta import parse_fasta
from aln_quality.aln_quality.parsers.var_file import (parse_var_file,
                                                      convert_var_to_aln)
from aln_quality.aln_quality.num_seq import get_var_pos


def test_convert_var_to_aln():
    var_file = ["ID1, sdfg SFG 0 DFG abc",
                "ID2, GFG g VBF bc"
                ]
    expected = {"ID1": "SDFGSFG-DFGABC--",
                "ID2": "----GFGGVBF---BC"
                }

    eq_(convert_var_to_aln(var_file), expected)


def test_parse_var_file():
    test_path = "tests/testdata/test.Var"
    var = parse_var_file(test_path)
    eq_(len(var['aln']), 2)
    eq_(var['ids'], ['1HVXA', '1E43A'])


@patch('aln_quality.aln_quality.parsers.fasta.open',
       mock_open(read_data=">ID1\nA-C\nDEF\n>ID2\nGHI\n"), create=True)
@patch('aln_quality.aln_quality.parsers.fasta.os.path.exists')
def test_parse_fasta(mock_path_exists):
    mock_path_exists.return_value = True
    aln = parse_fasta("path", ["ID1", "ID2"])
    expected = {"ID1": "A-CDEF", "ID2": "GHI"}
    eq_(aln, expected)


@patch('aln_quality.aln_quality.parsers.golden.os.path.exists')
@patch('aln_quality.aln_quality.parsers.golden.os.listdir')
@patch('aln_quality.aln_quality.parsers.golden.parse_var_file')
def test_parse_golden_alns(mock_parse, mock_listdir, mock_path_exists):
    mock_parse.side_effect = [{"ids": ["id1", "id2"], "aln": "test_aln",
                               "full": {'id1': 'seq1', 'id2': 'seq2'}},
                              {"ids": ["id2", "id3"], "aln": "test_aln",
                               "full": {'id1': 'seq1', 'id2': 'seq2'}}]
    mock_path_exists.return_value = True
    mock_listdir.return_value = ["file1.Var", "file2.var",
                                 "file3.fasta", "file4.Var"]
    eq_(2, len(ca.parse_golden_alns("test_dir")[0]))


def test_core_to_num_seq():

    test_seq = "--ABC--D"
    full_seq = "SABCSD"
    expected = ['-', '-', 2, 3, 4, '-', '-', 6]
    eq_(expected, num_seq.core_to_num_seq(test_seq, full_seq))

    # TODO: reduce it to a smaller test (must be case where the cores need to
    # be split up)
    test_seq = "--KSW-SYVRSQTPLFTVGEYWSYNKLHNY-"

    full_seq = "KSWGKWYVNTTNIDGFRLDAVKHIKFSFFPDWLSYVRSQTGKPLFTVGEYWSYDINKLHN" \
               "YIMKTNGTMSLFDAPLHNKFYTASK"
    expected = ["-", "-", 1, 2, 3, "-", 34, 35, 36, 37, 38, 39, 40, 43, 44, 45,
                46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, "-"]

    eq_(expected, num_seq.core_to_num_seq(test_seq, full_seq))


def test_get_next_core():
    seq = "--ASDFSSSDFH-SS-"

    start = 0
    r = num_seq.get_next_core(seq, start)
    expected_core = "ASDFSSSDFH"
    expected_start = 2
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)

    start = 12
    r = num_seq.get_next_core(seq, start)
    expected_core = "SS"
    expected_start = 13
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)


def test_get_var_pos():
    num_seq = ['-', '-', 2, '-', 5, 6, '-']
    full_seq = "ABCDEFG"

    expected = [1, 3, 4, 7]
    eq_(expected, get_var_pos(num_seq, full_seq))


def test_score_var_regions():
    golden = {"1": [1, 2, 3, '-', '-', 4, 5, 6, '-', '-'],
              "2": ['-', '-', '-', 1, 2, 3, 4, 5, '-', '-']
              }

    var = [1, 2, 3, 4]
    matrix = ca.score_var_regions(golden, '1', '2', var)
    eq_({"TN": 3, "FN": 1}, matrix)


def test_calc_stats():
    matrices = {'m1': {"TP": 2, "FN": 0, "TN": 4, "FP": 6}}

    expected = {'m1': {'specificity': 0.4, 'sensitivity': 1, 'ppv': 0.25,
                       'npv': 1}}
    eq_(expected, ca.calc_stats(matrices))


def test_core_to_num_seq_known_cores():
    aligned = "ABCDE--FGHI--JKL"
    core_indexes = [0, 5, 11]
    full_seq = "ABCDEGFGHIGJKL"
    expected = [1, 2, 3, 4, 5, '-', '-', 7, 8, 9, 10, '-', '-', 12, 13, 14]
    result = num_seq.core_to_num_seq_known_cores(aligned, full_seq,
                                                 core_indexes)
    eq_(result, expected)


@patch('aln_quality.aln_quality.num_seq.open',
       mock_open(read_data="1ABCA ABCD ABC ABC-- DFGJH"), create=True)
def test_get_core_indexes():
    result = num_seq.get_core_indexes('testfile')
    expected = [0, 4, 7, 12]
    eq_(result, expected)
