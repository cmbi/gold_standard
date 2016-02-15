from nose.tools import eq_, ok_
from mock import patch, mock_open

from aln_quality_script.aln_quality import calc_alignment_quality as ca


def test_convert_var_to_aln():
    var_file = ["ID1, sdfg SFG 0 DFG abc",
                "ID2, GFG g VBF bc"
                ]
    expected = {"ID1": "SDFGSFG-DFGABC--",
                "ID2": "----GFGGVBF---BC"
                }

    eq_(ca.convert_var_to_aln(var_file), expected)


def test_parse_var_file():
    test_path = "tests/testdata/test.Var"
    var = ca.parse_var_file(test_path)
    eq_(len(var['aln']), 2)
    eq_(var['ids'], ['1HVXA', '1E43A'])


@patch('aln_quality_script.aln_quality.calc_alignment_quality.open',
       mock_open(read_data=">ID1\nA-C\nDEF\n>ID2\nGHI\n"), create=True)
@patch('aln_quality_script.aln_quality.calc_alignment_quality.os.path.exists')
def test_parse_fasta(mock_path_exists):
    mock_path_exists.return_value = True
    aln = ca.parse_fasta("path", ["ID1", "ID2"])
    expected = {"ID1": "A-CDEF", "ID2": "GHI"}
    eq_(aln, expected)


def test_calc_confusion_matrix():
    golden_aln = {'id1': "--12--345-6",
                  'id2': "--12--3456-"}
    test_aln = {'id1': "--12--345-6",
                'id2': "-12--3-456-"}

    expected = {"TP": 4, "FN": 4, "TN": 2, "FP": 2}

    m = ca.calc_confusion_matrix(golden_aln, 'id1', test_aln['id1'], 'id2',
                                 test_aln['id2'])
    eq_(m, expected)


@patch('aln_quality_script.aln_quality.calc_alignment_quality.'
       'calc_confusion_matrix')
def test_calc_confusion_matrices(mock_calc):
    m1 = {"TP": 2, "FN": 0, "TN": 1, "FP": 1}
    m2 = {"TP": 3, "FN": 4, "TN": 7, "FP": 1}
    m3 = {"TP": 1, "FN": 3, "TN": 2, "FP": 0}
    mock_calc.side_effect = [m1, m2, m3]
    fs = frozenset
    test_aln = {"id1": "seq1", "id2": "seq2", "id3": "seq3"}
    golden_alns = {
        fs(['id1', 'id2']): "gold1",
        fs(['id1', 'id3']): "gold2",
        fs(['id2', 'id3']): "gold3"
    }
    expected_all = {"TP": 6, "FN": 7, "TN": 10, "FP": 2}
    result = ca.calc_confusion_matrices(golden_alns, test_aln)
    ok_(all(m in result[0].values() for m in [m1, m2, m3]))
    eq_(result[1], expected_all)


@patch('aln_quality_script.aln_quality.calc_alignment_quality.'
       'calc_confusion_matrix_3dm')
def test_calc_confusion_matrices_3dm(mock_calc):
    m1 = {"TP": 2, "FN": 0, "TN": 1, "FP": 1}
    m2 = {"TP": 3, "FN": 4, "TN": 7, "FP": 1}
    m3 = {"TP": 1, "FN": 3, "TN": 2, "FP": 0}
    mock_calc.side_effect = [m1, m2, m3]
    fs = frozenset
    test_aln = {"cores": {"id1": "seq1", "id2": "seq2", "id3": "seq3"},
                "var": {"id1": "var", "id2": "var", "id3": "var"}
                }
    golden_alns = {
        fs(['id1', 'id2']): "gold1",
        fs(['id1', 'id3']): "gold2",
        fs(['id2', 'id3']): "gold3"
    }
    expected_all = {"TP": 6, "FN": 7, "TN": 10, "FP": 2}
    result = ca.calc_confusion_matrices_3dm(golden_alns, test_aln)
    ok_(all(m in result[0].values() for m in [m1, m2, m3]))
    eq_(result[1], expected_all)


@patch('aln_quality_script.aln_quality.calc_alignment_quality.os.path.exists')
@patch('aln_quality_script.aln_quality.calc_alignment_quality.os.listdir')
@patch('aln_quality_script.aln_quality.calc_alignment_quality.parse_var_file')
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
    eq_(expected, ca.core_to_num_seq(test_seq, full_seq))

    # TODO: reduce it to a smaller test (must be case where the cores need to
    # be split up)
    test_seq = "--KSW-SYVRSQTPLFTVGEYWSYNKLHNY-"

    full_seq = "KSWGKWYVNTTNIDGFRLDAVKHIKFSFFPDWLSYVRSQTGKPLFTVGEYWSYDINKLHN" \
               "YIMKTNGTMSLFDAPLHNKFYTASK"
    expected = ["-", "-", 1, 2, 3, "-", 34, 35, 36, 37, 38, 39, 40, 43, 44, 45,
                46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, "-"]

    eq_(expected, ca.core_to_num_seq(test_seq, full_seq))


def test_get_next_core():
    seq = "--ASDFSSSDFH-SS-"

    start = 0
    r = ca.get_next_core(seq, start)
    expected_core = "ASDFSSSDFH"
    expected_start = 2
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)

    start = 12
    r = ca.get_next_core(seq, start)
    expected_core = "SS"
    expected_start = 13
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)


def test_get_var_pos():
    num_seq = ['-', '-', 2, '-', 5, 6, '-']
    full_seq = "ABCDEFG"

    expected = [1, 3, 4, 7]
    eq_(expected, ca.get_var_pos(num_seq, full_seq))


@patch('aln_quality_script.aln_quality.calc_alignment_quality.open',
       mock_open(read_data=">ID1\n--BCD--\n>ID2\n-BC-FG-\n"), create=True)
@patch('aln_quality_script.aln_quality.calc_alignment_quality.os.path.exists')
def test_parse_3dm_aln(mock_path_exists):
    mock_path_exists.return_value = True
    full_seq = {"ID1": "ABCDEF", "ID2": "BCDEFG"}
    expected = {"cores": {"ID1": ['-', '-', 2, 3, 4, '-', '-'],
                          "ID2": ['-', 1, 2, '-', 5, 6, '-']},
                "var": {"ID1": [1, 5, 6], "ID2": [3, 4]}}
    aln = ca.parse_3dm_aln("testpath", full_seq, ["ID1", "ID2"])
    eq_(expected, aln)


def test_score_var_regions():
    golden = {"1": [1, 2, 3, '-', '-', 4, 5, 6, '-', '-'],
              "2": ['-', '-', '-', 1, 2, 3, 4, 5, '-', '-']
              }

    var = [1, 2, 3, 4]
    matrix = ca.score_var_regions(golden, '1', '2', var)
    eq_({"TN": 3, "FN": 1}, matrix)


def test_calc_confusion_matrix_3dm():
    golden = {"1": [1, 2, 3, '-', '-', 4, 5, 6, '-', '-', 7, 8],
              "2": ['-', '-', '-', 1, 2, 3, 4, 5, '-', '-', 6, '-']
              }
    var1 = [1, 2, 3, 4]
    var2 = [1, 2]
    seq1 = ['-', 4, 5, 6, '-', '-', 8]
    seq2 = [3, 4, '-', 5, 6, '-', '-']
    expected = {"TN": 6, "FN": 4, "TP": 2, "FP": 2}
    m = ca.calc_confusion_matrix_3dm(golden, '1', seq1, var1, '2', seq2, var2)
    eq_(expected, m)


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
    result = ca.core_to_num_seq_known_cores(aligned, full_seq, core_indexes)
    eq_(result, expected)


@patch('aln_quality_script.aln_quality.calc_alignment_quality.open',
       mock_open(read_data="1ABCA ABCD ABC ABC-- DFGJH"), create=True)
def test_get_core_indexes():
    result = ca.get_core_indexes('testfile')
    expected = [0, 4, 7, 12]
    eq_(result, expected)
