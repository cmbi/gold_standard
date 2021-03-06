from nose.tools import eq_
from mock import mock_open, patch

from gold_standard_src.gold_standard.parsers.aln3SSP import parse_3SSP
from gold_standard_src.gold_standard.parsers.gold import parse_gold_pairwise, \
    fill_in_target, parse_gold_json
from gold_standard_src.gold_standard.parsers.fasta import parse_fasta
from gold_standard_src.gold_standard.parsers.var_file import (
    parse_var_file, convert_var_to_aln)


def test_convert_var_to_aln():
    var_file = ["ID1, sdfg SFG 0 DFG abc",
                "ID2, GFG g VBF bc"]
    expected = {"ID1": "SDFGSFG-DFGABC--",
                "ID2": "----GFGGVBF---BC"}
    eq_(convert_var_to_aln(var_file), expected)


def test_parse_var_file():
    test_path = "gold_standard_src/tests/testdata/test.Var"
    var = parse_var_file(test_path)
    eq_(len(var['alns']), 2)
    eq_(var['ids'], ['1HVXA', '1E43A'])


@patch('gold_standard_src.gold_standard.parsers.fasta.open', mock_open(
    read_data=">ID1\nA-C\nDEF\n>ID2\nGHI\n"), create=True)
@patch('gold_standard_src.gold_standard.parsers.fasta.os.path.exists')
def test_parse_fasta(mock_path_exists):
    mock_path_exists.return_value = True
    aln = parse_fasta("path", ["ID1", "ID2"])
    expected = {"ID1": "A-CDEF", "ID2": "GHI"}
    eq_(aln, expected)


@patch('gold_standard_src.gold_standard.parsers.gold.os.path.exists')
@patch('gold_standard_src.gold_standard.parsers.gold.os.listdir')
@patch('gold_standard_src.gold_standard.parsers.gold.parse_var_file')
def test_parse_golden_alns(mock_parse, mock_listdir, mock_path_exists):
    mock_parse.side_effect = [{"ids": ["id1", "id2"], "alns": "test_aln",
                               "full_seq": {'id1': 'seq1', 'id2': 'seq2'}},
                              {"ids": ["id2", "id3"], "alns": "test_aln",
                               "full_seq": {'id1': 'seq1', 'id2': 'seq2'}}]
    mock_path_exists.return_value = True
    mock_listdir.return_value = ["file1.Var", "file2.var",
                                 "file3.fasta", "file4.Var"]
    res = parse_gold_pairwise("test_dir")
    eq_(2, len(res['alns']))


@patch('gold_standard_src.gold_standard.parsers.aln3SSP.open', mock_open(
    read_data="ID01 A A-?CDEF\nID02 A GHI\n"), create=True)
@patch('gold_standard_src.gold_standard.parsers.aln3SSP.os.path.exists')
def test_aln_3SSP(mock_path_exists):
    mock_path_exists.return_value = True
    expected = {
        'ID01A': 'A--CDEF',
        'ID02A': 'GHI'
    }
    parsed = parse_3SSP('filename')
    eq_(parsed, expected)


def test_fill_in_target():
    # corvar = {
    #     "target": "1ABCA",
    #     "alns":
    #         {
    #             "cores": {
    #                 "1ABCA": ["A", "B", "C"],
    #                 "1ABCB": ["D", "E", "F"],
    #                 "1ABCC": ["G", "H", "I"]
    #             },
    #             "vars": {
    #                 "1ABCA": ["a", "b", "", "hej"],
    #                 "1ABCB": ["", "elk", "bam"],
    #                 "1ABCC": ["sth", "sth1", "sth2"]
    #             }
    #         }}

    # expected = {
    #     "target": "1ABCA",
    #     "alns":
    #         {
    #             "cores": {
    #                 "1ABCA": ["AAB", "B", "CHEJ"],
    #                 "1ABCB": ["-D-", "E", "F---"],
    #                 "1ABCC": ["-G-", "H", "I---"]
    #             },
    #             "vars": {
    #                 "1ABCA": ["", "", "", ""],
    #                 "1ABCB": ["", "elk", "bam"],
    #                 "1ABCC": ["sth", "sth1", "sth2"]
    #             }
    #         }}

    # result = fill_in_target(corvar)
    # eq_(result, expected)
    corvar = {
        "target": "1ABCA",
        "alns":
            {
                "cores": {
                    "1ABCA": [2, 4, 5],
                    "1ABCB": [1, "-", 5],
                    "1ABCC": [2, 4, 6]
                },
                "vars": {
                    "1ABCA": [1, 3, "-", 6, 7, 8],
                    "1ABCB": ["-", 2, 3, 4, "-", "-"],
                    "1ABCC": [1, 3, 5, "-", "-", "-"]
                }
            }}

    expected = {
        "target": "1ABCA",
        "alns":
            {
                "cores": {
                    "1ABCA": [1, 2, 3, 4, 5, 6, 7, 8],
                    "1ABCB": ["-", 1, "-", "-", 5, "-", "-", "-"],
                    "1ABCC": ["-", 2, "-", "4", 6, "-", "-", "-"]
                },
                "vars": {
                    "1ABCA": ["-", "-", "-", "-", "-", "-"],
                    "1ABCB": ["-", 2, 3, 4, "-", "-"],
                    "1ABCC": [1, 3, 5, "-", "-", "-"]
                }
            }}


def test_parse_gold_json():
    json_path = "gold_standard_src/tests/testdata/final_core.json"
    var_path = "gold_standard_src/tests/testdata/final_core.txt.Var"
    res = parse_gold_json(json_path, var_path)
    print res
    eq_(len(res["full_seq"]), 13)
    eq_(len(res["ids"]), 13)
    # n - 1 alignments (n equals 13)
    eq_(len(res["alns"]), 12)
