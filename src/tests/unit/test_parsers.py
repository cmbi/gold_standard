from nose.tools import eq_
from mock import mock_open, patch

from src.gold_standard.parsers.golden import parse_golden_alns
from src.gold_standard.parsers.fasta import parse_fasta
from src.gold_standard.parsers.var_file import (parse_var_file,
                                                convert_var_to_aln)


def test_convert_var_to_aln():
    var_file = ["ID1, sdfg SFG 0 DFG abc",
                "ID2, GFG g VBF bc"]
    expected = {"ID1": "SDFGSFG-DFGABC--",
                "ID2": "----GFGGVBF---BC"}
    eq_(convert_var_to_aln(var_file), expected)


def test_parse_var_file():
    test_path = "tests/testdata/test.Var"
    var = parse_var_file(test_path)
    eq_(len(var['aln']), 2)
    eq_(var['ids'], ['1HVXA', '1E43A'])


@patch('aln_quality.parsers.fasta.open', mock_open(
    read_data=">ID1\nA-C\nDEF\n>ID2\nGHI\n"), create=True)
@patch('aln_quality.parsers.fasta.os.path.exists')
def test_parse_fasta(mock_path_exists):
    mock_path_exists.return_value = True
    aln = parse_fasta("path", ["ID1", "ID2"])
    expected = {"ID1": "A-CDEF", "ID2": "GHI"}
    eq_(aln, expected)


@patch('aln_quality.parsers.golden.os.path.exists')
@patch('aln_quality.parsers.golden.os.listdir')
@patch('aln_quality.parsers.golden.parse_var_file')
def test_parse_golden_alns(mock_parse, mock_listdir, mock_path_exists):
    mock_parse.side_effect = [{"ids": ["id1", "id2"], "aln": "test_aln",
                               "full": {'id1': 'seq1', 'id2': 'seq2'}},
                              {"ids": ["id2", "id3"], "aln": "test_aln",
                               "full": {'id1': 'seq1', 'id2': 'seq2'}}]
    mock_path_exists.return_value = True
    mock_listdir.return_value = ["file1.Var", "file2.var",
                                 "file3.fasta", "file4.Var"]
    eq_(2, len(parse_golden_alns("test_dir")[0]))
