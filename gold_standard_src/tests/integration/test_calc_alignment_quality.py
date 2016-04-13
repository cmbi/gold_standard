from mock import mock_open, patch
from nose.tools import eq_

from gold_standard_src.calc_alignment_quality import calculate_aln_quality


@patch('gold_standard_src.gold_standard.parsers.fasta.open', mock_open(
    read_data=">ID11A\nA-C\nDEF\n>ID12B\nGHI\n"), create=True)
@patch('gold_standard_src.gold_standard.num_seq.open', mock_open(
    read_data="ID11B, A-C DEF \nID12B, GHI"), create=True)
@patch('gold_standard_src.gold_standard.parsers.var_file.open', mock_open(
    read_data="ID11B, 0 A-C 0 DEF 0\nID12B, 0 GHI 0"), create=True)
@patch('gold_standard_src.gold_standard.parsers.fasta.os.path.exists')
@patch('gold_standard_src.calc_alignment_quality.process_results')
def test_calc_alignment_quality_multi(mock_process_results, mock_path_exists):
    # setup mocks
    mock_path_exists.return_value = True
    mock_process_results.return_value = 0

    # define input vars
    input_paths = {
        'gold_dir': 'testpath',
        'gold_path': 'testpath',
        'aln_path': 'testpath',
        'final_core': ''
    }
    output = 'testpath'
    input_format = "3dm"
    multi = True

    # expected call to process_results
    expected = [
        {},                                    # pairwise score matrix
        {"TP": 0, "FP": 0, "FN": 0, "TN": 0},  # confusion matrix
        {},                                    # SP scores
        'testpath'                             # output path
    ]

    calculate_aln_quality(input_paths, output, input_format, multi)
    process_res_calls = list(mock_process_results.call_args[0])
    eq_(process_res_calls, expected)
