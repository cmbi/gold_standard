from mock import mock_open, patch
from nose.tools import eq_

from gold_standard_src.calc_alignment_quality import calculate_aln_quality


@patch('gold_standard_src.gold_standard.parsers.fasta.open', mock_open(
    read_data=">ID11A\nA-C\nDEF\n>ID12B\nABCGHI\n"), create=True)
@patch('gold_standard_src.gold_standard.num_seq.open', mock_open(
    read_data="ID11A, 0 A-C 0 DEF 0\nID12B, 0 ABC 0 GHI 0"), create=True)
@patch('gold_standard_src.gold_standard.parsers.var_file.open', mock_open(
    read_data="ID11A, 0 A-C 0 DEF 0\nID12B, 0 ABC 0 GHI 0"), create=True)
@patch('gold_standard_src.gold_standard.parsers.fasta.os.path.exists')
@patch('gold_standard_src.calc_alignment_quality.process_results')
def test_calc_alignment_quality_multi(mock_process_results, mock_path_exists):
    # TODO: make a real test

    # dummy test to check if all the mocks work
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
    expected_call = [
        {   # pairwise confusion matrices
            frozenset(['ID11A', 'ID12B']):
                {"TP": 10, "FP": 0, "FN": 0, "TN": 1}
        },
        {"TP": 10, "FP": 0, "FN": 0, "TN": 1},  # overall confusion matrix
        {frozenset(['ID11A', 'ID12B']): 1.0},  # SP scores; 1.0 is max
        'testpath'                             # output path
    ]

    expected_result = {
        'num_aln': {
            'var': {
                'ID11A': [],
                'ID12B': []
            },
            'cores': {
                'ID11A': [1, '-', 2, 3, 4, 5],
                'ID12B': [1, 2, 3, 4, 5, 6]
            }
        },
        'full': {
            'ID11A': 'ACDEF',
            'ID12B': 'ABCGHI'

        },
        'aa_aln': {
            'ID11A': 'A-CDEF',
            'ID12B': 'ABCGHI'
        },
        'wrong_cols': {
            'ID11A': {},
            'ID12B': {}
        },
        'core_indexes': None
    }

    calc_result = calculate_aln_quality(input_paths, output, input_format, multi)
    process_res_call = list(mock_process_results.call_args[0])
    eq_(process_res_call, expected_call)
    eq_(calc_result, expected_result)
