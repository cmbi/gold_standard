from mock import mock_open, patch
from nose.tools import eq_, ok_

from gold_standard_src.calc_alignment_quality import calculate_aln_quality_simple, calculate_aln_quality_complex
from gold_standard_src.gold_standard.html_handler import HtmlHandler


@patch('gold_standard_src.gold_standard.parsers.fatcat.open', mock_open(
    read_data="ID11 A A-CDEF\nID12 B ABCGHI\nID13 C ACDFGH\n"),
       create=True)
@patch('gold_standard_src.gold_standard.num_seq.open', mock_open(
    read_data="ID11A, 0 A-C 0 DEF 0\nID12B, 0 ABC 0 GHI 0\nID13C, "
    "0 ACD 0 FGH 0\n"), create=True)
@patch('gold_standard_src.gold_standard.parsers.var_file.open', mock_open(
    read_data="ID11A, 0 A-C 0 DEF 0\nID12B, 0 ABC 0 GHI 0\nID13C, "
    "0 ACD 0 FGH 0"), create=True)
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
    expected_call = [
        {   # pairwise confusion matrices
            frozenset(['ID11A', 'ID12B']):
                {"TP": 10, "FP": 0, "FN": 0, "TN": 1},
            frozenset(['ID11A', 'ID13C']):
                {"TP": 10, "FP": 0, "FN": 0, "TN": 1},
            frozenset(['ID13C', 'ID12B']):
                {"TP": 12, "FP": 0, "FN": 0, "TN": 0}
        },
        {"TP": 32, "FP": 0, "FN": 0, "TN": 2},  # overall confusion matrix
        {
            frozenset(['ID11A', 'ID12B']): 1.0,
            frozenset(['ID11A', 'ID13C']): 1.0,
            frozenset(['ID13C', 'ID12B']): 1.0,
        },  # SP scores; 1.0 is max
        'testpath',                            # output path
        3   # number of templates
    ]

    expected_result = {
        'num_aln': {
            'var': {'ID11A': [], 'ID12B': [], 'ID13C': []},
            'cores': {
                'ID11A': [1, '-', 2, 3, 4, 5],
                'ID12B': [1, 2, 3, 4, 5, 6],
                'ID13C': [1, 2, 3, 4, 5, 6]
            }
        },
        'full': {'ID11A': 'ACDEF', 'ID12B': 'ABCGHI', 'ID13C': 'ACDFGH'},
        'aa_aln': {'ID11A': 'A-CDEF', 'ID12B': 'ABCGHI', 'ID13C': 'ACDFGH'},
        'wrong_cols': {'ID11A': {}, 'ID12B': {}, 'ID13C': {}},
        'core_indexes': [0]
    }

    calc_result = calculate_aln_quality_simple(input_paths, output, input_format, multi)
    process_res_call = list(mock_process_results.call_args[0])
    eq_(process_res_call, expected_call)
    eq_(calc_result, expected_result)


def test_calc_alignment_quality_complex_dummy():
    # testcase 1 - basic, all good
    gold_json_path = "gold_standard_src/tests/testdata/complex_scoring/dummy_gold_aln.json"

    aln_path = "gold_standard_src/tests/testdata/complex_scoring/dummy_test_aln.txt"
    paths = {
        "gold_path": gold_json_path,
        "aln_path": aln_path
    }
    output = ""
    in_format = "3SSP"
    calc_result = calculate_aln_quality_complex(paths, output, in_format, write_json=True)

    expected_result = {
        'overall_score': 1.0,
        'per_residue_scores':
        {
            '2ABC': {1: (True, 0.0), 2: (True, 1.0), 3: (True, 1.0), 4: (True, 1.0), 5: (True, 0.0)},
            '3ABC': {1: (True, 1.0), 2: (True, 1.0), 3: (True, 1.0), 4: (True, 1.0)},
         }
    }
    eq_(calc_result, expected_result)

    # testcase 2
    aln_path = "gold_standard_src/tests/testdata/complex_scoring/dummy_test_aln_v2.txt"
    paths = {
        "gold_path": gold_json_path,
        "aln_path": aln_path
    }
    calc_result = calculate_aln_quality_complex(paths, output, in_format, write_json=True)
    expected_result = {
        'overall_score': 0.7142857,
        'per_residue_scores':
        {
            '2ABC': {1: (True, 0.0), 2: (True, 1.0), 3: (True, 1.0), 4: (False, -1.0), 5: (True, 0.0)},
            '3ABC': {1: (True, 1.0), 2: (True, 1.0), 3: (True, 1.0), 4: (True, 1.0)},
        }
    }
    eq_(calc_result["per_residue_scores"], expected_result["per_residue_scores"])
    eq_(calc_result["per_residue_scores"], expected_result["per_residue_scores"])

    ok_(abs(calc_result["overall_score"] - expected_result["overall_score"]) < 0.0001)


def test_calc_alignment_quality_complex():
    # testcase 1 - basic, all good
    gold_json_path = "gold_standard_src/tests/testdata/complex_scoring/gold_tauto.json"
    # gold_corvar_path = "gold_standard_src/tests/testdata/complex_scoring/gold_tauto.txt.Var"

    aln_path = "gold_standard_src/tests/testdata/complex_scoring/tautomerase_final_core.txt"
    paths = {
        "gold_path": gold_json_path,
        "aln_path": aln_path
    }
    # gold_in = parse_gold_json(gold_path, corvar_path)
    output = "gold_standard_src/tests/testdata/complex_scoring/tautomerase_final_core.txt_out"
    in_format = "3SSP"
    calc_result = calculate_aln_quality_complex(paths, output, in_format, write_json=True)

    # create html output
    hh = HtmlHandler()
    hh.write_html(calc_result, "tmp", complex_scoring=True)

    # test case 2 - two solutions for first residue of 3ABFA - can be aligned
    # either with 1 or 2 (each gets 0.5 score)
    # gold_json_path = "gold_standard_src/tests/testdata/complex_scoring/gold_tauto_m5.json"
    aln_path = "gold_standard_src/tests/testdata/complex_scoring/tautomerase_final_core_alternative.txt"
    paths = {
        "gold_path": gold_json_path,
        "aln_path": aln_path
    }
    # gold_in = parse_gold_json(gold_path, corvar_path)
    output = "gold_standard_src/tests/testdata/complex_scoring/tautomerase_final_core.txt_out"
    calc_result2 = calculate_aln_quality_complex(paths, output, in_format, write_json=True)
    # create html output
    hh = HtmlHandler()
    hh.write_html(calc_result2, "tmp2", complex_scoring=True)

    expected_overall = 0.82642
    ok_(abs(calc_result2["overall_score"] - expected_overall) < 0.0001)

    # check if all residues in 3M21A are incorrectcly aligned
    all_incorrect = [
        not i[0] for k, i in calc_result2["per_residue_scores"]["3M21A"].iteritems()
    ]
    ok_(all_incorrect)
