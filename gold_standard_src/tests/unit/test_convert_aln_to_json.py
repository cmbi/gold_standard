from nose.tools import eq_

from gold_standard_src.convert_final_core_var_to_json import convert_corvar_to_json_alignment


def test_convert_corvar_to_json_alignment():
    # testcase 1
    # final_core.txt.Var:
    #   1ABCA   0   TKM l APL 0
    #   1DEFA   ltk TK- 0 TGM l
    # grounded cores:
    #   1ABCA  1, 2, 3 | 5, 6, 7
    #   1DEFA  4, 5, - | 6, 7, 8

    corvar_data = {
        "alns": {
            "cores": {
                "1ABCA": [1, 2, 3, 5, 6, 7],
                "1DEFA": [4, 5, "-", 6, 7, 8]
            },
            "var": {
                "1ABCA": [4],
                "1DEFA": [1, 2, 3]
            }
        },
        "target": "1ABCA"
    }
    score_modifiers = {'a': 1, 'u': -1, 'd': 0, 'm0': 0, 'm1': 0.1, 'm2': 0.2}

    expected = {
        "target": "1ABCA",
        "score_modifiers": score_modifiers,
        "alignments": {
            "1DEFA": {
                # first var region
                1: {"*": "u"},
                2: {"*": "u"},
                3: {"*": "u"},

                # first core
                4: {1: "a"},
                5: {2: "a"},
                # a gap

                # second core
                6: {5: "a"},
                7: {6: "a"},
                8: {7: "a"},
            }}}

    result = convert_corvar_to_json_alignment(corvar_data)
    eq_(result["alignments"], expected["alignments"])
    eq_(result, expected)
