from nose.tools import eq_

from gold_standard_src import template_corealn_check as tcc


def test_merge_corvar():
    aligned_templates = {
        "1ABCA": "0 FDT 0 ABCD 0 EF 0 GH 0 BCDF 0",
        "1ABCB": "0 FDS a ABCE 0 EF 0 GG a ACDF 0",
        "1ABCC": "d YDT 0 ABCD 0 EF 0 GH 0 BCDF 0",
    }

    tcc.merge_corvar(aligned_templates)

    expected = {
        "1ABCA": "0 FDT 0 ABCDEFGH 0 BCDF 0",
        "1ABCB": "0 FDS a ABCEEFGG a ACDF 0",
        "1ABCC": "d YDT 0 ABCDEFGH 0 BCDF 0",
    }

    for seq_id, seq in expected.iteritems():
        if aligned_templates[seq_id] != expected[seq_id]:
            print aligned_templates[seq_id]
            print expected[seq_id]
    eq_(aligned_templates, expected)

    aligned_templates = {
        "1ABCA": "0 FDT 0 ABCD 0 EF 0 GH 0 BCDF 0",
        "1ABCB": "0 FDS a ABC- 0 EF 0 GG a ACDF 0",
        "1ABCC": "d YDT 0 ABCD 0 EF 0 GH 0 BCDF 0",
    }

    tcc.merge_corvar(aligned_templates)

    expected = {
        "1ABCA": "0 FDT 0 ABCD 0 EFGH 0 BCDF 0",
        "1ABCB": "0 FDS a ABC- 0 EFGG a ACDF 0",
        "1ABCC": "d YDT 0 ABCD 0 EFGH 0 BCDF 0",
    }

    for seq_id, seq in expected.iteritems():
        if aligned_templates[seq_id] != expected[seq_id]:
            print aligned_templates[seq_id]
            print expected[seq_id]
    eq_(aligned_templates, expected)


def test_get_newcorvar():
    aligned_regions = [
        '-ABCDEFGHJIKL---',
        'ABFGHJKLMNOPTRGH'
    ]

    # expected
    exp_left_var = 'a'
    exp_core = 'BFGHJKLMNOPT'
    exp_right_var = 'rgh'

    core, left_var, right_var = tcc.get_newcorvar(aligned_regions)
    eq_(core, exp_core)
    eq_(left_var, exp_left_var)
    eq_(right_var, exp_right_var)

    aligned_regions = [
        '-ABCD-FGHJIKL---',
        'ABFGHJKLMNOPTRGH'
    ]
    core, left_var, right_var = tcc.get_newcorvar(aligned_regions)
    eq_(core, "")
    eq_(left_var, "")
    eq_(right_var, "abfghjklmnoptrgh")
