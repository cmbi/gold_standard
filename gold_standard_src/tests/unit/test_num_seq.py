from nose.tools import eq_, raises
from mock import mock_open, patch

import gold_standard_src.gold_standard.num_seq as ns


def test_core_to_num_seq():

    test_seq = "--ABC--D"
    full_seq = "SABCSD"
    expected = ['-', '-', 2, 3, 4, '-', '-', 6]
    eq_(expected, ns.core_to_num_seq(test_seq, full_seq)[0])

    test_seq = "--KSW-SYVSQTPLFTVNKLHNY-"

    full_seq = "KSWGKWYVNTTNIDGFRLDAVKHIKFSFFPDWLSYVRSQTGKPLFTVYDINKLHN" \
               "YIMKTNGTMSLFDAPLHNKFYTASK"
    expected = ["-", "-", 1, 2, 3, "-", 34, 35, 36, 38, 39, 40, 43, 44, 45,
                46, 47, 51, 52, 53, 54, 55, 56, "-"]

    eq_(expected, ns.core_to_num_seq(test_seq, full_seq)[0])


def test_get_var_pos():
    num_seq = ['-', '-', 2, '-', 5, 6, '-']
    full_seq = "ABCDEFG"

    expected = [1, 3, 4, 7]
    eq_(expected, ns.get_var_pos(num_seq, full_seq))


def test_get_next_core():
    seq = "--ASDFSSSDFH-SS-"

    start = 0
    r = ns.get_next_core(seq, start)
    expected_core = "ASDFSSSDFH"
    expected_start = 2
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)

    start = 12
    r = ns.get_next_core(seq, start)
    expected_core = "SS"
    expected_start = 13
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)

    aln_seq = 'GRH--wGhwGH'
    expected_start = 5
    expected_core = "WGH"
    res = ns.get_next_core(aln_seq, 3)
    eq_(res['core'], expected_core)
    eq_(res['core_start'], expected_start)

    aln_seq = 'GRhwGH'
    expected_start = 3
    res = ns.get_next_core(aln_seq, 3)
    eq_(res['core'], expected_core)
    eq_(res['core_start'], expected_start)


@patch('gold_standard_src.gold_standard.num_seq.open',
       mock_open(read_data="1ABCA ABCD ABC ABC-- DFGJH"), create=True)
def test_get_core_indexes():
    result = ns.get_core_indexes('testfile')
    expected = [0, 4, 7, 12]
    eq_(result, expected)


def test_core_to_num_seq_known_cores():
    aligned = "ABCDE--FGHI--JKL"
    core_indexes = [0, 5, 11]
    full_seq = "ABCDEGFGHIGJKL"
    expected = [1, 2, 3, 4, 5, '-', '-', 7, 8, 9, 10, '-', '-', 12, 13, 14]
    result = ns.core_to_num_seq_known_cores(aligned, full_seq, core_indexes)
    eq_(result, expected)


def test_split_core():
    add_index = 0
    core = "DAS"
    full_seq = "SDWAS"
    new_cores = ns.split_core(core, full_seq, add_index)
    expected_cores = [
        {'pos': 1, 'seq': 'D'},
        {'pos': 3, 'seq': 'AS'}
    ]
    eq_(new_cores, expected_cores)

    core = "DASTGHMTGGG"
    full_seq = "SDWASMGTMGTTGHKLMKLMWMTGVGVG"
    expected_cores = [
        {'pos': 1, 'seq': 'D'},
        {'pos': 3, 'seq': 'AS'},
        {'pos': 11, 'seq': 'TGH'},
        {'pos': 21, 'seq': 'MTG'},
        {'pos': 25, 'seq': 'G'},
        {'pos': 27, 'seq': 'G'}
    ]
    new_cores = ns.split_core(core, full_seq, add_index)
    eq_(new_cores, expected_cores)


@raises(Exception)
def test_split_core_exception():
    core = "ASYTGHMTG"
    full_seq = "ASMGTMGTYGHKLMKLMWMTG"
    ns.split_core(core, full_seq)


def test_corvar_to_num():
    corvar_line = '0 ADFG bc ABC dc AFC 0 AGH 0'
    expected = {
        'cores': [1, 2, 3, 4, 7, 8, 9, 12, 13, 14, 15, 16, 17],
        'var': [5, 6, 10, 11],
        'full': 'ADFGBCABCDCAFCAGH'
    }
    num_seq = ns.corvar_to_num(corvar_line)
    eq_(num_seq, expected)


def test_core_aln_to_num():
    full_seq = {
        '1': 'ABCDEFGHIJKL'
    }

    aln_dict = {
        '1': 'ABDEFGHJKL'
    }

    expected = {
        'cores': {'1': [1, 2, 4, 5, 6, 7, 8, 10, 11, 12]},
        'var': {'1': [3, 9]}
    }

    num_aln = ns.core_aln_to_num(aln_dict, full_seq, core_indexes=None)[0]
    eq_(num_aln, expected)
