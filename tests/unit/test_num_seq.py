from nose.tools import eq_
from mock import mock_open, patch

from aln_quality.num_seq import (core_to_num_seq_known_cores, core_to_num_seq,
                                 get_core_indexes, get_next_core, get_var_pos)


def test_core_to_num_seq():

    test_seq = "--ABC--D"
    full_seq = "SABCSD"
    expected = ['-', '-', 2, 3, 4, '-', '-', 6]
    eq_(expected, core_to_num_seq(test_seq, full_seq))

    test_seq = "--KSW-SYVSQTPLFTVGEYWSYNKLHNY-"

    full_seq = "KSWGKWYVNTTNIDGFRLDAVKHIKFSFFPDWLSYVRSQTGKPLFTVGEYWSYDINKLHN" \
               "YIMKTNGTMSLFDAPLHNKFYTASK"
    expected = ["-", "-", 1, 2, 3, "-", 34, 35, 36, 38, 39, 40, 43, 44, 45,
                46, 47, 48, 49, 50, 51, 52, 53, 56, 57, 58, 59, 60, 61, "-"]

    eq_(expected, core_to_num_seq(test_seq, full_seq))


def test_get_var_pos():
    num_seq = ['-', '-', 2, '-', 5, 6, '-']
    full_seq = "ABCDEFG"

    expected = [1, 3, 4, 7]
    eq_(expected, get_var_pos(num_seq, full_seq))


def test_get_next_core():
    seq = "--ASDFSSSDFH-SS-"

    start = 0
    r = get_next_core(seq, start)
    expected_core = "ASDFSSSDFH"
    expected_start = 2
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)

    start = 12
    r = get_next_core(seq, start)
    expected_core = "SS"
    expected_start = 13
    eq_(r["core"], expected_core)
    eq_(r["core_start"], expected_start)


@patch('aln_quality.num_seq.open',
       mock_open(read_data="1ABCA ABCD ABC ABC-- DFGJH"), create=True)
def test_get_core_indexes():
    result = get_core_indexes('testfile')
    expected = [0, 4, 7, 12]
    eq_(result, expected)


def test_core_to_num_seq_known_cores():
    aligned = "ABCDE--FGHI--JKL"
    core_indexes = [0, 5, 11]
    full_seq = "ABCDEGFGHIGJKL"
    expected = [1, 2, 3, 4, 5, '-', '-', 7, 8, 9, 10, '-', '-', 12, 13, 14]
    result = core_to_num_seq_known_cores(aligned, full_seq, core_indexes)
    eq_(result, expected)
