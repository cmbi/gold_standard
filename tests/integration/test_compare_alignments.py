import os

from nose.tools import eq_, ok_

from aln_quality.compare_alignments import run_comparison


def test_compare_alignments():
    aln1 = "tests/testdata/comp_aln1.fasta"
    aln2 = "tests/testdata/comp_aln2.fasta"
    full = "tests/testdata/comp_full.fasta"
    outprefix = "tests/testdata/comp_out"
    run_comparison(aln1, aln2, outprefix, full)
    out1_path = outprefix + '1.html'
    out2_path = outprefix + '2.html'
    out1_exp_path = outprefix + '1_expected.html'
    out2_exp_path = outprefix + '2_expected.html'

    ok_(os.path.exists(out1_path))
    with open(out1_path) as a:
        out1 = a.read()
    with open(out1_exp_path) as a:
        out1_exp = a.read()
    eq_(out1, out1_exp)
    os.remove(out1_path)

    ok_(os.path.exists(out2_path))
    with open(out2_exp_path) as a:
        out2_exp = a.read()
    with open(out2_path) as a:
        out2 = a.read()
    eq_(out2, out2_exp)
    os.remove(out2_path)
