import os

from nose.tools import ok_

from aln_quality.compare_alignments import run_comparison


def test_compare_alignments():
    aln1 = "tests/testdata/comp_aln1.fasta"
    aln2 = "tests/testdata/comp_aln2.fasta"
    full = "tests/testdata/comp_full.fasta"
    outprefix = "tests/testdata/comp_out"
    run_comparison(aln1, aln2, outprefix, full)
    out1 = outprefix + '1.html'
    out2 = outprefix + '2.html'
    ok_(os.path.exists(out1))
    os.remove(out1)
    ok_(os.path.exists(out2))
    os.remove(out2)
