import os

from nose.tools import eq_, ok_

from compare_alignments import run_comparison


def test_compare_alignments_dummy():
    aln1 = "gold_standard_src/tests/testdata/comp_aln1.fasta"
    aln2 = "gold_standard_src/tests/testdata/comp_aln2.fasta"
    full = "gold_standard_src/tests/testdata/comp_full.fasta"
    outprefix = "gold_standard_src/tests/testdata/comp_out"
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
    result_excerpt = out1.splitlines()[3:]
    expected_excerpt = out1_exp.splitlines()[3:]
    eq_(result_excerpt, expected_excerpt)
    os.remove(out1_path)

    ok_(os.path.exists(out2_path))
    with open(out2_exp_path) as a:
        out2_exp = a.read()
    with open(out2_path) as a:
        out2 = a.read()
    result_excerpt = out2.splitlines()[3:]
    expected_excerpt = out2_exp.splitlines()[3:]
    eq_(result_excerpt, expected_excerpt)
    os.remove(out2_path)


def test_compare_alignments_real():
    aln1 = "gold_standard_src/tests/testdata/filtered_amylase_mafft.fasta"
    aln2 = "gold_standard_src/tests/testdata/filtered_amylase_muscle.fasta"
    full = "gold_standard_src/tests/testdata/filtered_amylase_plain.fasta"
    outprefix = "gold_standard_src/tests/testdata/comp_mafft_muscle_out"
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
    result_excerpt = out1.splitlines()[3:]
    expected_excerpt = out1_exp.splitlines()[3:]
    eq_(result_excerpt, expected_excerpt)
    os.remove(out1_path)

    ok_(os.path.exists(out2_path))
    with open(out2_exp_path) as a:
        out2_exp = a.read()
    with open(out2_path) as a:
        out2 = a.read()
    result_excerpt = out2.splitlines()[3:]
    expected_excerpt = out2_exp.splitlines()[3:]
    eq_(result_excerpt, expected_excerpt)
    os.remove(out2_path)
