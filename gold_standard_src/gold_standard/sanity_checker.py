def check_pairwise_score(sequences, var_regs, result, id1, id2):
    res_only = [x for x in sequences[id1] + sequences[id2] + var_regs[id1] +
                var_regs[id2] if x != '-']
    if len(res_only) != sum(result['matrix'].values()):
        raise Exception("Sum of values in the confusion matrix({}) should be "
                        "equal to the total number of residues({})".format(
                            len(res_only), sum(result['matrix'].values())))
