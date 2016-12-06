import logging


_log = logging.getLogger(__name__)


def check_pairwise_score_3dm(sequences, var_regs, result, id1, id2):
    res_only = [x for x in sequences[id1] + sequences[id2] + var_regs[id1] +
                var_regs[id2] if x != '-']
    _log.info("%s %s", id1, id2)
    _log.info(result['matrix'])
    if len(res_only) != sum(result['matrix'].values()):
        raise Exception("Sum of values in the confusion matrix({}) should be "
                        "equal to the total number of residues({})".format(
                            len(res_only), sum(result['matrix'].values())))


def check_pairwise_score(seq1, seq2, matrix):
    res_only = [x for x in seq1 + seq2 if x != '-']
    if len(res_only) != sum(matrix.values()):
        raise Exception("Sum of values in the confusion matrix({}) should be "
                        "equal to the total number of residues({})".format(
                            len(res_only), sum(matrix.values())))
