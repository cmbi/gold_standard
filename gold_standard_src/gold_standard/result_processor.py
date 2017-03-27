import logging
from numpy import sqrt

_log = logging.getLogger("__main__")


def calc_stats(confusion_matrices):
    _log.info("Calculating stats")
    stats = {}
    for m_id, m in confusion_matrices.iteritems():
        stats[m_id] = {'specificity': None, 'sensitivity': None,
                       'ppv': None, 'npv': None}
        if m['TN'] + m['FP'] != 0:
            stats[m_id]['specificity'] = float(m['TN']) / (m['TN'] + m['FP'])
        if m['TP'] + m['FN'] != 0:
            stats[m_id]['sensitivity'] = float(m['TP']) / (m['TP'] + m['FN'])
        if m['TP'] + m['FP'] != 0:
            stats[m_id]['ppv'] = float(m['TP']) / (m['TP'] + m['FP'])
        if m['TN'] + m['FN'] != 0:
            stats[m_id]['npv'] = float(m['TN']) / (m['TN'] + m['FN'])
        up = (m['TP'] * m['TN']) - (m['FP'] * m['FN'])
        down = (m['TP'] + m['FP']) * (m['TP'] + m['FN']) * \
            (m['TN'] + m['FP']) * (m['TN'] + m['FN'])
        stats[m_id]['mcc'] = up / sqrt(down)

    return stats


def process_results(matrices, full_matrix, sp_scores, output, tmpl_no):
    _log.info("Processing the results")
    out_txt = ""

    # FULL MATRIX #
    out_txt += "#### RESULTS ####\n"
    # sensitivity, specificity, ppv, npv
    out_txt += ' '.join(["{}: {}".format(k, v)
                         for k, v in full_matrix.iteritems()]) + '\n'
    # conf matrix rates (e.g. TP / total number of aa)
    total = float(sum(full_matrix.values()))
    rates = {key + 'r': val / total for key, val in full_matrix.iteritems()}
    out_txt += ' '.join(["%s: %.3f" % (k, v)
                         for k, v in rates.iteritems()]) + '\n'
    # FP, TP, FN, TN values
    full_stats = calc_stats({'full': full_matrix})['full']
    out_txt += ''.join(["{}: {}\n".format(k, v)
                        for k, v in full_stats.iteritems()]) + '\n'

    out_txt += 'aligned templates: {}\n'.format(tmpl_no)
    # average SP score
    sp_score = sum(sp_scores.values()) / len(sp_scores)
    out_txt += "SP score: {}\n".format(sp_score)
    _log.debug("Overall SP score: %f", sp_score)
    _log.debug("Conf matrix: %s", str(full_matrix))

    # PAIRWISE stats
    stats = calc_stats(matrices)
    for s_id, s in stats.iteritems():
        header = s_id if (s_id is str) else ' '.join(list(s_id))
        out_txt += "# {}\n".format(header)
        # sensitivity, specificity, ppv, npv
        out_txt += ' '.join(["{}: {}".format(k, v)
                             for k, v in matrices[s_id].iteritems()]) + '\n'
        # FP, TP, FN, TN values
        out_txt += ''.join(["{}: {}\n".format(k, v)
                            for k, v in s.iteritems()]) + '\n'
        # SP score
        out_txt += "SP score: {}\n".format(sp_scores[s_id])

    with open(output, 'w') as out:
        out.write(out_txt)
    _log.info("Created the output file: %s", output)
