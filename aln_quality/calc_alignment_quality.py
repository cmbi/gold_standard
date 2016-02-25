#!/usr/bin/python
import argparse
import logging
import os
import re

from custom_exceptions import CustomException
from html_handler import aln_to_html
from paths import CSS, TEMPLATE
from var_parser import convert_var_to_aln


fs = frozenset

FORMAT = '%(asctime)s:%(levelname)s:%(funcName)s - %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
_log = logging.getLogger(__name__)


def split_core(core, full_seq, add_index=0):
    """
    Recursively splits up a core into multiple cores
    :param core: core to split up
    :param full_seq: sequence in which we want find the cores
    :param add_index: only changed when called from inside this function,
        indicates by how many residues was the full_seq trimmed to adjust
        the position accordingly
    :return: list of dicts cores [{"core": string [core's aa seq],
                                   "pos": int [position in the full sequence]}]
    """
    _log.debug("Splitting up a core: {}\n full seq: {}".format(core, full_seq))
    new_cores = []
    for i in xrange(1, len(core)):
        if full_seq.find(core[:-i]) != -1 and full_seq.find(core[-i:]) != -1:
            new_cores = [core[:-i], core[-i:]]
            new_cores = [{"seq": core[:-i],
                          "pos": full_seq.find(core[:-i]) + add_index},
                         {"seq": core[-i:],
                          "pos": full_seq.find(core[-i:]) + add_index}]
            break
        elif full_seq.find(core[:-i]) != -1:
            position = full_seq.find(core[:-i])
            new_cores.append({"seq": core[:-i],
                              "pos": position + add_index})
            new_cores.extend(split_core(core[-i:], full_seq))

            break
    _log.debug("Split up core {} in two cores: {}".format(core, new_cores))
    if not new_cores:
        raise CustomException("Didn't find a way to split up the "
                              "core {}".format(core))
    return new_cores


def core_to_num_seq_known_cores(aligned_seq, full_seq, core_indexes):
    """
    convert aligned seq to grounded seq when core lengths are known
    """
    grounded = []
    for i, c_i in enumerate(core_indexes):
        if i < len(core_indexes) - 1:
            next_start = core_indexes[i + 1]
        else:
            next_start = len(aligned_seq)
        core = aligned_seq[c_i:next_start]
        degapped = re.sub('-', '', core)
        seq_index = full_seq.find(degapped)
        if seq_index < 0:
            raise CustomException("Didn't find the {} core in the full "
                                  "sequence".format(core))
        res_count = 0
        for j, res in enumerate(core):
            if res == '-':
                grounded.append('-')
            else:
                grounded.append(seq_index + res_count + 1)
                res_count += 1
    return grounded


def core_to_num_seq(aligned_seq, full_seq):
    """
    this function returns a grounded seq for a normally aligned seq
    :param aligned_seq:
    :param full_seq:

    alig_seq_a = '--------------DFMRI----------------HVRVYT--------------------'
    alig_seq_b = '---------------------LVKKIA----NIKY-------'

    full_seq_a = 'IVAFCLYKYFPFGGLQRDFMRIASTVAARGHHVRVYTQSWEG'
    full_seq_b = 'NLLHVAVLAFPFGTHAAPLLSLVKKIATEAPKVTFSFFCTTTTNDTLFSRSNEFLPNIKYY'
    """
    # 1-based!!!
    grounded_seq = []
    start = 0
    finished = False
    while not finished:
        c = get_next_core(aligned_seq, start)
        core = c["core"]
        core_aligned_start = c["core_start"]
        if core == '':
            finished = True
            continue
        # fill in the gaps to the next core
        grounded_seq += '-' * (core_aligned_start - start)
        start = core_aligned_start + len(core)
        # position of the core in the full sequence
        core_full_start = full_seq.find(core)
        is_single_core = core_full_start != -1
        # the next 'core' actually consists of multiple cores, thus we need to
        # split them up
        if not is_single_core:
            cores = split_core(core, full_seq)
            cores = sorted(cores, key=lambda x: x["pos"])
            for c in cores:
                grounded_seq.extend([c['pos'] + i + 1
                                     for i in range(len(c['seq']))])

        else:
            grounded_seq.extend([core_full_start + i + 1
                                 for i in range(len(core))])
    # fill in the c-terminal gaps
    grounded_seq += '-' * (len(aligned_seq) - len(grounded_seq))
    return grounded_seq


def get_next_core(aligned_seq, start):
    core = ""
    core_start = 0
    it = start
    found = False
    in_core = False

    while not found and it < len(aligned_seq):
        res = aligned_seq[it]
        if not in_core and res != '-':
            in_core = True
            core += res
            core_start = it
        elif res == '-' and in_core:
            found = True
        elif in_core:
            core += res
        it += 1
    return {"core": core, "core_start": core_start}


def merge_dicts(dict1, dict2):
    """
    Sum up values with the same keys
    """
    return {k: dict1.get(k, 0) + dict2.get(k, 0)
            for k in set(dict1) | set(dict2)}


def aln_seq_to_num(aln_seq):
    num_seq = []
    count = 1
    for i in aln_seq:
        if i != '-':
            num_seq.append(count)
            count += 1
        else:
            num_seq += '-'
    return num_seq


def parse_var_file(file_path):
    _log.debug("Parsing var file: {}".format(file_path))
    with open(file_path) as a:
        var_file = a.read().splitlines()
    ids = [i.split(',')[0] for i in var_file]
    aln = convert_var_to_aln(var_file)
    full = {}
    full[ids[0]] = re.sub('-', '', aln[ids[0]])
    full[ids[1]] = re.sub('-', '', aln[ids[1]])
    num_aln = {seq_id: aln_seq_to_num(seq) for seq_id, seq in aln.iteritems()}
    return {'ids': ids, 'aln': num_aln, 'full': full}


def parse_golden_alns(golden_dir):
    _log.info("Getting golden standard alignments")
    if not os.path.exists(golden_dir):
        raise CustomException("No such directory: {}".format(golden_dir))
    # get all ".Var" files in the given directory
    var_list = filter(lambda x: x.endswith(".Var"), os.listdir(golden_dir))
    _log.info("Got {} var files".format(len(var_list)))
    golden_alns = {}
    golden_ids = set()
    full_seq = {}
    for v in var_list:
        var = parse_var_file(os.path.join(golden_dir, v))
        golden_ids.update(var['ids'])
        aln_ids = fs(var['ids'])
        golden_alns[aln_ids] = var['aln']
        full_seq.update(var['full'])
    _log.info("Finished parsing .Var files")
    return golden_alns, golden_ids, full_seq


def get_var_pos(num_seq, full_seq):
    """
    Returns position missing in the core num seq
    """
    full_num = range(1, len(full_seq) + 1)
    var_pos = [i for i in full_num if i not in num_seq]
    return var_pos


def get_core_indexes(final_core_file):
    with open(final_core_file) as a:
        final_core = a.read().splitlines()[0].split()
    if len(final_core[0]) == 5:
        cores = final_core[1:]
    elif len(final_core[0]) == 4 and len(final_core[1]) == 1:
        cores = final_core[2:]
    else:
        raise CustomException("final_core file has incorrect format")
    indexes = []
    prev_end = -1
    for i in cores:
        indexes.append(prev_end + 1)
        prev_end += len(i)
    return indexes


def aln_3dm_to_num(aln_dict, full_seq, golden_ids, final_core):
    """
    :param aln_path: path to the core alignment file
    :param full_seq_path: path to the full plain sequences in fasta format
    :param final_core: path to the final_core.txt file
    :return: dict with core alignment (num), and var - a dict of the positions
        in the variable regions
    """
    _log.info("Parsing 3DM alignment")
    aln_3dm = {"cores": {}, "var": {}}
    if final_core:
        core_indexes = get_core_indexes(final_core)
    else:
        core_indexes = None
    for seq_id, seq in aln_dict.iteritems():
        if final_core:
            aln_3dm["cores"][seq_id] = core_to_num_seq_known_cores(
                seq, full_seq[seq_id], core_indexes)
        else:
            aln_3dm["cores"][seq_id] = core_to_num_seq(seq, full_seq[seq_id])
        aln_3dm["var"][seq_id] = get_var_pos(aln_3dm["cores"][seq_id],
                                             full_seq[seq_id])
    return aln_3dm


def parse_fasta(aln_path, golden_ids):
    _log.info("Parsing FASTA: {}".format(aln_path))
    if not os.path.exists(aln_path):
        raise CustomException("FASTA file doesn't exist: {}".format(aln_path))
    with open(aln_path) as a:
        aln_file = a.read().splitlines()
    aln_dict = {}
    seq_id = ""
    for l in aln_file:
        if l.startswith(">"):
            if '|' in l:
                seq_id = l.lstrip('>').split('|')[0]
            else:
                seq_id = l.lstrip('>')
            if seq_id in aln_dict.keys():
                raise CustomException("Sequence {} is duplicated".format(
                    seq_id))
            if seq_id in golden_ids:
                aln_dict[seq_id] = ""
        elif seq_id and seq_id in golden_ids:
            aln_dict[seq_id] += l.upper().rstrip("*")
    return aln_dict


def get_aligned_res(res_num, query_id, id2, golden_aln):
    aln_pos = golden_aln[query_id].index(res_num)
    return golden_aln[id2][aln_pos]


def calc_pairwise_score(golden_aln, id1, seq1, id2, seq2):
    if len(seq1) != len(seq2):
        raise CustomException("Aligned sequences {} and {} are not of the same "
                              "length")

    matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    sp_score = 0
    wrong_cols = {
        id1: {},
        id2: {}
    }

    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                sp_score += 2
                matrix["TP"] += 2
            else:
                sp_score -= 2
                matrix["FP"] += 2
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq1[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if res2_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq2[i] != '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        # if both are gaps do nothing

    # check output sanity
    res_only = filter(lambda x: x != "-", seq1 + seq2)
    if len(res_only) != sum(matrix.values()):
        raise CustomException("Sum of values in the confusion matrix({}) should"
                              " be equal to the total number of"
                              " residues({})".format(len(res_only),
                                                     sum(matrix.values())))
    sp_max = get_max_sp_score(golden_aln)
    sp_score = float(sp_score) / sp_max

    return {"matrix": matrix, "SP": sp_score, "wrong_cols": wrong_cols}


def calc_sums(aln, id1, id2):
    sum_pos = 0
    sum_neg = 0
    for i in range(len(aln[id1])):
        if ((aln[id1][i] == "-" and aln[id2][i] != "-") or
                (aln[id1][i] != "-" and aln[id2][i] == "-")):
            sum_neg += 1
        elif aln[id1][i] != "-" and aln[id2][i] != "-":
            sum_pos += 2
    return {"pos": sum_pos, "neg": sum_neg}


def score_var_regions(golden_aln, id1, id2, var1):
    matrix = {"FN": 0, "TN": 0}
    for p in var1:
        aln_res = get_aligned_res(p, id1, id2, golden_aln)
        if aln_res == '-':
            matrix["TN"] += 1
        else:
            matrix["FN"] += 1
    return matrix


def get_max_sp_score(golden_aln):
    seq1, seq2 = golden_aln.values()
    id1, id2 = golden_aln.keys()
    score = 0
    for i in range(len(golden_aln[id1])):
        if seq1[i] != '-' and seq2[i] != '-':
            score += 2
        elif seq1[i] != '-' or seq2[i] != '-':
            score += 1
    return score


def calc_pairwise_score_3dm(golden_aln, id1, seq1, var1, id2, seq2, var2):
    if len(seq1) != len(seq2):
        raise CustomException("Aligned sequences {} and {} are not of the same "
                              "length: {} and {}".format(seq1, seq2, len(seq1),
                                                         len(seq2)))
    matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    sp_score = 0
    wrong_cols = {
        id1: {},
        id2: {}
    }
    # score core regions
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                sp_score += 2
                matrix["TP"] += 2
            else:
                sp_score -= 2
                matrix["FP"] += 2
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq1[i] != '-' and seq2[i] == '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if res2_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        elif seq2[i] != '-' and seq1[i] == '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                sp_score += 1
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
                wrong_cols[id1][i] = 1
                wrong_cols[id2][i] = 1
        # if both are gaps do nothing
    # score variable regions in seq1
    var_matrix = score_var_regions(golden_aln, id1, id2, var1)
    matrix = merge_dicts(matrix, var_matrix)

    # score variable regions in seq2
    var_matrix = score_var_regions(golden_aln, id2, id1, var2)
    matrix = merge_dicts(matrix, var_matrix)

    # check output sanity
    res_only = filter(lambda x: x != "-", seq1 + seq2 + var1 + var2)
    if len(res_only) != sum(matrix.values()):
        raise CustomException("Sum of values in the confusion matrix({}) should"
                              " be equal to the total number of"
                              " residues({})".format(len(res_only),
                                                     sum(matrix.values())))

    sp_max = get_max_sp_score(golden_aln)
    sp_score = float(sp_score) / sp_max

    return {"matrix": matrix, "SP": sp_score, "wrong_cols": wrong_cols}


def calc_scores_3dm(golden_alns, test_aln):
    _log.info("Calculating confusion matrices [3DM mode]")
    matrices_dict = {}
    full_matrix = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    sp_scores = {}
    wrong_cols = {seq_id: {} for seq_id in test_aln["cores"].keys()}
    for id1, seq1 in test_aln["cores"].iteritems():
        for id2, seq2 in test_aln["cores"].iteritems():
            id_set = fs([id1, id2])
            if (id1 != id2 and id_set not in matrices_dict.keys() and
                    id_set in golden_alns.keys()):
                _log.debug("Calculating confusion matrix for sequences {} "
                           "and {}".format(id1, id2))
                scores = calc_pairwise_score_3dm(
                    golden_alns[id_set], id1, seq1, test_aln["var"][id1],
                    id2, seq2, test_aln["var"][id2])
                wrong_cols[id1] = merge_dicts(wrong_cols[id1],
                                              scores["wrong_cols"][id1])
                wrong_cols[id2] = merge_dicts(wrong_cols[id2],
                                              scores["wrong_cols"][id2])
                matrices_dict[id_set] = scores['matrix']
                # add the values to the matrix with overall scores
                full_matrix = merge_dicts(full_matrix, scores['matrix'])
                sp_scores[id_set] = scores['SP']
    return {'pairwise': matrices_dict, 'full': full_matrix, 'SP': sp_scores,
            'wrong_cols': wrong_cols}


def calc_scores(golden_alns, test_aln):
    _log.info("Calculating confusion matrices")
    matrices_dict = {}
    # full matrix holds scores summed up from all sequence pairs
    full_matrix = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
    sp_scores = {}
    wrong_cols = set()
    for id1, seq1 in test_aln.iteritems():
        for id2, seq2 in test_aln.iteritems():
            id_set = fs([id1, id2])
            # check if not the same sequence and if not already calculated
            if (id1 != id2 and id_set not in matrices_dict.keys() and
                    id_set in golden_alns.keys()):
                scores = calc_pairwise_score(golden_alns[id_set], id1, seq1,
                                             id2, seq2)
                wrong_cols[id1] = merge_dicts(wrong_cols[id1],
                                              scores["wrong_cols"][id1])
                wrong_cols[id2] = merge_dicts(wrong_cols[id2],
                                              scores["wrong_cols"][id2])
                matrices_dict[id_set] = scores['matrix']
                # add the values to the matrix with overall scores
                full_matrix = merge_dicts(full_matrix, scores['matrix'])
                sp_scores[id_set] = scores['SP']
    return {'pairwise': matrices_dict, 'full': full_matrix, 'SP': sp_scores,
            'wrong_cols': wrong_cols}


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
    return stats


def process_results(matrices, full_matrix, sp_scores, output):
    _log.info("Processing the results")
    out_txt = ""

    # FULL MATRIX #
    out_txt += "#### RESULTS ####\n"
    # sensitivity, specificity, ppv, npv
    out_txt += ' '.join(["{}: {}".format(k, v)
                         for k, v in full_matrix.iteritems()]) + '\n'
    # FP, TP, FN, TN values
    full_stats = calc_stats({'full': full_matrix})['full']
    out_txt += ''.join(["{}: {}\n".format(k, v)
                        for k, v in full_stats.iteritems()]) + '\n'
    # average SP score
    out_txt += "SP score: {}\n".format(sum(sp_scores.values()) / len(sp_scores))

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
    _log.info("Created the output file: {}".format(output))


def calculate_aln_quality(golden_dir, test_aln_path, output, in3dm, html,
                          final_core=None):
    golden_alns, golden_ids, full_seq = parse_golden_alns(golden_dir)
    if in3dm:
        _log.info("Calculating alignment quality in 3DM mode")
        aln_dict = parse_fasta(test_aln_path, golden_ids)
        num_aln_dict = aln_3dm_to_num(aln_dict, full_seq, golden_ids,
                                      final_core)
        scores = calc_scores_3dm(golden_alns, num_aln_dict)
    else:
        _log.info("Calculating alignment quality")
        aln_dict = parse_fasta(test_aln_path, golden_ids)
        num_aln_dict = {seq_id: aln_seq_to_num(seq)
                        for seq_id, seq in aln_dict.iteritems()}
        scores = calc_scores(golden_alns, num_aln_dict)
    process_results(scores['pairwise'], scores['full'], scores['SP'], output)
    if html:
        return {
            'wrong_cols': scores["wrong_cols"],
            'aln': aln_dict,
        }


def write_html(quality_data, outname):
    outtxt = aln_to_html(quality_data)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    css_full_path = os.path.join(script_dir, CSS)
    with open(os.path.join(script_dir, TEMPLATE)) as a:
        template_fmt = a.read()
    with open(outname, 'w') as out:
        out.write(template_fmt.format(css_full_path, outtxt))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates quality of the"
                                                 " multiple sequence alignment"
                                                 " based on the provided golden"
                                                 " standard pairwise "
                                                 " alignments")
    parser.add_argument("golden_dir")
    parser.add_argument("test_aln_path")
    parser.add_argument("output")
    parser.add_argument("--html", action="store_true")
    parser.add_argument("--in3dm", default=False, action="store_true")
    parser.add_argument("-d", "--debug", default=False, action="store_true")
    parser.add_argument("--final_core", help="final core file")

    args = parser.parse_args()

    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)

    try:
        quality_data = calculate_aln_quality(args.golden_dir,
                                             args.test_aln_path, args.output,
                                             args.in3dm, args.html,
                                             args.final_core)
        if args.html:
            write_html(quality_data, args.output + ".html")
    except CustomException as e:
        _log.error("{}".format(e.message))
        exit(1)
