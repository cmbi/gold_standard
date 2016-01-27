#!/usr/bin/python
import argparse
import logging
import os
import re

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
            new_add_index = position + len(core[:-i])
            remaining_seq = full_seq[new_add_index:]
            new_cores.extend(split_core(core[-i:], remaining_seq,
                                        new_add_index))
            break
    _log.debug("Split up core {} in two cores: {}".format(core, new_cores))
    if not new_cores:
        _log.error("Didn't find a way to split up the core {}".format(core))
        raise Exception
    print new_cores, full_seq, add_index
    return new_cores


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
            print cores
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


def convert_var_to_aln(var_file):
    id1 = var_file[0].split(',')[0]
    id2 = var_file[1].split(',')[0]
    var1 = var_file[0].split(',')[1]
    var2 = var_file[1].split(',')[1]
    # filter out the '0' core separators and split in segments
    var1 = re.sub('0', '', var1).split()
    var2 = re.sub('0', '', var2).split()
    aln_seq1 = ""
    aln_seq2 = ""
    i = 0
    j = 0
    while i < len(var1) or j < len(var2):
        if i < len(var1) and var1[i].islower():
            # var region in sequence 1 -> insert gaps in sequence 2
            aln_seq1 += var1[i].upper()
            aln_seq2 += "-" * len(var1[i])
            i += 1
        elif j < len(var2) and var2[j].islower():
            # var region in sequence 2 -> insert gaps in sequence 1
            aln_seq2 += var2[j].upper()
            aln_seq1 += "-" * len(var2[j])
            j += 1
        elif (i < len(var1) and j < len(var2) and
                var1[i].isupper() and var2[j].isupper()):
            # both segments are cores
            if len(var2[j]) != len(var1[i]):
                _log.error("Core {} and {} have different lengths".format(
                    var1[i], var2[j]))
                raise Exception
            aln_seq2 += var2[j]
            aln_seq1 += var1[i]
            i += 1
            j += 1
        else:
            _log.error("Incorrect Var file: check number of cores")
            raise Exception
    return {id1: aln_seq1, id2: aln_seq2}


def parse_var_file(file_path):
    _log.debug("Parsing var file: {}".format(file_path))
    with open(file_path) as a:
        var_file = a.read().splitlines()
    ids = [i.split(',')[0] for i in var_file]
    aln = convert_var_to_aln(var_file)
    num_aln = {seq_id: aln_seq_to_num(seq) for seq_id, seq in aln.iteritems()}
    return {'ids': ids, 'aln': num_aln}


def get_golden_alns(golden_dir):
    _log.info("Getting golden standard alignments")
    if not os.path.exists(golden_dir):
        _log.error("No such directory: {}".format(golden_dir))
        raise Exception
    # get all ".Var" files in the given directory
    var_list = filter(lambda x: x.endswith(".Var"), os.listdir(golden_dir))
    _log.info("Got {} var files".format(len(var_list)))
    golden_alns = {}
    golden_ids = set()
    for v in var_list:
        var = parse_var_file(os.path.join(golden_dir, v))
        golden_ids.update(var['ids'])
        aln_ids = fs(var['ids'])
        golden_alns[aln_ids] = var['aln']
    _log.info("Finished parsing .Var files")
    return golden_alns, golden_ids


def get_var_pos(num_seq, full_seq):
    """
    Returns position missing in the core num seq
    """
    full_num = range(1, len(full_seq) + 1)
    var_pos = [i for i in full_num if i not in num_seq]
    return var_pos


def parse_3dm_aln(aln_path, full_seq, golden_ids):
    """
    :param aln_path: path to the core alignment file
    :param full_seq_path: path to the full plain sequences in fasta format
    :return: dict with core alignment (num), and var - a dict of the positions
        in the variable regions
    """
    _log.info("Parsing 3DM alignment")
    aln_3dm = {"cores": {}, "var": {}}
    aln_fasta = parse_fasta(aln_path, golden_ids)
    for seq_id, seq in aln_fasta.iteritems():
        aln_3dm["cores"][seq_id] = core_to_num_seq(seq, full_seq[seq_id])
        aln_3dm["var"][seq_id] = get_var_pos(aln_3dm["cores"][seq_id],
                                             full_seq[seq_id])
    return aln_3dm


def parse_fasta(aln_path, golden_ids):
    _log.info("Parsing FASTA: {}".format(aln_path))
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
                _log.error("Sequence {} is duplicated".format(seq_id))
                raise Exception
            if seq_id in golden_ids:
                aln_dict[seq_id] = ""
        elif seq_id and seq_id in golden_ids:
            aln_dict[seq_id] += l
    return aln_dict


def get_aligned_res(res_num, query_id, id2, golden_aln):
    aln_pos = golden_aln[query_id].index(res_num)
    return golden_aln[id2][aln_pos]


def calc_confusion_matrix(golden_aln, id1, seq1, id2, seq2):
    if len(seq1) != len(seq2):
        _log.error("Aligned sequences {} and {} are not of the same length")
        raise Exception

    matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}

    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                matrix["TP"] += 1
            else:
                matrix["FP"] += 1
        elif seq1[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if res2_gold == '-':
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
        elif seq2[i] != '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
        # if both are gaps do nothing
    return matrix


def score_var_regions(golden_aln, id1, id2, var1):
    matrix = {"FN": 0, "TN": 0}
    for p in var1:
        aln_res = get_aligned_res(p, id1, id2, golden_aln)
        if aln_res == '-':
            matrix["TN"] += 1
        else:
            matrix["FN"] += 1
    return matrix


def calc_confusion_matrix_3dm(golden_aln, id1, seq1, var1, id2, seq2, var2):
    if len(seq1) != len(seq2):
        _log.error("Aligned sequences {} and {} are not of the same length")
        raise Exception

    matrix = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    # score core regions
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if seq2[i] == res2_gold:
                matrix["TP"] += 1
            else:
                matrix["FP"] += 1
        elif seq1[i] != '-':
            res2_gold = get_aligned_res(seq1[i], id1, id2, golden_aln)
            if res2_gold == '-':
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
        elif seq2[i] != '-':
            res1_gold = get_aligned_res(seq2[i], id2, id1, golden_aln)
            if res1_gold == '-':
                matrix["TN"] += 1
            else:
                matrix["FN"] += 1
        # if both are gaps do nothing

    # score variable regions in seq1
    var_matrix = score_var_regions(golden_aln, id1, id2, var1)
    matrix = merge_dicts(matrix, var_matrix)

    # score variable regions in seq2
    var_matrix = score_var_regions(golden_aln, id2, id1, var2)
    matrix = merge_dicts(matrix, var_matrix)
    return matrix


def calc_confusion_matrices_3dm(golden_alns, test_aln):
    _log.info("Calculating confusion matrices [3DM mode]")
    matrices_dict = {'all': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}}
    for id1, seq1 in test_aln["cores"].iteritems():
        for id2, seq2 in test_aln["cores"].iteritems():
            id_set = fs([id1, id2])
            if (id1 != id2 and id_set not in matrices_dict.keys() and
                    id_set in golden_alns.keys()):
                _log.debug("Calculating confusion matrix for sequences {} "
                           "and {}".format(id1, id2))
                matrices_dict[id_set] = calc_confusion_matrix_3dm(
                    golden_alns[id_set], id1, seq1, test_aln["var"][id1],
                    id2, seq2, test_aln["var"][id2])
                # add the values to the matrix with overall scores
                matrices_dict['all'] = merge_dicts(matrices_dict['all'],
                                                   matrices_dict[id_set])
    return matrices_dict


def calc_confusion_matrices(golden_alns, test_aln):
    _log.info("Calculating confusion matrices")
    # one of the matrices ('all') holds scores summed up from all sequence pairs
    matrices_dict = {'all': {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}}
    for id1, seq1 in test_aln.iteritems():
        for id2, seq2 in test_aln.iteritems():
            id_set = fs([id1, id2])
            # check if not the same sequence and if not already calculated
            if (id1 != id2 and id_set not in matrices_dict.keys() and
                    id_set in golden_alns.keys()):
                matrices_dict[id_set] = calc_confusion_matrix(
                        golden_alns[id_set], id1, seq1, id2, seq2)
                # add the values to the matrix with overall scores
                matrices_dict['all'] = merge_dicts(matrices_dict['all'],
                                                   matrices_dict[id_set])
    return matrices_dict


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


def process_results(matrices, output):
    _log.info("Processing the results")
    stats = calc_stats(matrices)
    with open(output, 'w') as out:
        for s_id, s in stats.iteritems():
            newlines = ''.join(["{}: {}\n".format(k, v)
                                for k, v in s.iteritems()])
            out.write("# {}\n".format(' '.join(list(s_id))))
            out.write(' '.join(["{}: {}".format(k, v)
                                for k, v in matrices[s_id].iteritems()]) + '\n')
            out.write(newlines)


def calculate_aln_quality(golden_dir, test_aln_path, output, in3dm,
                          full_seq_path):
    golden_alns, golden_ids = get_golden_alns(golden_dir)
    if in3dm:
        _log.info("Calculating alignment quality in 3DM mode")
        full_seq = parse_fasta(full_seq_path, golden_ids)
        num_aln_dict = parse_3dm_aln(test_aln_path, full_seq, golden_ids)
        confusion_matrices = calc_confusion_matrices_3dm(golden_alns,
                                                         num_aln_dict)
    else:
        _log.info("Calculating alignment quality")
        aln_dict = parse_fasta(test_aln_path, golden_ids)
        num_aln_dict = {seq_id: aln_seq_to_num(seq)
                        for seq_id, seq in aln_dict.iteritems()}
        confusion_matrices = calc_confusion_matrices(golden_alns, num_aln_dict)
    process_results(confusion_matrices, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates quality of the"
                                                 " multiple sequence alignment"
                                                 " based on the provided golden"
                                                 " standard pairwise "
                                                 " alignments")
    parser.add_argument("golden_dir")
    parser.add_argument("test_aln_path")
    parser.add_argument("output")
    parser.add_argument("--in3dm", default=False, action="store_true")
    parser.add_argument("-f", "--full_seq")
    parser.add_argument("-d", "--debug", default=False, action="store_true")

    args = parser.parse_args()

    # change logging level in debug mode
    if args.debug:
        _log.setLevel(logging.DEBUG)

    # check arguments
    if args.in3dm and not args.full_seq:
        parser.error("In the 3DM mode you need to provide the path to the full "
                     "sequences file (--full_seq)")
    try:
        calculate_aln_quality(args.golden_dir, args.test_aln_path, args.output,
                              args.in3dm, args.full_seq)
    except Exception as e:
        _log.error("Exiting: {}".format(e.message))
        exit(1)
