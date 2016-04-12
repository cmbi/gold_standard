"""
all numerical sequences are 1-based
"""
import logging
import re


_log = logging.getLogger(__name__)


def corvar_to_num(corvar_line):
    aln = {'cores': [], 'var': [], 'full': ''}
    # remove numbers and whitespaces form the corvar line
    sequence = re.sub(r'[0-9\s]', '', corvar_line)
    aln['full'] = re.sub('-', '', sequence).upper()
    for i, res_i in enumerate(sequence):
        if res_i.isupper() or res_i == '-':
            # add 1 because it needs to be 1-based
            aln['cores'].append(i + 1)
        elif res_i.islower():
            aln['var'].append(i + 1)
        else:
            raise Exception("Incorrect character ({}) in the corvar line "
                            "({})".format(res_i, corvar_line))
    return aln


def core_aln_to_num(aln_dict, full_seq, core_indexes, golden_ids=None):
    """
    :param aln_path: path to the core alignment file
    :param full_seq_path: path to the full plain sequences in fasta format
    :param final_core: path to the final_core.txt file
    :return: dict with core alignment (num), and var - a dict of the positions
    in the variable regions
    """
    _log.info("Converting 3DM alignment to grounded sequences")
    aln_3dm = {"cores": {}, "var": {}}
    for seq_id, seq in aln_dict.iteritems():
        if golden_ids and seq_id not in golden_ids:
            continue
        if core_indexes:
            aln_3dm["cores"][seq_id] = core_to_num_seq_known_cores(
                seq, full_seq[seq_id], core_indexes)
        else:
            aln_3dm["cores"][seq_id] = core_to_num_seq(seq, full_seq[seq_id])
        aln_3dm["var"][seq_id] = get_var_pos(aln_3dm["cores"][seq_id],
                                             full_seq[seq_id])
    return aln_3dm


def get_var_pos(num_seq, full_seq):
    """
    Returns position missing in the core num seq
    """
    full_num = range(1, len(full_seq) + 1)
    var_pos = [i for i in full_num if i not in num_seq]
    return var_pos


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


def get_core_indexes(final_core_file):
    with open(final_core_file) as a:
        final_core = a.read().splitlines()[0].split()
    if len(final_core[0]) == 5:
        cores = final_core[1:]
    elif len(final_core[0]) == 4 and len(final_core[1]) == 1:
        cores = final_core[2:]
    else:
        raise Exception("final_core file has incorrect format")
    indexes = []
    prev_end = -1
    for i in cores:
        indexes.append(prev_end + 1)
        prev_end += len(i)
    return indexes


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
            raise Exception("Didn't find the {} core in the full "
                            "sequence".format(core))
        res_count = 0
        for res in core:
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
    prev_core = 0
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
        core_full_start = full_seq[prev_core:].find(core) + prev_core
        is_single_core = core_full_start != (-1 + prev_core)
        # the next 'core' actually consists of multiple cores, thus we need to
        # split them up
        if not is_single_core:
            cores = split_core(core, full_seq[prev_core:])
            cores = sorted(cores, key=lambda x: x["pos"])
            for c in cores:
                grounded_seq.extend([c['pos'] + i + 1 + prev_core
                                     for i in range(len(c['seq']))])
            prev_core = cores[-1]['pos'] + prev_core + len(cores[-1]['seq'])

        else:
            grounded_seq.extend([core_full_start + i + 1
                                 for i in range(len(core))])
            prev_core = core_full_start + len(core)
    # fill in the c-terminal gaps
    grounded_seq += '-' * (len(aligned_seq) - len(grounded_seq))
    return grounded_seq


def get_next_core(aligned_seq, start):
    """
    Returns index (from the beginning of the sequence) and sequence of the next
    core (index >= start)
    """
    core = ""
    core_start = 0
    it = start
    found = False
    in_core = False

    while not found and it < len(aligned_seq):
        res = aligned_seq[it]
        if not in_core and res != '-':
            in_core = True
            core += res.upper()
            core_start = it
        elif (res == '-' or res.islower()) and in_core:
            found = True
            if res.islower():
                core += res.upper()
        elif in_core:
            core += res.upper()
        it += 1
    return {"core": core, "core_start": core_start}


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
    _log.debug("Splitting up a core: %s\n full seq: %s\
            add_index: %s", core, full_seq, add_index)
    new_cores = []
    for i in xrange(1, len(core)):
        if full_seq[add_index:].find(core[:-i]) != -1 and \
                full_seq[add_index:].find(core[-i:]) != -1:
            new_cores = [core[:-i], core[-i:]]
            new_cores = [
                {"seq": core[:-i],
                 "pos": full_seq[add_index:].find(core[:-i]) + add_index},
                {"seq": core[-i:],
                 "pos": full_seq[add_index:].find(core[-i:]) + add_index}]
            break
        elif full_seq[add_index:].find(core[:-i]) != -1:
            position = full_seq[add_index:].find(core[:-i]) + add_index
            new_cores.append({"seq": core[:-i],
                              "pos": position})
            add_index = position + len(core[:-i])
            new_cores.extend(split_core(core[-i:], full_seq,
                                        add_index))
            break
    _log.debug("Split up core %s in two cores: %s", core, new_cores)
    if not new_cores:
        raise Exception("Didn't find a way to split up the core {}\
                full sequence[add_index:]: {}".format(
                    core, full_seq))
    return new_cores
