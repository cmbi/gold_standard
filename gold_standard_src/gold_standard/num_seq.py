"""
all numerical sequences are 1-based
"""
import logging
import re


_log = logging.getLogger(__name__)


class ParsingError(Exception):
    pass


def corvar_to_num(corvar_line):
    aln = {'cores': [], 'var': [], 'full': ''}
    # remove numbers and whitespaces form the corvar line
    sequence = re.sub(r'[0-9\s]', '', corvar_line)
    aln['full'] = re.sub('-', '', sequence).upper()
    count = 1
    for res_i in sequence:
        # add 1 to count because it needs to be 1-based
        if res_i.isupper():
            aln['cores'].append(count)
            count += 1
        elif res_i == '-':
            aln['cores'].append('-')
        elif res_i.islower():
            aln['var'].append(count)
            count += 1
        else:
            raise ParsingError("Incorrect character ({}) in the corvar line "
                               "({})".format(res_i, corvar_line))
    return aln


def core_aln_to_num(aln_dict, full_seq, golden_ids=None):
    """
    :param aln_path: path to the core alignment file
    :param full_seq_path: path to the full plain sequences in fasta format
    :param final_core: path to the final_core.txt file
    :return: dict with core alignment (num), and var - a dict of the positions
    in the variable regions
    """
    _log.info("Converting 3DM alignment to grounded sequences")
    aln_3dm = {"cores": {}, "var": {}}
    core_indexes = set()
    for seq_id, seq in aln_dict.iteritems():
        if golden_ids and seq_id not in golden_ids:
            continue
        try:
            aln_3dm["cores"][seq_id], core_indexes_tmp = core_to_num_seq(
                seq, full_seq[seq_id])
            core_indexes = core_indexes.union(set(core_indexes_tmp))
            aln_3dm["var"][seq_id] = get_var_pos(aln_3dm["cores"][seq_id],
                                                 full_seq[seq_id])
        except ParsingError as e:
            msg = "There was an error processing sequence %s.\nfull sequence:\n%s\ncore sequence:\n%s" % ( seq_id, full_seq[seq_id], seq)
            e.message = msg + "\n" + e.message
            raise e

    return aln_3dm, list(core_indexes)


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
        raise ParsingError("final_core file has incorrect format")
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
            raise ParsingError("Didn't find the {} core in the full "
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
    if aligned_seq.count("-") == len(aligned_seq):
        return list(aligned_seq), [0]
    # 1-based!!!
    grounded_seq = []
    start = 0
    finished = False
    prev_core = 0
    all_cores = []
    if aligned_seq.count("-") == len(aligned_seq):
        grounded_seq = list(aligned_seq)
        finished = True
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
        try:
            cores = split_core(core, full_seq[prev_core:])
        except ParsingError as e:
            _log.error("Found these cores: %s", all_cores)
            raise e
        cores = sorted(cores, key=lambda x: x["pos"])
        all_cores.extend(cores)
        for c in cores:
            grounded_seq.extend([c['pos'] + i + 1 + prev_core
                                 for i in range(len(c['seq']))])
        prev_core = cores[-1]['pos'] + prev_core + len(cores[-1]['seq'])
    # fill in the c-terminal gaps
    grounded_seq += '-' * (len(aligned_seq) - len(grounded_seq))
    new_core_indexes = get_core_indexes_from_grounded(grounded_seq)
    return grounded_seq, new_core_indexes


def get_core_indexes_from_grounded(grounded_seq):
    prev = get_first_num(grounded_seq) - 1
    core_indexes = [0]
    for i, res in enumerate(grounded_seq):
        if prev == '-' or res == '-':
            prev = res
            continue
        if prev != int(res) - 1:
            core_indexes.append(i)
        prev = int(res)
    return core_indexes


def get_first_num(grounded_seq):
    first_num = 1
    for i in grounded_seq:
        if i != '-':
            first_num = i
            return first_num
    # no cores
    return len(grounded_seq)


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

    if full_seq[add_index:].find(core) != -1 and len(core) == 1:
        return [{'pos': full_seq[add_index:].find(core) + add_index,
                 'seq': core}]

    for i in xrange(1, len(core)):
        core1_pos = full_seq[add_index:].find(core[:-i])
        if core1_pos == -1 and core.find("X") != -1:
            tmp_core = core.replace("X", "G")
            core1_pos = full_seq[add_index:].find(tmp_core[:-i])

        core2_pos = full_seq[add_index + core1_pos + len(core[:-i]):].find(
            core[-i:])

        if core1_pos != -1 and core2_pos != -1:
            # means all cores are split up, can exit the function
            new_cores = [
                {"seq": core[:-i],
                 "pos": core1_pos + add_index},
                {"seq": core[-i:],
                 "pos": core2_pos + add_index + core1_pos + len(core[:-i])}]
            break

        elif core1_pos != -1:
            # core1 was found, core2 not, so core 2 needs to be split up
            # further - call split_core on core2 now
            position = core1_pos + add_index
            new_cores.append({"seq": core[:-i],
                              "pos": position})
            add_index = position + len(core[:-i])
            new_cores.extend(split_core(core[-i:], full_seq,
                                        add_index))
            break

    _log.debug("Split up core %s in two cores: %s", core, new_cores)
    if not new_cores:
        raise ParsingError("Didn't find a way to split up the core {}\n"
                        "full sequence[{}:]: {} newcores: {}\n".format(
                            core, add_index, full_seq[add_index:], new_cores))
    return new_cores
