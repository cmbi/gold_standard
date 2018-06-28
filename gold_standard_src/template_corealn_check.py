"""
Module responsible for filtering out templates and incorrect cores
from final_core.txt
"""
import argparse
import glob
import itertools
import json
import logging
import os
import random
import re
import subprocess
from collections import Counter
from copy import deepcopy
from numpy import isclose

from gold_standard_src.gold_standard.parsers.fasta import parse_fasta, write_fasta
from gold_standard_src.gold_standard.sbst_mat import get_blosum_score


# setup logging
logger = logging.getLogger('3DM.' + __name__)
fs = frozenset

MAFFT = 'mafft'


def parse_var_file(var_path):
    """
    Parse var file, key: seq id, value: corevar seq (string)
    :param var_path:
    :return:
    """

    with open(var_path) as a:
        inlines = a.read().splitlines()
    target_id = inlines[0].split(",")[0]
    var_aln = {}
    strcts_order = []
    for i in inlines:
        seq_id = i.split(",")[0]
        strcts_order.append(seq_id)
        sequence = i.split(",")[1]
        var_aln[seq_id] = sequence.strip()
    return var_aln, target_id, strcts_order


def calc_similarity(seq1, seq2):
    if isinstance(seq1, list):
        seq1 = "".join(seq1)
    if isinstance(seq2, list):
        seq2 = "".join(seq2)

    # max score: what would be the score if the templ seq (seq1) was aligned with itself
    max_sim_score = sum([get_blosum_score(i, i) for i in seq1.replace("-", "")])

    score = 0.0
    for i, res_1 in enumerate(seq1):
        if res_1 != "-":
            score += get_blosum_score(res_1, seq2[i])
    similarity = score / max_sim_score

    return similarity


def calc_identity(seq1, seq2):
    if isinstance(seq1, list):
        seq1 = "".join(seq1)
    if isinstance(seq2, list):
        seq2 = "".join(seq2)

    matching = 0
    for i, res_1 in enumerate(seq1):
        if res_1 == seq2[i] and res_1 != "-":
            matching += 1
    return float(matching) / len(seq1.replace("-", ""))


def convert_to_cores_by_index(incorrect_cores):
    """
    :param incorrect_cores: dict,
        {fs(seq_id1, seq_id2):
            {core_index:
                ((wif_core1, wif_core2), (mafft_core1, mafft_core2))
        }}
    :return:
        {core_index:
            {fs(seq_id1, seq_id2):
                ((wif_core1, wif_core2), (mafft_core1, mafft_core2))
        }}
    """
    new_core_dict = {}
    for seq_ids, cores_dict in incorrect_cores.iteritems():
        for core_index, core_alns in cores_dict.iteritems():
            if core_index not in new_core_dict:
                new_core_dict[core_index] = {}
            new_core_dict[core_index][seq_ids] = core_alns
    return new_core_dict


def merge_sets(sets_list):
    """
    Flatten list of sets to one level list
    [set(1,2,3), set(1,4)] => [1,2,3, 1, 4]
    """
    flat_list = []
    for i in sets_list:
        flat_list += i
    return flat_list


def fix_cores(incorrect_cores, sources, release_years, resolutions, aligned_templates):
    """
    Removes incorrect cores (== replace.them with gaps)
    """
    cores_by_index = convert_to_cores_by_index(incorrect_cores)
    new_core_alignments = deepcopy(cores_by_index)
    tmpl_no = len(sources)
    for core_index, core_alignments in cores_by_index.iteritems():
        not_corrected = core_alignments.keys()
        while not_corrected:
            # 1. first to go should be sequence templates,
            # 2. then the ones that are most abundant in the error list
            # 3. if there are only 2 templates or multiple templates are 'the most abundant'
            # than the older one is removed
            all_ids = merge_sets(not_corrected)

            ids_counter = Counter(all_ids)
            seq_id_to_remove = max(ids_counter, key=lambda x: ids_counter[x])  # pylint: disable=W0640
            freq_count = ids_counter[seq_id_to_remove]

            # if any template is very abundant in the incorrect_cores dict we will remove it
            if ids_counter.values().count(freq_count) > 1:
                # there is more than one template with the highest abundance
                frequent_ids = [seq_id for seq_id, c in ids_counter.iteritems() if c == freq_count]
                seq_id_to_remove = get_worst_template(frequent_ids, sources, release_years, resolutions)

            for seq_ids, core_alns in core_alignments.iteritems():
                if seq_id_to_remove in seq_ids:
                    core_len = len(core_alns[0][seq_id_to_remove])
                    new_aln = '-' * core_len
                    if tmpl_no == 2:
                        # if there are only 2 templates both cores are removed
                        mafft_alns = new_core_alignments[core_index][seq_ids][1]
                        new_core_alignments[core_index][seq_ids] = ([new_aln, new_aln], mafft_alns)
                        second_id = list(seq_ids.difference([seq_id_to_remove]))[0]
                        aligned_templates[seq_id_to_remove][core_index] = ''
                        aligned_templates[second_id][core_index] = ''
                    else:
                        aligned_templates[seq_id_to_remove][core_index] = new_aln
                        new_core_alignments[core_index][seq_ids][0][seq_id_to_remove] = new_aln
                    new_not_corrected = not_corrected[::]
                    for ids_pair in not_corrected:
                        if seq_id_to_remove in ids_pair:
                            new_not_corrected.remove(ids_pair)
                    not_corrected = new_not_corrected
    aligned_templates = {
        seq_id: [c for c in cores if c] for seq_id, cores in aligned_templates.iteritems()
    }
    return new_core_alignments, aligned_templates


def check_aln_coverage(aligned_cores):
    """
    Check if alignment coverage is good enough
    """
    # length of the templ seq (non-aligned)
    tmpl_seq_len = len(aligned_cores[0].replace('-', ''))

    tmpl_aln_seq = aligned_cores[0]
    aln_seq2 = aligned_cores[1]
    # tmpl residues aligned
    covered_positions = 0
    for i in range(len(tmpl_aln_seq)):
        if tmpl_aln_seq[i] != "-" and aln_seq2[i] != "-":
            covered_positions += 1

    if covered_positions >= 6:
        return True
    if float(covered_positions) / tmpl_seq_len < 0.5:
        # check coverage percentage only if length of the shorter core
        # is lower than 6
        return False
    return True


def get_newcorvar(aligned_regs):
    """
    Get core and var of the 2nd sequence based on the alignmnet
    :param aligned_regs:
    :return:
    """
    tmpl_core = aligned_regs[0].replace("-", "")
    core2 = ""
    var2 = ""
    right_var = ""
    left_var = ""
    for i in range(len(aligned_regs[0])):
        if aligned_regs[0][i] != '-':
            if right_var:
                raise RuntimeError

            core2 += aligned_regs[1][i]
        elif core2:
            right_var += aligned_regs[1][i].lower()
        else:
            left_var += aligned_regs[1][i].lower()
    if not right_var:
        right_var = "0"
    if not left_var:
        left_var = "0"
    return core2, left_var, right_var


def check_corevar(corevar1, corevar2, mafft_identity_cutoff):
    corevar1 = corevar1.split()
    corevar2 = corevar2.split()
    for i, reg_i in enumerate(corevar1):
        if i % 2 == 0:
            # this is a var region, skip it
            continue

        if corevar1[i].count('-'):
            # only take template cores if there are no gaps
            continue

        if not corevar2[i].replace('-', ''):
            # gaps-only core region
            # TODO: take the var region to the right and try to align it to the
            # core in the template (corevar1)
            # core of the template
            seq1 = corevar1[i]
            # var region to the right
            seq2 = corevar2[i + 1].upper()
            if seq2 == '0':
                continue

            aligned_regs = run_mafft_alignment(seq1, seq2)
            # if check_aln_coverage(aligned_regs) and calc_identity(aligned_regs[0], aligned_regs[1]) > mafft_identity_cutoff:
            if check_aln_coverage(aligned_regs) and calc_similarity(aligned_regs[0], aligned_regs[1]) > mafft_identity_cutoff:
                new_core, left_var, right_var = get_newcorvar(aligned_regs)
                if corevar2[i - 1] == "0":
                    corevar2[i - 1] = left_var
                elif left_var != "0":
                    corevar2[i - 1] += left_var

                corevar2[i] = new_core
                corevar2[i + 1] = right_var
    return " ".join(corevar2)


def run_mafft_alignment(core1, core2):
    """
    Align the two cores with mafft, return a list of two strings (aligned cores)
    """
    fastapath = write_fasta({'core1': core1, 'core2': core2})
    sp_args = [MAFFT, '--anysymbol', '--op', '4', fastapath]
    output = subprocess.check_output(sp_args, stderr=subprocess.PIPE)

    if os.path.exists(fastapath):
        os.remove(fastapath)
    alignment = parse_fasta(output)
    aligned_cores = [alignment["core1"], alignment["core2"]]
    return aligned_cores


def write_core_errors_to_file(aligned_templates, errors, final_core_path):
    """
    Write file listing all encountered core errors
    """
    outlines = []
    if not errors['errors']:
        outlines.append("No core errors found")
    else:
        logger.warning("Found errors in file: %s", final_core_path)
    for seq_ids, incorrect_cores in errors['errors'].iteritems():
        s = list(seq_ids)
        seq_id1 = s[0]
        seq_id2 = s[1]
        newlines = "{} {}\n{} {}".format(
            seq_id1, aligned_templates[seq_id1],
            seq_id2, aligned_templates[seq_id2],
        )
        outlines.append(newlines + "\n")
        for core_index, core_pair in incorrect_cores.iteritems():
            outlines.append("core number: {}".format(core_index + 1))
            outlines.append("WHATIF:")
            outlines.append("\n".join(core_pair[0].values()))
            outlines.append("MAFFT:")
            outlines.append("\n".join(core_pair[1].values()))
            outlines.append("whatif identity: {}\nMAFFT identity: {}\n".format(core_pair[2], core_pair[3]))
            outlines.append("")
    if not errors['full_seq_errors']:
        outlines.append("No full seq errors found")
    else:
        outlines.append("Full sequence errors:")
        outlines.append(json.dumps(errors['full_seq_errors']), indent=2)
        outlines.append("")

    output_path = final_core_path + ".core_check_errors.log"
    with open(output_path, 'w') as o:
        o.write("\n".join(outlines) + "\n")
    logger.info("Output written to: %s", output_path)


def write_core(aligned_templates, outpath, templates_order):
    """
    Write template cores to file
    """
    outlines = []
    for seq_id in templates_order:
        if seq_id in aligned_templates:
            outlines.append("{} {}".format(seq_id, " ".join(aligned_templates[seq_id].split()[1::2])))
    with open(outpath, 'w') as o:
        o.write("\n".join(outlines) + '\n')


def write_corevar(aligned_templates, outpath, templates_order):
    """
    Write template cores to file
    """
    outlines = []
    for seq_id in templates_order:
        if seq_id in aligned_templates:
            outlines.append("{}, {}".format(seq_id, aligned_templates[seq_id]))
    with open(outpath, 'w') as o:
        o.write("\n".join(outlines) + '\n')


def merge_two_cores(aligned_templates, var_index):
    """
    Will merge 2 cores neighbouring with region of index 'var_index'
    """
    for seq_id, seq in aligned_templates.iteritems():
        regions = seq.split()
        assert regions[var_index] == "0"
        new_seq = " ".join(
            regions[:var_index - 1]  + [regions[var_index - 1] + regions[var_index + 1]] + regions[var_index + 2:]
        )
        aligned_templates[seq_id] = new_seq


def merge_corvar(aligned_templates, merged=False):
    """
    If no templates have residues in a var and there are no neighbouring gaps merge neighbouring cores
    """
    var_index = -1

    # iterate over all vars except first and last - there are no 2 neighbouring
    # cores to merge there obviously
    for i in range(2, len(aligned_templates.values()[0].split()) - 2, 2):
        cannot_merge = False

        for seq_id, seq in aligned_templates.iteritems():
            regions = seq.strip().split()
            if regions[i] != '0' or regions[i - 1].endswith("-") or regions[i + 1].startswith("-"):
                cannot_merge = True
                break

        if not cannot_merge:
            var_index = i
            break


    if not cannot_merge:
        # no residues at this position, merge these cores and rerun
        # merge_corevar
        merge_two_cores(aligned_templates, var_index)
        merge_corvar(aligned_templates, merged=True)

    # did not find any cores to merge, quit function
    return merged


def run_check(corevar_path, tmpl_identity=0.4, tmpl_id="",
              mafft_identity=0.6, write_log=False, diff_check=True,
              outvar="", outfinal=""):
    """
    Checks correctness of core alignments in the final_core alignment by comparison
        to MAFFT alignments and removed the possibly incorrect cores
    """
    cutoffs = {'mafft': mafft_identity, 'tmpl': tmpl_identity}

    # parse final_core.txt.Var
    aligned_templates, target_id, strcts_order = parse_var_file(corevar_path)
    if tmpl_id and tmpl_id in strcts_order:
        target_id = tmpl_id

    full_sequences = get_full_sequences(aligned_templates)

    # check cores
    check_result = check_template_cores(
            aligned_templates, target_id, cutoffs['tmpl'], cutoffs['mafft'],
            write_log)

    filled_in = check_result["changed"]

    # merge cores with no var regions in between
    # merged = False
    tmp = deepcopy(aligned_templates)
    merged = merge_corvar(aligned_templates)

    if not (filled_in or merged):
        print "No changes to the input alignment"
        return
    if merged:
        print "Merged cores"
        assert tmp != aligned_templates

    print "Changed %d out of %d. target id: %s" % (check_result["changed"], len(aligned_templates), target_id)
    # sanity check to make sure that the changes are correct
    new_full_sequences = get_full_sequences(aligned_templates)

    if new_full_sequences != full_sequences:
        print "#### NEW ####"
        print  new_full_sequences
        print "#### OLD ####"
        print  full_sequences
        raise RuntimeError("Incorrect output: %s" % str(aligned_templates))

    if not outvar:
        # if output var not provided overwrite input var
        outvar = corevar_path
    write_corevar(check_result["new_templates"], outvar, strcts_order)
    if not outfinal and outvar.endswith(".Var"):
        outfinal = outvar[:-4]
    if outfinal:
        write_core(check_result["new_templates"], outfinal, strcts_order)


def get_full_sequences(aligned_templates):
    """
    Get full sequences from the corvar sequences
    """
    full_sequences = {}

    for seq_id, seq in aligned_templates.iteritems():
        full_seq = seq.replace(" ", "").replace("0", "").replace("-", "").upper()
        full_sequences[seq_id] = full_seq
    return full_sequences


def check_template_cores(aligned_templates, tmpl_id, tmpl_identity_cutoff=0.5,
                         mafft_identity_cutoff=0.6, write_log=False):
    """
    Checks correctness of core alignments in the final_core alignment by comparison
        to MAFFT alignments and removed the possibly incorrect cores
    """
    checked_pairs = set()
    all_pairs = (len(aligned_templates) * (len(aligned_templates) - 1)) / 2

    # if more than 50 templates than only do a 1 vs all check, not all vs all
    tmpl_cores = aligned_templates[tmpl_id]
    checks_to_run = len(aligned_templates) - 1

    seq_id1 = tmpl_id
    corevar1 = tmpl_cores
    changed = 0
    # run per-template core checks

    new_aligned_templates = {tmpl_id: tmpl_cores}
    core_seq1 = corevar1.split(" ")[1::2]

    for seq_id2, corevar2 in aligned_templates.iteritems():
        core_seq2 = corevar2.split(" ")[1::2]

        # check if this templates pair wasn't already checked
        ids_set = fs([seq_id1, seq_id2])
        if seq_id1 == seq_id2 or ids_set in checked_pairs:
            continue
        checked_pairs.add(ids_set)

        # run the check
        full_identity = calc_similarity(core_seq1, core_seq2)
        if 0.1 < full_identity < tmpl_identity_cutoff:
            # only compare highly identical templates
            new_aligned_templates[seq_id2] = corevar2
            continue

        # run core-by-core check
        newcorevar = check_corevar(corevar1, corevar2, mafft_identity_cutoff)
        if newcorevar != corevar2:
            changed += 1

        new_aligned_templates[seq_id2] = newcorevar.replace("  ", " ")

        # check if all pairs have been checked
        if len(checked_pairs) == all_pairs:
            break

    return {
        'new_templates': new_aligned_templates,
        'old_templates': aligned_templates,
        "changed": changed
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("corvar", help="corevar file")
    parser.add_argument("--tmpl_identity", help="template identity cutoff",
                        default=0.2, type=float)
    parser.add_argument("--mafft_identity", help="mafft identity cutoff",
                        default=0.2, type=float)
    parser.add_argument('--write_log', action='store_true', default=False)
    parser.add_argument('--nodiffcheck', action='store_true', default=False)

    parser.add_argument('--outvar')
    parser.add_argument('--outfinal')
    parser.add_argument('--tmpl_id')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    run_check(corevar_path=args.corvar, tmpl_identity=args.tmpl_identity,
              mafft_identity=args.mafft_identity, tmpl_id=args.tmpl_id,
              write_log=args.write_log, diff_check=not args.nodiffcheck, outvar=args.outvar,
              outfinal=args.outfinal)
