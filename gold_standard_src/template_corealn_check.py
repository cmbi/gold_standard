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

import fasta
from gold_standard.parsers.var_file import parse_var_file

# setup logging
logger = logging.getLogger('3DM.' + __name__)
fs = frozenset

MAFFT = get_path('mafft')


def calc_identity(seq1, seq2):
    matching = 0
    for i, res_1 in enumerate(seq1):
        if res_1 == seq2[i] and res_1 != "-":
            matching += 1
    avg_len = float(len(seq1.replace("-", "")) + len(seq2.replace("-", ""))) / 2
    return float(matching) / avg_len


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
    corelen1 = len(aligned_cores[0].replace('-', ''))
    corelen2 = len(aligned_cores[1].replace('-', ''))

    # check which core is the shorter one and which is the longer one
    if corelen1 <= corelen2:
        shortlen = corelen1
        longlen = corelen2
    else:
        shortlen = corelen2
        longlen = corelen1

    covered_positions = shortlen - (len(aligned_cores[0]) - longlen)
    if covered_positions >= 6:
        return True
    if float(covered_positions) / shortlen < 0.75:
        # check coverage percentage only if length of the shorter core
        # is lower than 6
        return False
    return True

def check_corevar(corevar1, corevar2, mafft_identity_cutoff, whatif_identity_cutoff,
                seq_id1, seq_id2, diff_check=True):
    for i, reg_i in enumerate(corevar1):
        if i % 2 == 0:
            # this is a var region, skip it
            continue

        if not corevar2.replace('-', ''):
            # gaps-only core region
            # TODO: take the var region to the left and try to align it to the
            # core in the template (corevar1)
            seq1 = corevar2[i + 1].upper()
            seq2 = corevar1[i]
            aligned_regs = run_mafft_alignment(seq1, seq2)


def check_cores(cores1, cores2, mafft_identity_cutoff, whatif_identity_cutoff,
                seq_id1, seq_id2, diff_check=True):
    """
    Checks if there any potential errors in the alignment of two cores
    compares original alignment identity to mafft alignment identity
    :return: dictionary of incorrect cores,
        {core index:
         ({id1: wif_core1, id2: wif_core2},
         {id1: mafft_core1, id2: mafft_core2})
        }
    """
    incorrect_cores = {}
    for i, core_i1 in enumerate(cores1):
        core_i1_nogaps = core_i1.replace('-', '')
        core_i2 = cores2[i]
        core_i2_nogaps = core_i2.replace('-', '')
        if not (core_i1_nogaps and core_i2_nogaps):
            # one of the cores contains only gaps
            continue

        # check identity - only try to align if the original identity is lower than cutoff
        whatif_identity = aliDef.calculateCoreIdentity(core_i1, core_i2)
        if whatif_identity < whatif_identity_cutoff:
            # align cores
            aligned_cores = run_mafft_alignment(core_i1_nogaps, core_i2_nogaps)
            if not aligned_cores:
                # MAFFT created an empty alignment
                continue

            # check if mafft alignment coverage is good enough
            coverage_ok = check_aln_coverage(aligned_cores)
            if not coverage_ok:
                # mafft alignment not good enough to detect incorrect cores
                continue

            # compare mafft identity with whatif identity
            mafft_identity = aliDef.calculateCoreIdentity(aligned_cores[0], aligned_cores[1])
            if mafft_identity > mafft_identity_cutoff or \
                    (diff_check and mafft_identity > whatif_identity + 0.35):
                # whatif alignment incorrect because either mafft identity was higher than the cutoff or
                # mafft idenrtity is more than 0.35 higher than whatif identity
                incorrect_cores[i] = (
                    {seq_id1: core_i1, seq_id2: core_i2},
                    {seq_id1: aligned_cores[0], seq_id2: aligned_cores[1]},
                    whatif_identity, mafft_identity)

    return incorrect_cores


def run_mafft_alignment(core1, core2):
    """
    Align the two cores with mafft, return a list of two strings (aligned cores)
    """
    fastapath = fasta.write_fasta({
        'core1': core1, 'core2': core2})
    sp_args = [MAFFT, '--op', '4', fastapath]
    output = subprocess.check_output(sp_args, stderr=subprocess.PIPE)

    if os.path.exists(fastapath):
        os.remove(fastapath)
    aligned_cores = fasta.parse_fasta_text(output).values()
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


def get_source_by_acc_string(acc):
    """
    Check if sequence source is pdb based on the sequence id
    """
    if len(acc) == 5 and acc[0].isdigit():
        return 'pdb'
    return 'other'


def write_core_file(templates_order, aligned_templates, outpath):
    """
    Write template cores to file
    """
    outlines = []
    for seq_id in templates_order:
        if seq_id in aligned_templates:
            outlines.append("{}   {}".format(seq_id, " ".join(aligned_templates[seq_id])))
    with open(outpath, 'w') as o:
        o.write("\n".join(outlines) + '\n')


def run_check(corevar_path, tmpl_identity=0.4,
              mafft_identity=0.6, whatif_identity=0.4, write_log=False, diff_check=True):
    """
    Checks correctness of core alignments in the final_core alignment by comparison
        to MAFFT alignments and removed the possibly incorrect cores
    """
    cutoffs = {'mafft': mafft_identity, 'whatif': whatif_identity, 'tmpl': tmpl_identity}

    # parse final_core.txt
    aligned_templates, templates_order = parse_corevar_file(corevar_path)

    # check cores
    check_result = check_template_cores(
            aligned_templates, global_settings, cutoffs['tmpl'], cutoffs['mafft'],
            cutoffs['whatif'], write_log, diff_check=diff_check)


def check_template_cores(aligned_templates, global_settings, tmpl_identity_cutoff=0.5,
                         mafft_identity_cutoff=0.6, whatif_identity_cutoff=0.5,
                         write_log=False, diff_check=True):
    """
    Checks correctness of core alignments in the final_core alignment by comparison
        to MAFFT alignments and removed the possibly incorrect cores
    """
    checked_pairs = set()
    all_pairs = (len(aligned_templates) * (len(aligned_templates) - 1)) / 2
    possible_core_errors = {}
    possible_full_seq_errors = set()
    counter = 0

    # if more than 50 templates than only do a 1 vs all check, not all vs all
    tmpl_id = aligned_templates.keys()[0]
    tmpl_cores = aligned_templates[tmpl_id]
    checks_to_run = len(aligned_templates) - 1

    seq_id1 = tmpl_id
    corevar1 = tmpl_cores
    # run per-template core checks
    full_seq1 = "".join(cores1).replace("0", "").upper()
    core_seq1 = corevar1.split(" ")[1::2]

    for seq_id2, corevar2 in aligned_templates.iteritems():
        # full_seq2 = "".join(cores2).replace("0", "").upper()
        core_seq2 = corevar2.split(" ")[1::2]

        # check if this templates pair wasn't already checked
        ids_set = fs([seq_id1, seq_id2])
        if seq_id1 == seq_id2 or ids_set in checked_pairs:
            continue
        checked_pairs.add(ids_set)

        # log progress
        counter += 1
        if counter % 10 == 0:
            logger.info("Running check %s out of %s", counter, checks_to_run)

        # run the check
        full_identity = calc_identity(core_seq1, core_seq2)
        if 0.1 < full_identity < tmpl_identity_cutoff:
            # only compare highly identical templates
            continue

        # run core-by-core check
        incorrect_cores = check_corevar(corevar1, corevar2, mafft_identity_cutoff, whatif_identity_cutoff,
                                        seq_id1, seq_id2, diff_check=diff_check)
        if incorrect_cores:
            possible_core_errors[ids_set] = incorrect_cores

        # check if all pairs have been checked
        if len(checked_pairs) == all_pairs:
            break

    # write lgos to file
    if write_log and possible_core_errors:
        logpath = "template_core_check.log"
        errors = {'errors': possible_core_errors, 'full_seq_errors': possible_full_seq_errors}
        write_core_errors_to_file(aligned_templates, errors, logpath)

    # if no errors were found return
    if (possible_core_errors or possible_full_seq_errors):
        print "Found {} core errors".format(len(possible_core_errors))

    return {
        'new_templates': new_aligned_templates,
        'old_templates': aligned_templates,
        'errors': possible_core_errors,
        'full_seq_errors': possible_full_seq_errors
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("corvar", help="corevar file")
    parser.add_argument("--tmpl_identity", help="template identity cutoff",
                        default=0.6, type=float)
    parser.add_argument("--mafft_identity", help="mafft identity cutoff",
                        default=0.7, type=float)
    parser.add_argument("--whatif_identity", help="whatif identity cutoff",
                        default=0.5, type=float)
    parser.add_argument('--write_log', action='store_true', default=False)
    parser.add_argument('--nodiffcheck', action='store_true', default=False)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    run_check(corevar_path=args.corevar, tmpl_identity=args.tmpl_identity,
              mafft_identity=args.mafft_identity, whatif_identity=args.whatif_identity,
              write_log=args.write_log, diff_check=not args.nodiffcheck)
