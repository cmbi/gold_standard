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


def filter_templates(template_list, global_settings, family_id):
    """
    This module analyses the templates blast hits too identify very similar templates
    These need to be removed because they might be structural problems or different prot states
    :return:
    """
    logger.debug("\nTemplate list: %s", template_list)
    # calculate structural core identity for all pairs
    template_core_id_dict = {}
    template_aa_count_dict = {}
    final_core_file = '%sfamilyID_%s/final_core.txt' % (global_settings['coredir'], family_id)
    template_corseq_dict = superP.parse_core_file(final_core_file, seq_format='core_seq')[0]
    for template, core_seq in template_corseq_dict.iteritems():
        template_core_id_dict[template] = {}

        template_aa_count_dict[template] = len(re.sub('[^A-Z]', '', core_seq))

        for aln_template, aln_core_seq in template_corseq_dict.iteritems():
            core_identity = round(aliDef.calculateCoreIdentity(core_seq, aln_core_seq) * 100, 3)
            template_core_id_dict[template][aln_template] = core_identity

    # parse the blast results of each template in the other template
    template_blast_sim_dict = {}
    blast_db_name = 'pdb'
    for template in template_list:
        subfam_specific_blastdir = '%s%s/%s/' % (global_settings['blastresultdir'], family_id, template)
        blast_file = '%spdb_blast_out' % subfam_specific_blastdir
        r_dict = blasts.find_prot_in_blast_result(blast_file, blast_db_name, template_list)
        template_blast_sim_dict[template] = r_dict

    # are any similarities higher than 98%
    high_sim_count = 0
    suspicious_template_list = []
    high_pairs_list = []

    # do not use itertools.combinations here since the blast hit data is direction dependent
    for template, compare_template in itertools.product(template_list, template_list):
        if template == compare_template or compare_template not in template_blast_sim_dict[template]:
            # skip itself, should be 100%
            continue

        similarity = template_blast_sim_dict[template][compare_template]
        core_identity = template_core_id_dict[template][compare_template]
        difference = abs(similarity - core_identity)

        logging.debug('%s vs %s similarity: blast: %.2f, final core: %.2f, diff: %.2f',
                      template, compare_template, similarity, core_identity, difference)

        # requierements for finding suspicious templates, key: blast similarity,
        # value: blast similarity vs core identity difference
        similarity_diffs = {95: 10, 90: 25, 80: 30, 60: 35, 50: 40}
        for blast_sim_cutoff, diff_cutoff in similarity_diffs.iteritems():
            if similarity > blast_sim_cutoff and difference > diff_cutoff:
                high_sim_count += 1
                logging.debug('Suspicious template pair: %s-%s\t%.2f\t%.2f\tdiff: %.2f',
                              template, compare_template, similarity, core_identity, difference)

                # Add both the templates to the suspicious list (not a unique list)
                suspicious_template_list.append(template)
                suspicious_template_list.append(compare_template)
                high_pairs_list.append([template, compare_template])
                break

    if not high_sim_count:
        logging.debug('No templates have a similarity >90 and a difference >25, no template filtering needed')
        return template_list

    logging.debug('%s too high similarities encountered for these templates: %s',
                  high_sim_count, ', '.join(list(set(suspicious_template_list))))

    # identify the wrong structure(s)

    if not high_sim_count:
        logging.debug('No templates have a similarity >90 and a difference >25, no template filtering needed')
        return template_list

    logging.debug('%s too high similarities encountered for these templates: %s',
                  high_sim_count, ', '.join(list(set(suspicious_template_list))))

    # identify the wrong structure(s)
    remove_template = None

    count_templates_dict = {}
    for template in set(suspicious_template_list):
        count = suspicious_template_list.count(template)
        if count not in count_templates_dict:
            count_templates_dict[count] = [template]
        else:
            count_templates_dict[count].append(template)

    max_count = sorted(count_templates_dict.keys(), reverse=True)[0]
    max_count_templates = count_templates_dict[max_count]

    if len(max_count_templates) == 1:
        # if multiple templates are closely related pick the one that is wrong
        # might be the one that occurs most in the suspicious_template_list
        remove_template = max_count_templates[0]
        logger.info('Remove %s', remove_template)
    else:
        # decide which template to remove
        logger.info('Remove one of these: %s', ', '.join(max_count_templates))

        # decide on
        # pub_date, resolution, seq_length
        solr = solr_talk.SolrPDB(global_settings)
        pdb_list = list(set([pdb[:-1] for pdb in max_count_templates]))
        pdb_info_dict = solr.get_pdb_information(pdb_list)

        # loop over the pdbs to find possible pdbs to remove
        oldest_template = ''
        oldest_template_year = None
        closest_year_diff = 100

        bad_resolution_list = []

        shortest_seq = ''
        shortest_seq_length = None
        closest_seq_diff = 1000

        obsolete_templates = []

        for template in max_count_templates:
            pdb = template[:-1]

            # look at release year
            if pdb_info_dict[pdb]['status'].find('obsolete') != -1:
                year = 1980
                obsolete_templates.append(template)
            else:
                # if the pdb info dict doesn't contain a year for
                # this pdb id set year to 1980
                year = pdb_info_dict[pdb].get('year', 1980)

            if not oldest_template_year or year < oldest_template_year:
                oldest_template = template
                if not oldest_template_year:
                    oldest_template_year = year
                else:
                    closest_year_diff = abs(oldest_template_year - year)
                    oldest_template_year = year

            if template != oldest_template:
                # not the oldest PDB
                year_diff = abs(oldest_template_year - year)
                if year_diff < closest_year_diff:
                    closest_year_diff = year_diff

            # look at resolution
            if 'resolution' not in pdb_info_dict[pdb] or pdb_info_dict[pdb]['resolution'] > 3:
                bad_resolution_list.append(template)

            # look at seq length
            core_aa_count = template_aa_count_dict[template]
            if not shortest_seq or core_aa_count < shortest_seq_length:
                shortest_seq = template
                if not shortest_seq_length:
                    shortest_seq_length = core_aa_count
                else:
                    closest_seq_diff = abs(shortest_seq_length - core_aa_count)
                    shortest_seq_length = core_aa_count

            if template != shortest_seq:
                # not the shortest PDB
                seq_diff = abs(shortest_seq_length - core_aa_count)
                if seq_diff < closest_seq_diff:
                    closest_seq_diff = seq_diff

        # decide on this data which template should be removed, check order:
        #   1. obsolete
        #   2. release year + resolution
        #   3. resolution
        #   4. closest sequence length difference vs shortest sequence length rate
        #   5. release year
        #   6. closest sequence length difference

        if len(obsolete_templates) == 1:
            remove_template = obsolete_templates[0]

        elif closest_year_diff >= 4 or (closest_year_diff >= 2 and oldest_template in bad_resolution_list):
            remove_template = oldest_template

        elif len(bad_resolution_list) == 1:
            remove_template = bad_resolution_list[0]

        elif (closest_seq_diff + 0.1) / (shortest_seq_length + 0.1) >= 0.30:
            # closest is at least 30% larger
            remove_template = shortest_seq

        elif closest_year_diff >= 1:
            remove_template = oldest_template
        elif closest_seq_diff >= 1:
            remove_template = shortest_seq

        if not remove_template:
            logger.warning('Could not decide which template to remove yet')
            return template_list
        logger.info('Removing %s', remove_template)

    # copy template list
    new_template_list = template_list[:]
    if remove_template:
        new_template_list.remove(remove_template)

        if len(new_template_list) > 1:
            # call recursively to check if other templates should be removed as well
            logger.info('Template removed, starting new iteration')
            new_template_list = filter_templates(new_template_list, global_settings, family_id)

    return new_template_list


def remove_useless_gaps(aligned_templates):
    """
    Remove positions on which all sequences have gaps (that might happen as a
    result of fixing cores and/or removing incorrect sequences
    """
    new_aligned_templates = deepcopy(aligned_templates)
    positions_to_remove = {}  # key: core index, value: list of residue indices
    # first find positions to remove (per core)
    for core_index, core_seq in enumerate(aligned_templates.values()[0]):
        for res_index in range(len(core_seq)):
            found_residue = False
            for seq_id, cores in aligned_templates.iteritems():
                if cores[core_index][res_index] != '-':
                    found_residue = True
                    break

            if not found_residue:
                # less than two sequences have residues on this position
                if core_index not in positions_to_remove:
                    positions_to_remove[core_index] = []
                positions_to_remove[core_index].append(res_index)

    # remove the positions
    for core_index in sorted(positions_to_remove, reverse=True):
        for res_index in sorted(positions_to_remove[core_index], reverse=True):
            for seq_id, seq in new_aligned_templates.iteritems():
                new_aligned_templates[seq_id][core_index] = \
                    seq[core_index][:res_index] + seq[core_index][res_index + 1:]

    # remove empty cores
    for seq_id, seq in new_aligned_templates.iteritems():
        new_aligned_templates[seq_id] = [core for core in seq if core]
    return new_aligned_templates


def remove_too_gappy_seqs(aligned_templates, min_res_cutoff):
    """
    Remove templates with number of aligned residues lower than 'min_res_cutoff'
    """
    remove_templates = set()
    for seq_id, sequence in aligned_templates.iteritems():
        seq_str = "".join(sequence)
        res_number = len(seq_str) - seq_str.count('-')
        if res_number < min_res_cutoff:
            remove_templates.add(seq_id)
    if remove_templates:
        aligned_templates = {seq_id: sequence for seq_id, sequence in aligned_templates.iteritems()
                             if seq_id not in remove_templates}
    return aligned_templates


def get_worst_template(seq_ids_list, sources, release_years, resolutions):
    """
    Find the worst template from seq_ids_list, checks:
        1. source: keeps the one that's has the best source according to the
        following order:
            other < genbank < trembl < swiss-prot < pdb
        (Thus pdb templates are the most preferred ones)
        If the best source isn't pdb, and there are multiple templates with
        the best source then choose a random one out of these

        2. release year: keeps the newest one, if there are multiple with the most
        recent release year than chooses a random one
    """
    # first check if there is any sequence template
    source_order = ['other', 'genbank', 'trembl', 'swiss-prot', 'pdb']
    worst_templates = []
    worst_src = source_order[-1]
    worst_src_index = len(source_order) - 1
    for seq_id in seq_ids_list:
        src = sources[seq_id]
        src_index = source_order.index(src)
        if src == worst_src:
            worst_templates.append(seq_id)
        elif src_index < worst_src_index:
            # encountered a new worst source
            worst_templates = [seq_id]
            worst_src = src
            worst_src_index = src_index
    if worst_src != 'pdb':
        # means all templates are from pdb
        worst_tmpl = random.choice(worst_templates)
        return worst_tmpl

    # if there are structures with no resolution (so not X-ray), reduce the seq id list to only these
    no_resolution = [seq_id for seq_id in seq_ids_list if seq_id not in resolutions]
    if no_resolution and len(no_resolution) == 1:
            return no_resolution[0]

    # all templates are from pdb, determine the worst one based on
    # release year (if a pdb doesn't have a release year remove this one
    no_release = [seq_id for seq_id in seq_ids_list if seq_id not in release_years]
    if no_release:
        # TODO: figure out what to do here, for now
        # if no release year just take a random one
        # MAYBE we should jsut decide they're both incorrect?
        worst_tmpl = random.choice(no_release)
    elif no_resolution:
        worst_tmpl = min(no_resolution, key=lambda x: release_years[x])
    else:
        # same issue here as above -> what to do when two pdbs have
        # the oldest release year (now it just chooses the first one in list)
        worst_tmpl = max(seq_ids_list, key=lambda x: resolutions[x])
        if resolutions.values().count(resolutions[worst_tmpl]) > 1:
            # multiple templates have the same highest resolution
            seq_ids_filtered = [
                seq_id for seq_id in seq_ids_list if
                isclose(resolutions[seq_id], resolutions[worst_tmpl])
            ]
            worst_tmpl = min(seq_ids_filtered, key=lambda x: release_years[x])
    return worst_tmpl


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


def fix_full_seq_errors(full_seq_errors, core_errors, sources, release_years, resolutions, aligned_templates):
    """
    Remove incorrectly aligned sequences
    """
    full_seq_errors_tmp = deepcopy(full_seq_errors)
    new_aligned_templates = deepcopy(aligned_templates)
    while full_seq_errors_tmp:
        ids_set = full_seq_errors_tmp.pop()
        tmpl_id_rm = get_worst_template(ids_set, sources, release_years, resolutions)
        if tmpl_id_rm in new_aligned_templates:
            del new_aligned_templates[tmpl_id_rm]
        # remove id from full seq errors
        for i in full_seq_errors:
            if tmpl_id_rm in i and i in full_seq_errors_tmp:
                full_seq_errors_tmp.remove(i)
        # remove id from core errors
        dict_tmp = deepcopy(core_errors)
        for i in dict_tmp:
            if tmpl_id_rm in i:
                del core_errors[i]

    return new_aligned_templates


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


def update_cutoffs(cutoffs):
    """
    Ask user if he wants to adjust the identity cutoff values;
    update them accordingly
    """
    answer = raw_input("Do you want to change identity cut-offs? (default: no) ")
    if answer.lower() not in ['y', 'yes']:
        # leave cutoffs unchanged
        return cutoffs

    while answer:
        answer = raw_input("Which identity cut-off do you want to change? ({})".format(
                ", ".join(cutoffs.keys())))
        if answer in cutoffs:
            new_value = float(raw_input("Give new value: "))
            if new_value > 1:
                print "New value needs to be <= 1"
            else:
                cutoffs[answer] = new_value
        else:
            print "'{}' is not a valid option".format(answer)
    return cutoffs


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

    ## write output files
    #if check_result['errors'] or check_result['full_seq_errors'] or write_log:

    #    if check_result['errors'] or check_result['full_seq_errors']:
    #        outpath = final_core_path + "_new"
    #        # write new core file
    #        write_core_file(templates_order, check_result['new_templates'], outpath)

    #    # write errors to error file
    #    write_core_errors_to_file(aligned_templates, check_result, corevar_path)

    ## finished
    #if manual:
    #    raw_input("done")


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
    if len(aligned_templates) > 70:
        tmpl_id = aligned_templates.keys()[0]
        aligned_templates_to_check = {tmpl_id: aligned_templates[tmpl_id]}
        checks_to_run = len(aligned_templates) - 1
    else:
        aligned_templates_to_check = aligned_templates
        checks_to_run = (len(aligned_templates) * (len(aligned_templates) - 1)) / 2

    # run per-template core checks
    for seq_id1, cores1 in aligned_templates_to_check.iteritems():
        full_seq1 = "".join(cores1).replace("0", "").upper()
        for seq_id2, cores2 in aligned_templates.iteritems():
            full_seq2 = "".join(cores2).replace("0", "").upper()

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
            full_identity = aliDef.calculateCoreIdentity(full_seq1, full_seq2)
            if 0.1 < full_identity < tmpl_identity_cutoff:
                # only compare highly identical templates
                continue
            elif full_identity <= 0.1:
                # identity suspiciously low, check full seq alignment
                full_seq_mafft_identity_cutoff = 0.45
                incorrect_sequences = check_cores(
                    [full_seq1], [full_seq2], full_seq_mafft_identity_cutoff, whatif_identity_cutoff,
                    seq_id1, seq_id2, diff_check=diff_check)
                if incorrect_sequences:
                    possible_full_seq_errors.add(ids_set)

            # run core-by-core check
            incorrect_cores = check_cores(cores1, cores2, mafft_identity_cutoff, whatif_identity_cutoff,
                                          seq_id1, seq_id2, diff_check=diff_check)
            if incorrect_cores:
                possible_core_errors[ids_set] = incorrect_cores

            # check if all pairs have been checked
            if len(checked_pairs) == all_pairs:
                break

    # write lgos to file
    if "coredir" in global_settings and write_log:
        logpath = "template_core_check.log"
        errors = {'errors': possible_core_errors, 'full_seq_errors': possible_full_seq_errors}
        write_core_errors_to_file(aligned_templates, errors, logpath)

    # if no errors were found return
    if not (possible_core_errors or possible_full_seq_errors):
        logger.info("No errors found")
        if manual:
            raw_input("done")
        return {
            'new_templates': aligned_templates,
            'old_templates': aligned_templates,
            'errors': possible_core_errors,
            'full_seq_errors': possible_full_seq_errors
        }

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
