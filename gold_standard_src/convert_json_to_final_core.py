import argparse
import json

from gold_standard.parsers import var_file as vf


def convert_pairwise_aln_to_json(target_cores, core_aln_strct2, var_regions_strct2):
    """
    converts a pairwise alignment to json
    :param target_cores: grounded target core sequence
    :param core_aln_strct2: grounded sequence of the 2nd structure
    :param var_regions_strct2: list of residues in the variable regions (2nd structure)
    :return: json alignment
    """
    all_strct2_residues = filter(lambda x: x != '-', sorted(core_aln_strct2 + var_regions_strct2))

    pairwise_json_aln = {}
    for res_num in all_strct2_residues:
        if res_num == '-':
            continue
        if res_num in core_aln_strct2:
            index = core_aln_strct2.index(res_num)
            pairwise_json_aln[res_num] = {target_cores[index]: 'a'}
        else:
            pairwise_json_aln[res_num] = {'*': 'u'}
    return pairwise_json_aln


def convert_corvar_to_json_alignment(corvar_data):
    """
    Convert corvar data to a json alignment
    :param corvar_data: dictionary: {
        ids: list of template ids
        alns: {
            cores: dictionary of grounded aligned cores
            vars: dictionary of grounded aligned cores
        }
        full: dictionary of full sequences
        target: id of the target structure
    }

    :return: dictionary of alignment all-vs-target:
        {
        target: target id,
        # dictionary of modifiers for each alignment category
        # for now the categories are:
        #   'a': aligned
        #   'u': unaligned
        #   'd' don't know
        #   'm1', 'm2', ..., 'mN' - multiple solutions
        score_modifiers: {'a': 1, 'u': -1, 'd': 0, 'm0': 0, 'm1': 0.1, 'm2': 0.2}
        }
        alignments: {
            strct_id_1: {
                # residue number 1 from strct 1
                #  should be aligned with residue 2 from the target seq
                1: {2: 'a'},

                # residue number 2 from strct 1
                #  can be aligned with residue 3 (for score m1) from the target seq
                #  or residue 4 (for score m2)
                2: {3: 'm1', 4: 'm2'},

                # give the "unaligned" penalty score if residue 3 is aligned with anything
                3: {'*': 'u'},

                # give the "don't know" penalty score if residue 4 is aligned with anything
                4: {'*': 'd'},
            }
        }
    """
    score_modifiers = {'a': 1, 'b': 0.8, 'c': 0.5, 'd': 0, 'u': -1}

    var_regions = corvar_data["alns"]["var"]

    # retrieve target alignment
    target_id = corvar_data["target"]
    target_cores = corvar_data["alns"]["cores"][target_id]

    # convert alignments
    json_alignments = {}
    for strct_id, core_aln in corvar_data["alns"]["cores"].iteritems():
        # if strct_id == target_id:
        #     # don't compare the target to itself
        #     continue

        json_pairwise_alignment = convert_pairwise_aln_to_json(target_cores, core_aln, var_regions[strct_id])

        json_alignments[strct_id] = json_pairwise_alignment

    result = {
        "alignments": json_alignments, "target": corvar_data["target"], "score_modifiers": score_modifiers
    }
    return result


def get_best_solution(solutions):
    if len(solutions) == 1:
        return solutions.keys()[0]
    else:
        best_score = 0
        best_res_number = None

        for res_number, score in solutions.iteritems():
            score = float(score.lstrip("m")) / 10
            if score > best_score:
                best_score = score
                best_res_number = res_number
    return best_res_number


def convert_json_to_final_core(json_aln, full_sequences, master_id):
    inverted_json_aln = {}

    # invert the aln dict: {master position: {seq id: test seq position}}
    prev_positions = {seq_id: 0 for seq_id in full_sequences}  # last encountered position from each strct
    lowercase = {seq_id: [] for seq_id in full_sequences}
    for seq_id, aligned_seq in json_aln.iteritems():
        for res_from, solutions in aligned_seq.iteritems():
            if "*" in solutions:
                continue

            res_to = get_best_solution(solutions)
            if res_to not in inverted_json_aln:
                inverted_json_aln[res_to] = {}
            inverted_json_aln[res_to][seq_id] = res_from

            prev_pos = prev_positions[seq_id]
            if int(prev_pos) != int(res_from) - 1:
                lowercase[seq_id].extend([prev_pos, res_from])
            prev_positions[seq_id] = res_from

    aln_dict = {
        master_id: full_sequences[master_id]
    }

    for i in range(1, len(full_sequences[master_id])):
        for seq_id in full_sequences.keys():
            if seq_id == master_id:
                continue

            curr_pos_aln = inverted_json_aln[str(i + 1)]
            if seq_id not in curr_pos_aln:
                residue = "-"
            else:
                res_num = curr_pos_aln[seq_id]
                residue = full_sequences[seq_id][int(res_num) - 1]
                if res_num in lowercase[seq_id]:
                    print "make lower"
                    residue = residue.lower()

            if seq_id not in aln_dict:
                aln_dict[seq_id] = ""
            aln_dict[seq_id] += residue
    return aln_dict


def run_convert(inpath, incorvar, outpath):
    """
    Convert final_core.txt.Var to a json alignment file
    """
    # parse input
    corvar_data = vf.parse_var_file(incorvar, multi=True)
    full_sequences = corvar_data['full_seq']
    master_id = corvar_data['target']

    with open(inpath) as a:
        json_aln = json.load(a)["alignments"]

    aln_dict = convert_json_to_final_core(json_aln, full_sequences, master_id)
    with open(outpath, "w") as o:
        o.write(
                "\n".join(
                        ["{}  {}".format(seq_id, sequence) for seq_id, sequence in aln_dict.iteritems()]
                ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert final_core.json to a final_core alignment file")
    parser.add_argument("-i", "--injson")
    parser.add_argument("-c", "--incorvar")
    parser.add_argument("-o", "--outpath")

    args = parser.parse_args()
    run_convert(args.injson, args.incorvar, args.outpath)
