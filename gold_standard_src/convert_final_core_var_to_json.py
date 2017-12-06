import argparse
import json
import os

from gold_standard.parsers import var_file as vf


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
        score_modifiers: {'a': 1, 'u': -1, , 'd': 0, 'm0': 0, 'm1': 0.1, 'm2': 0.2}
        }
        alignments: {
            strct_id_1: {

            }
        }

    """


def run_convert(inpath, outpath):
    """
    Convert final_core.txt.Var to a json alignment file
    """
    # parse input
    corvar_data = vf.parse_var_file(inpath, multi=True)

    # convert to dict
    json_alignment = convert_corvar_to_json_alignment(corvar_data)

    # write to json file
    with open(outpath, "w") as o:
        json.dump(json_alignment, o)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert final_core.txt.Var to a json alignment file")
    parser.add_argument("inpath")
    parser.add_argument("outpath")

    args = parser.parse_args()
    run_convert(args.inpath, args.outpath)
