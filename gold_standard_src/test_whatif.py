import argparse
import logging
import os
import sys

import calc_alignment_quality as ca
from gold_standard.html_handler import HtmlHandler

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def get_input_paths(project_name):
    project_path= os.path.join('families/from_scop/', project_name)
    if not os.path.exists(project_path):
        project_path= os.path.join('families/manual/', project_name)
        if not os.path.exists(project_path):
            print "No such project. Exiting"
            sys.exit()

    input_paths = {
        "gold_path": os.path.join(project_path, 'final_core.txt.Var'),
        "aln_path": os.path.join('WIFALI', project_name + '.3SSP'),
        "final_core": "",
        "gold_dir": ""
    }
    return input_paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("project_name")
    args = parser.parse_args()
    input_paths = get_input_paths(args.project_name)
    output = os.path.join('WIFALI', args.project_name + "_test_result")
    quality_data = ca.calculate_aln_quality(
        input_paths, output, 'fatcat', True)
    hh = HtmlHandler(var=False, var_short=False)
    hh.write_html(quality_data, output)
