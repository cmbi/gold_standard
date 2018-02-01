import logging
import os
from copy import deepcopy

from .error_types import ParserError
from .var_file import parse_var_file

_log = logging.getLogger(__name__)
fs = frozenset


def parse_gold_pairwise(gold_dir):
    _log.info("Getting gold standard alignments")
    if not os.path.exists(gold_dir):
        raise ParserError("No such directory: {}".format(gold_dir))
    # get all ".Var" files in the given directory
    var_list = [x for x in os.listdir(gold_dir) if x.endswith(".Var")]
    _log.info("Got %s var files", len(var_list))
    gold_alns = {}
    gold_ids = set()
    full_seq = {}
    for v in var_list:
        var = parse_var_file(os.path.join(gold_dir, v))
        gold_ids.update(var['ids'])
        aln_ids = fs(var['ids'])
        gold_alns[aln_ids] = var['alns']
        full_seq.update(var['full_seq'])
    _log.info("Finished parsing .Var files")
    return {
        'alns': gold_alns,
        'ids': gold_ids,
        'full_seq': full_seq
    }


def parse_gold_multi(gold_path):
    corvar = parse_var_file(gold_path, multi=True)
    corvar = fill_in_target(corvar)
    return {'alns': corvar['alns'], 'full_seq': corvar['full_seq'],
            'ids': corvar['alns']['cores'].keys()}


def fill_in_target(corvar):
    new_cores = corvar["alns"]["cores"]
    new_vars = corvar["alns"]["var"]
    print new_vars
    target_id = corvar["target"]

    target_vars = corvar["alns"]["var"][target_id]

    # check if all vars are empty
    if not any(target_vars):
        print "skipping"
        return corvar

    for i, var_i in enumerate(target_vars):
        if var_i:
            # var not empty add it to previous (or next if it's the first core)
            # and add gaps in all other sequences
            new_vars[target_id][i] = ""
            if i == 0:
                # new_cores[target_id][i] = var_i.upper() + new_cores[target_id][i]
                new_cores[target_id][i] = var_i + new_cores[target_id][i]
                core = "next"
            else:
                new_cores[target_id][i - 1] = new_cores[target_id][i - 1] + var_i
                core = "prev"

            for seq_id in corvar["alns"]["cores"]:
                if seq_id == target_id:
                    continue
                if core == "next":
                    new_cores[seq_id][i] = "-" * len(var_i) + new_cores[seq_id][i]
                else:
                    new_cores[seq_id][i - 1] = new_cores[seq_id][i - 1] + "-" * len(var_i)
    new_corvar = deepcopy(corvar)
    new_corvar["alns"]["cores"] = new_cores
    new_corvar["alns"]["var"] = new_vars
    return new_corvar
