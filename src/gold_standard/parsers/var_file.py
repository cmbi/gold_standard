import logging
import re

from src.gold_standard.num_seq import aln_seq_to_num
from src.gold_standard.parsers.error_types import ParserError

_log = logging.getLogger(__name__)


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
                raise ParserError("Core {} and {} have different "
                                  "lengths".format(var1[i], var2[j]))
            aln_seq2 += var2[j]
            aln_seq1 += var1[i]
            i += 1
            j += 1
        else:
            raise ParserError("Incorrect Var file: check number of cores")
    return {id1: aln_seq1, id2: aln_seq2}


def parse_var_file(file_path):
    _log.debug("Parsing var file: %s", file_path)
    with open(file_path) as a:
        var_file = a.read().splitlines()
    ids = [i.split(',')[0] for i in var_file]
    aln = convert_var_to_aln(var_file)
    full = {}
    full[ids[0]] = re.sub('-', '', aln[ids[0]])
    full[ids[1]] = re.sub('-', '', aln[ids[1]])
    num_aln = {seq_id: aln_seq_to_num(seq) for seq_id, seq in aln.iteritems()}
    return {'ids': ids, 'aln': num_aln, 'full': full}