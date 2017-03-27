import logging
import os
import re

from gold_standard_src.gold_standard.parsers.error_types import ParserError


_log = logging.getLogger(__name__)


def parse_csv_alignment(inpath, gold_ids, with_identities=False):
    """
    Parses CSV file containing aligned seuqences, first line is the header,
    second line is the template sequence:
    All lines contain sequence and sequence id:
        "sequence_id", "alnsequence"
    """
    _log.info("Parsing the input alignment [csv]")
    if not os.path.exists(inpath):
        raise ParserError("File not found: {}".format(inpath))
    # read the file
    with open(inpath) as a:
        in_csv = a.read().splitlines()
        num_aln = {'cores': {}, 'var': {}}
    core_indexes = []
    core_indexes_i = []
    aa_aln = {}
    for line in in_csv:
        if "alnsequence" in line:
            continue
        # split up the line (comma delimiter) and remove quotes
        # and join it again in one line
        # only take fields 1,2 cause field one is proteinid
        tmp_line = ','.join([re.sub('"', '', x) for x in line.split(',')])
        tmp_line_split = tmp_line.split(',')
        if len(tmp_line_split) == 4:
            corvar_line = tmp_line_split[3]
            seq_id = ''.join(tmp_line_split[1:3]).upper()
        elif len(tmp_line_split) == 2:
            corvar_line = tmp_line_split[1]
            seq_id = tmp_line_split[0].upper()
        new_seq_aa, new_seq_num, core_indexes_i = csv_corvar_to_num(corvar_line)
        # new_seq, core_indexes_i = process_sequence(tmp_line)
        if not core_indexes:
            core_indexes = core_indexes_i
        # assert core_indexes == core_indexes_i
        if seq_id in gold_ids:
            # if new_seq.keys()[0] in gold_ids:
            # num_aln.update(new_seq)
            num_aln['cores'][seq_id] = new_seq_num['cores']
            num_aln['var'][seq_id] = new_seq_num['var']
            aa_aln[seq_id] = new_seq_aa
    return aa_aln, num_aln, core_indexes


def csv_corvar_to_num(corvar_line):
    aln = {'cores': [], 'var': [], 'full': ''}
    # remove numbers and whitespaces form the corvar line
    sequence = re.sub(r'[0-9\s]', '', corvar_line)
    seq_split = corvar_line.split()
    aln['full'] = re.sub('-', '', sequence).upper()
    count = 1
    ex = False
    core_indexes = set()
    # total length of cores up to now
    # cores_len = 0
    prev = 'var'
    aa_seq = []
    for j, seg_j in enumerate(seq_split):
        for res_i in seg_j:
            if j % 2 == 0:
                # if res_i.isupper():
                #     ex = True
                if res_i != '-' and res_i != '0':
                    # even segments are vars
                    aln['var'].append(count)
                    count += 1
                prev = 'var'
            else:
                if prev == 'var':
                    core_indexes.add(len(aln['cores']))
                # if res_i.islower():
                #     ex = True
                if res_i == '-':
                    aln['cores'].append('-')
                else:
                    aln['cores'].append(count)
                    count += 1
                aa_seq.append(res_i)
                prev = 'core'
            if ex:
                raise ParserError("Core and var regions are not in the "
                                  "right order")
    return aa_seq, aln, core_indexes
