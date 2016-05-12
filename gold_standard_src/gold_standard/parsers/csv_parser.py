import logging
import os
import re

from gold_standard_src.gold_standard.parsers.error_types import ParserError


_log = logging.getLogger(__name__)


def parse_csv_alignment(inpath, gold_ids):
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
    # iterate through the file ([1:] because first line is the header)
    in_aln = {}
    start = 1 if "alnsequence" in in_csv[0] else 0
    for line in in_csv[start:]:
        # split up the line (comma delimiter) and remove quotes
        # and join it again in one line
        # only take fields 1,2 cause field one is proteinid
        tmp_line = ', '.join([re.sub('"', '', x) for x in line.split(',')[1:]])
        new_seq = process_sequence(tmp_line)
        if new_seq.keys()[0] in gold_ids:
            in_aln.update(new_seq)
    return in_aln


def process_sequence(line):
    """
    Process sequence line from the 3DM alignment
    e.g.
        - input line: 1ABCA, ABC dfg FGH 0 IJK
        - output:
            {'1ABCA': {'cores': ['ABC', 'FGH', 'IJK'], 'vars': ['', 'dfg', '']
            }}
    :param line: sequence line
    :return: sequence dictionary
    """
    seq_id = line.split(',')[0]
    next_seg = 'var'
    if line.split()[1].isupper():
        # if the first segment is a core add an empty var segment
        next_seg = 'core'
    sequence = ""
    # start at 1, because the 0th segment is the ID
    for segment in line.split()[1:]:
        if next_seg == 'core':
            sequence += segment
            next_seg = 'var'
        else:
            next_seg = 'core'
    return {seq_id: sequence}
