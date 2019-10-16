import logging
import os
import re
import tempfile

from .error_types import ParserError

_log = logging.getLogger(__name__)


def parse_fasta(aln_path, golden_ids=None):
    with open(aln_path) as a:
        aln_lines = a.read().splitlines()

    aln_dict = {}
    seq_id = ""
    found_fasta_headers = False
    for l in aln_lines:
        if l.startswith(">"):
            found_fasta_headers = True
            if '|' in l:
                seq_id = l.lstrip('>').split('|')[0]
            else:
                # remove the fasta header symbol
                seq_id = l.lstrip('>')
                # remove whitespace characters
                seq_id = re.sub(r'\s', '', seq_id)
            if seq_id in aln_dict.keys():
                raise ParserError("Sequence {} is duplicated".format(seq_id))
            if golden_ids is None or seq_id in golden_ids:
                aln_dict[seq_id] = ""
        elif seq_id and (golden_ids is None or seq_id in golden_ids):
            aln_dict[seq_id] += l.upper().rstrip("*")
    if not found_fasta_headers:
        raise ParserError("Incorrect format, no fasta headers found.")

    return aln_dict


def write_fasta(sequences, name=None, prefix=''):
    """
    Write sequences to a temporary fasta file
    :param sequences: sequences dict, key: id, value: seq string
    :return: tmp path
    """
    if name:
        outfile = open(name, 'w')
    else:
        outfile = tempfile.NamedTemporaryFile(prefix=prefix, suffix=".fasta", delete=False)

    _log.debug("Writing %s sequences to a FASTA file: %s", len(sequences), outfile.name)
    # TODO: comment out the two lines below after debugging is finished
    k = sequences.keys()[0]
    _log.debug("Sequence %s: %s", k, sequences[k])
    fasta_sequences = '\n'.join(
        ['>{}\n{}'.format(s_id, s) for s_id, s in sequences.iteritems()]
    ) + "\n"
    with outfile as f:
        f.write(fasta_sequences)
    return outfile.name
