import os

from paths import CSS, TEMPLATE


def aln_to_html(aa_aln, wrong):
    html_out = ""
    aln_length = len(aa_aln)
    for seq_id, seq in aa_aln.iteritems():
        html_sequence = "{}    ".format(seq_id)
        for r, res in enumerate(seq):
            if res != "-" and res != " ":
                if r in wrong[seq_id].keys():
                    level = get_level(wrong[seq_id][r], aln_length)
                    new_res = "<span class=featWRONG{}>{}</span>".format(level,
                                                                         res)
                else:
                    new_res = "<span class=featOK>{}</span>".format(res)
            else:
                new_res = res
            html_sequence += new_res
        html_out += html_sequence + "\n"
    return html_out


def get_level(number, aln_length):
    n = float(number) / aln_length
    if n >= 0.8:
        return 5
    elif n >= 0.6:
        return 4
    elif n >= 0.4:
        return 3
    elif n >= 0.2:
        return 2
    elif n >= 0:
        return 1


def aln_to_html_var(num_aln, wrong, full_seq):
    html_out = ""
    var_lengths = get_max_var_lengths(num_aln)
    for seq_id, seq in num_aln.iteritems():
        html_seq = make_html_seq_var(seq, full_seq[seq_id], var_lengths,
                                     wrong[seq_id])
        html_sequence = "{}    {}".format(seq_id, html_seq)
        html_out += html_sequence + "\n"
    return html_out


def make_html_seq_var():
    # TODO: implement
    pass


def get_max_var_lengths(num_aln):
    max_lengths = []
    for v in range(len(num_aln.values()[0])):
        lengths = [len(num_aln[s_id]["vars"][v]) for s_id in num_aln.keys()]
        max_lengths.append(max(lengths))
    return max_lengths


def write_html(aln, wrong_cols, outname):
    outtxt = aln_to_html(aln, wrong_cols)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    css_full_path = os.path.join(script_dir, CSS)
    with open(os.path.join(script_dir, TEMPLATE)) as a:
        template_fmt = a.read()
    with open(outname + ".html", 'w') as out:
        out.write(template_fmt.format(css_full_path, outtxt))
