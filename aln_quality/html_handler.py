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


def write_html(aln, wrong_cols, outname):
    outtxt = aln_to_html(aln, wrong_cols)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    css_full_path = os.path.join(script_dir, CSS)
    with open(os.path.join(script_dir, TEMPLATE)) as a:
        template_fmt = a.read()
    with open(outname + ".html", 'w') as out:
        out.write(template_fmt.format(css_full_path, outtxt))


def seq(num_aln, aa_aln, full_seq, max_lengths, wrong):
    aa_aln_corvar = make_corvar(full_seq, num_aln)
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
    return html_sequence


def make_corvar(full_seq, num_aln):
    corvar_aln = {seq_id: {"cores": [], "vars": []} for
                  seq_id in num_aln.keys()}
    for seq_id, seq in num_aln.iteritems():
        for c in seq["cores"]:
            new_core = ""
            for res in c:
                new_res = full_seq[seq_id][res - 1]
                new_core += new_res
            corvar_aln[seq_id]["cores"].append(new_core)
        for c in seq["vars"]:
            new_var = ""
            for res in c:
                new_res = full_seq[seq_id][res - 1]
                new_var += new_res
            corvar_aln[seq_id]["vars"].append(new_var)
    return corvar_aln
