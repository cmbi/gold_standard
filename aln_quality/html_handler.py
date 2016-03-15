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


def aln_to_html_var(num_aln, aa_aln, wrong, full_seq):
    html_out = ""
    aln_length = len(aa_aln)
    aa_aln_corvar = make_corvar(full_seq, num_aln)
    var_lengths = get_max_var_lengths(num_aln)
    for seq_id, seq in aa_aln_corvar.iteritems():
        html_seq = make_html_var_seq(
            seq, wrong[seq_id], var_lengths, aln_length)
        html_sequence = "{}    {}".format(seq_id, html_seq)
        html_out += html_sequence + "\n"
    return html_out


def get_max_var_lengths(num_aln):
    max_lengths = []
    for v in range(len(num_aln.values()[0])):
        lengths = [len(num_aln[s_id]["vars"][v]) for s_id in num_aln.keys()]
        max_lengths.append(max(lengths))
    return max_lengths


def write_html(aa_aln, wrong_cols, outname, var=False, num_aln={}, full_seq={}):
    if var:
        outtxt = aln_to_html(aa_aln, wrong_cols)
    else:
        outtxt = aln_to_html_var(num_aln, aa_aln, wrong_cols, full_seq)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    css_full_path = os.path.join(script_dir, CSS)
    with open(os.path.join(script_dir, TEMPLATE)) as a:
        template_fmt = a.read()
    with open(outname + ".html", 'w') as out:
        out.write(template_fmt.format(css_full_path, outtxt))


def make_html_var_seq(corvar_seq, wrong, max_lengths, aln_length):
    html_seq = ""
    r_index = 0
    for c, core in enumerate(corvar_seq["cores"]):
        var = corvar_seq["vars"][c]
        html_seq += var + " " * (max_lengths[c] - len(var))
        r_index += len(var)
        for r, res in enumerate(core):
            if res != "-" and res != " ":
                if r in wrong.keys():
                    level = get_level(wrong[r], aln_length)
                    new_res = "<span class=featWRONG{}>{}</span>".format(
                        level, res)
                else:
                    new_res = "<span class=featOK>{}</span>".format(res)
            else:
                new_res = res
            r_index += 1
            html_seq += new_res
    var = corvar_seq["vars"][-1]
    html_seq += var + " " * (max_lengths[-1] - len(var))
    return html_seq


def make_corvar(full_seq, num_aln):
    corvar_aln = {seq_id: {"cores": [], "vars": []} for
                  seq_id in num_aln.keys()}
    for seq_id, seq in num_aln.iteritems():
        for c in seq["cores"]:
            new_core = ""
            for res in c:
                if res == '-':
                    new_res = '-'
                else:
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
