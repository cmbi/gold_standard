def aln_to_html(quality_data):
    aa_aln = quality_data["aln"]
    wrong = quality_data["wrong_cols"]
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
                    # new_res = "<span class=feat{}>{}</span>".format(res.upper(),
                    #                                                 res)
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
