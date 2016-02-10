def aln_to_html(quality_data):
    aa_aln = quality_data["aln"]
    wrong = quality_data["wrong_cols"]
    html_out = ""
    for seq_id, seq in aa_aln.iteritems():
        html_sequence = ""
        for r, res in enumerate(seq):
            if res != "-" and res != " ":
                if r in wrong:
                    new_res = "<span class=featWRONG>{}</span>".format(res)
                else:
                    new_res = "<span class=feat{}>{}</span>".format(res.upper(),
                                                                    res)
            else:
                new_res = res
            html_sequence += new_res
        html_out += html_sequence + "\n"
    return html_out
