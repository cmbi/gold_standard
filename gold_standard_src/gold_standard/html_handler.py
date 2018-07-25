import copy
import logging
import os

from .paths import TEMPLATE


_log = logging.getLogger("__main__")


class HtmlHandler(object):
    def __init__(self, var=False, var_short=False, long_len=20, short=False, pairwise=False):
        self.var = var
        self.var_short = var_short
        self.short = short
        self.long_len = long_len
        self.pairwise = pairwise

    def write_html(self, quality_data, outname, complex_scoring=False):
        if self.var or self.var_short:
            outtxt = self.aln_to_html_var(quality_data)
        elif self.pairwise:
            outtxt = self.aln_to_html_pairwise(
                    quality_data['aa_aln'], quality_data["gold_aln"], quality_data["full"],
                    quality_data['wrong_cols'], quality_data["order"])
        elif not complex_scoring:
            outtxt = self.aln_to_html(
                quality_data['aa_aln'], quality_data['wrong_cols'], quality_data["order"])
        else:
            outtxt = self.complex_aln_to_html(
                    quality_data['aa_aln'], quality_data['wrong_cols'], quality_data["order"])

        script_dir = os.path.dirname(os.path.abspath(__file__))
        tmpl_full_path = '/'.join(list(os.path.split(script_dir)[:-1]) +
                                  [TEMPLATE])
        with open(tmpl_full_path) as a:
            template_fmt = a.read()

        css = """
        <style>
        monospacediv {
            font-family: monospace;
        }
        .featWRONG{
            font-family: monospace;
            background: #FF0000;
        }
        .featWRONG5{
            font-family: monospace;
            background: #FF2200;
        }
        .featWRONG4{
            font-family: monospace;
            background: #FF9900;
        }
        .featWRONG3{
            font-family: monospace;
            background: #FFAA00;
        }
        .featWRONG2{
            font-family: monospace;
            background: #FFCC00;
        }
        .featWRONG1{
            font-family: monospace;
            background: #FFFF00;
        }
        .featOK{
            font-family: monospace;
            background: #AAFF00;
        }
        .asteriskBlank{
            opacity: 1.0;
            font-family: monospace;
        }
        .asteriskFull{
            font-family: monospace;
        }
        </style>
        """
        with open(outname + ".html", 'w') as out:
            out.write(template_fmt.format(css, outtxt))

    def aln_to_html_var(self, quality_data):
        # html_out = ""
        html_out = "<div class=monospacediv>\n<br>"
        aln_length = len(quality_data['aa_aln'])
        num_aln_c = self.split_cores(quality_data['num_aln'],
                                     quality_data['core_indexes'])
        num_aln_v = self.split_vars(num_aln_c)

        aa_aln_corvar = self.make_corvar(quality_data['full'], num_aln_v)
        var_lengths = self.get_max_var_lengths(num_aln_v)
        for seq_id in quality_data["order"]:
            seq = aa_aln_corvar[seq_id]
            html_seq = self.make_html_var_seq(
                seq, quality_data['wrong_cols'][seq_id], var_lengths,
                aln_length)
            html_sequence = "{}    {}".format(seq_id, html_seq)
            html_out += html_sequence + "\n"
        return html_out

    def get_max_var_lengths(self, num_aln):
        max_lengths = []
        for v in range(len(num_aln['var'].values()[0])):
            lengths = [len(num_aln['var'][s_id][v]) for
                       s_id in num_aln['var'].keys()]
            m_len = max(lengths)
            if self.short and m_len > self.long_len:
                # if max var is longer than allowed
                m_len = self.long_len + len(str(m_len - self.long_len)) + 2

            max_lengths.append(m_len)
        return max_lengths

    def make_html_var_seq(self, corvar_seq, wrong, max_lengths, aln_length):
        html_seq = ""
        r_index = 0
        for c, core in enumerate(corvar_seq["cores"]):
            var = corvar_seq["var"][c].lower()
            if self.short:
                var = self.make_short_var(var)
            var = " " + var + " " * (max_lengths[c] - len(var)) + " "
            html_seq += var
            for res in core:
                if res != "-" and res != " ":
                    if r_index in wrong.keys():
                        level = self.get_level(wrong[r_index], aln_length)
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                            level, res)
                    else:
                        new_res = "<span class=featOK>{}</span>".format(res)
                else:
                    new_res = res
                r_index += 1
                html_seq += new_res
        var = corvar_seq["var"][-1].lower()
        if self.short:
            var = self.make_short_var(var)
        html_seq += " " + var + " " * (max_lengths[-1] - len(var))
        return html_seq

    def make_short_var(self, var):
        if len(var) > self.long_len:
            edge = self.long_len / 2
            var_start = var[:edge]
            var_end = var[-edge:]
            short_var = var_start + " " + str(len(var) - self.long_len) + \
                " " + var_end
        else:
            short_var = var
        return short_var

    def aln_to_html_pairwise(self, aa_aln, gold_aln, full, wrong, order):
        html_out = "<div class=monospacediv>\n<br>"
        aln_length = len(aa_aln)

        _log.info("Creating pairiwse html")

        for seq_id in order:
            seq = aa_aln[seq_id]
            html_sequence = "<b>TEST</b> {}    ".format(seq_id)
            html_gold_sequence = "<b>GOLD</b> {}    ".format(seq_id)
            asterisk_line = "<span class=asteriskBlank>         {}</span>".format(" " * len(seq_id))
            gold_seq = gold_aln["cores"][seq_id]
            if len(gold_seq) != len(seq):
                msg = "Sequences not of equal length ({} vs {}):\n{}\n{}".format(len(gold_seq), len(seq), gold_seq, seq)
                _log.error(msg)
                raise RuntimeError(msg)
            for r, res in enumerate(seq):
                if gold_seq[r] != "-":
                    gold_aa_index = gold_seq[r] - 1
                    gold_aa = full[seq_id][gold_aa_index]
                else:
                    gold_aa = "-"
                if res != "-" and res != " ":
                    if r in wrong[seq_id].keys():
                        level = self.get_level(wrong[seq_id][r], aln_length)
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                                level, res)
                        new_gold_res = "<span class=featWRONG{}>{}</span>".format(
                                level, gold_aa)

                    else:
                        new_res = "<span class=featOK>{}</span>".format(res)
                        new_gold_res = "<span class=featOK>{}</span>".format(gold_aa)
                else:
                    new_res = res
                    new_gold_res = gold_aa

                if new_res.upper() != new_gold_res.upper() and new_res != "-":
                    asterisk_line += "<span class=asteriskFull>*</span>"
                else:
                    asterisk_line += "<span class=asteriskBlank> </span>"

                html_sequence += new_res
                html_gold_sequence += new_gold_res
            html_out += asterisk_line + "\n"
            html_out += html_sequence + "\n"
            html_out += html_gold_sequence + "\n"
            html_out += "<br>"
        html_out += "</div>"
        _log.info("Finished creating pairwise html")
        return html_out

    def complex_aln_to_html(self, aa_aln, wrong, order):
        # html_out = ""
        html_out = "<div class=monospacediv>\n<br>"
        aln_length = len(aa_aln)

        for seq_id in order:
            seq = aa_aln[seq_id]
            html_sequence = "{}    ".format(seq_id)
            for r, res in enumerate(seq):
                if res != "-" and res != " ":
                    score = wrong[seq_id][r + 1]
                    if score[0] and score[1] == 1:
                        if seq_id == "3M21A":
                            print score
                            raise
                        new_res = "<span class=featOK>{}</span>".format(res)
                    else:
                        level = self.get_level_cmplx(score[1])
                        print "level", level, "from", score[1]
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                            level, res)
                else:
                    new_res = res
                html_sequence += new_res
            html_out += html_sequence + "\n"
        return html_out

    def aln_to_html(self, aa_aln, wrong, order):
        # html_out = ""
        html_out = "<div class=monospacediv>\n<br>"
        aln_length = len(aa_aln)
        for seq_id in order:
            seq = aa_aln[seq_id]
            html_sequence = "{}    ".format(seq_id)
            for r, res in enumerate(seq):
                if res != "-" and res != " ":
                    if r in wrong[seq_id].keys():
                        level = self.get_level(wrong[seq_id][r], aln_length)
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                            level, res)
                    else:
                        new_res = "<span class=featOK>{}</span>".format(res)
                else:
                    new_res = res
                html_sequence += new_res
            html_out += html_sequence + "\n"
        return html_out

    @staticmethod
    def get_level_cmplx(number):
        n = 1 - float(number)
        if n >= 0.8:
            return 5
        elif n >= 0.6:
            return 4
        elif n >= 0.4:
            return 3
        elif n >= 0.2:
            return 2
        elif n > 0:
            return 1
        else:
            return 0

    @staticmethod
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

    @staticmethod
    def make_corvar(full_seq, num_aln):
        corvar_aln = {seq_id: {'cores': [], 'var': []} for
                      seq_id in num_aln['cores'].keys()}
        for seq_id, cores in num_aln['cores'].iteritems():
            for c in cores:
                new_core = ""
                for res in c:
                    if res == '-':
                        new_res = '-'
                    else:
                        new_res = full_seq[seq_id][res - 1]
                    new_core += new_res
                corvar_aln[seq_id]["cores"].append(new_core)
            for v in num_aln['var'][seq_id]:
                new_var = ""
                for res in v:
                    if res != '-':
                        new_res = full_seq[seq_id][res - 1]
                        new_var += new_res
                corvar_aln[seq_id]["var"].append(new_var)
        return corvar_aln

    @staticmethod
    def split_cores(num_aln, core_indexes):
        new_num_aln = copy.deepcopy(num_aln)
        for seq_id in num_aln['cores'].keys():
            new_num_aln['cores'][seq_id] = []
            for i, c_i in enumerate(core_indexes):
                if i < len(core_indexes) - 1:
                    next_c = core_indexes[i + 1]
                else:
                    # add last core
                    next_c = len(num_aln['cores'][seq_id])
                new_core = num_aln['cores'][seq_id][c_i:next_c]
                new_num_aln['cores'][seq_id].append(new_core)
        return new_num_aln

    def split_vars(self, num_aln):
        new_num_aln = copy.deepcopy(num_aln)
        for seq_id, cores in num_aln['cores'].iteritems():
            new_num_aln['var'][seq_id] = []
            prev = 1
            for core in cores:
                new = self.get_first_res(core)
                if new == -1:
                    new_num_aln['var'][seq_id].append([])
                else:
                    new_var = range(prev, new)
                    new_num_aln['var'][seq_id].append(new_var)
                    prev = self.get_first_res(core[::-1]) + 1
            if prev in num_aln['var'][seq_id]:
                p_index = num_aln['var'][seq_id].index(prev)
                new_num_aln['var'][seq_id].append(
                    num_aln['var'][seq_id][p_index:])
            else:
                new_num_aln['var'][seq_id].append([])

        return new_num_aln

    @staticmethod
    def get_first_res(core):
        for r in core:
            if isinstance(r, int):
                return r
        return -1
