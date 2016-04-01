import copy
import os

from src.gold_standard.paths import CSS, TEMPLATE


class HtmlHandler(object):
    def __init__(self, var=False, var_short=False, long_len=20, short=False):
        self.var = var
        self.var_short = var_short
        self.short = short
        self.long_len = long_len

    def write_html(self, aa_aln, wrong_cols, outname, num_aln=None,
                   full_seq=None, core_indexes=None):
        if self.var or self.var_short:
            outtxt = self.aln_to_html_var(num_aln, aa_aln, wrong_cols, full_seq,
                                          core_indexes)
        else:
            outtxt = self.aln_to_html(aa_aln, wrong_cols)
        script_dir = os.path.dirname(os.path.realpath(__file__))
        css_full_path = '/'.join(list(os.path.split(script_dir)[:-2]) + [CSS])
        tmpl_full_path = '/'.join(list(os.path.split(script_dir)[:-2]) +
                                  [TEMPLATE])
        with open(tmpl_full_path) as a:
            template_fmt = a.read()
        with open(outname + ".html", 'w') as out:
            out.write(template_fmt.format(css_full_path, outtxt))

    def aln_to_html_var(self, num_aln, aa_aln, wrong, full_seq, core_indexes):
        html_out = ""
        aln_length = len(aa_aln)

        num_aln_c = self.split_cores(num_aln, core_indexes)
        num_aln_v = self.split_vars(num_aln_c)

        aa_aln_corvar = self.make_corvar(full_seq, num_aln_v)
        var_lengths = self.get_max_var_lengths(num_aln_v)
        for seq_id, seq in aa_aln_corvar.iteritems():
            html_seq = self.make_html_var_seq(
                seq, wrong[seq_id], var_lengths, aln_length)
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

    def aln_to_html(self, aa_aln, wrong):
        html_out = ""
        aln_length = len(aa_aln)
        for seq_id, seq in aa_aln.iteritems():
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
