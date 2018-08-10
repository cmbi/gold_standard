import copy
import logging
import os

from .paths import TEMPLATE


_log = logging.getLogger("__main__")


class HtmlHandler(object):
    def __init__(self, long_len=20):
        self.long_len = long_len

    def write_html(self, quality_data, outname, mode="cores"):
        if mode in ["var", "var_short"]:
            outtxt = self.aln_to_html_var(quality_data, mode)

        elif mode == "pairwise":
            outtxt = self.aln_to_html_pairwise(
                    quality_data['aa_aln'], quality_data["gold_aln"], quality_data["full"],
                    quality_data['wrong_cols'], quality_data["order"])

        elif mode == "pairwise_complex":
            outtxt = self.aln_to_html_pairwise_complex(quality_data)

        elif mode == "cores":
            outtxt = self.aln_to_html(
                quality_data['aa_aln'], quality_data['wrong_cols'], quality_data["order"])

        elif mode == "cores_complex":
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
        .featMeh{
            font-family: monospace;
            background: #FFEF96;
        }
        .asteriskBlank{
            opacity: 1.0;
            font-family: monospace;
        }
        .asteriskFull{
            font-family: monospace;
        }
        .noFeat{
            font-family: monospace;
        }
        </style>
        """
        with open(outname + ".html", 'w') as out:
            out.write(template_fmt.format(css, outtxt))

    def aln_to_html_var(self, quality_data, mode):
        short_var = mode == "var_short"
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"
        aln_length = len(quality_data['aa_aln'])
        num_aln_c = self.split_cores(quality_data['num_aln'],
                                     quality_data['core_indexes'])
        num_aln_v = self.split_vars(num_aln_c)

        aa_aln_corvar = self.make_corvar(quality_data['full'], num_aln_v)
        var_lengths = self.get_max_var_lengths(num_aln_v, short_var)
        for seq_id in quality_data["order"]:
            if seq_id not in quality_data["gold_aln"]["cores"]:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue
            seq = aa_aln_corvar[seq_id]
            html_seq = self.make_html_var_seq(
                seq, quality_data['wrong_cols'][seq_id], var_lengths,
                aln_length, short_var)
            html_sequence = "{}    {}".format(seq_id, html_seq)
            html_out += html_sequence + "\n"
        return html_out

    def get_max_var_lengths(self, num_aln, short_var):
        max_lengths = []
        for v in range(len(num_aln['var'].values()[0])):
            lengths = [len(num_aln['var'][s_id][v]) for
                       s_id in num_aln['var'].keys()]
            m_len = max(lengths)
            if short_var and m_len > self.long_len:
                # if max var is longer than allowed
                m_len = self.long_len + len(str(m_len - self.long_len)) + 2

            max_lengths.append(m_len)
        return max_lengths

    def make_html_var_seq(self, corvar_seq, wrong, max_lengths, aln_length, short_var=False):
        html_seq = ""
        r_index = 0
        for c, core in enumerate(corvar_seq["cores"]):
            var = corvar_seq["var"][c].lower()
            if short_var:
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
                    new_res = "<span class=noFeat>" + res + "</span>"
                r_index += 1
                html_seq += new_res
        var = corvar_seq["var"][-1].lower()
        if short_var:
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

    def aln_to_html_pairwise_complex(self, quality_data):
        aa_aln = quality_data["aa_aln"]
        gold_aln = quality_data["gold_aln"]
        full = quality_data["full"]
        wrong = quality_data["wrong_cols"]
        order = quality_data["order"]
        target_id = quality_data["target_id"]

        target_seq = full[target_id]

        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"
        aln_length = len(aa_aln)

        _log.info("Creating pairiwse html")
        for seq_id in order:
            if seq_id not in gold_aln:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue
            for pos, res in enumerate(target_seq, start=1):
                # find out which residue is aligned with this in the gold aln
                pass

                # find out which residue is aligned with this in the test aln

        for seq_id in order:
            if seq_id not in gold_aln:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue

            seq = aa_aln[seq_id]
            html_sequence = "<b>TEST</b> {}    ".format(seq_id)
            html_gold_sequence = "<b>GOLD</b> {}    ".format(seq_id)
            asterisk_line = "<span class=asteriskBlank>         {}</span>".format(" " * len(seq_id))

            pairwise_gold_aln = gold_aln[seq_id]
            print seq_id, "running..."
            for position, score in wrong[seq_id].iteritems():
                print seq_id, position, len(full[seq_id])
                res = full[seq_id][position - 1]
                # if res != "-":
                #     real_pos += 1
                # else:
                #     continue

                # full_seq_index = self.get_index_in_full_seq()
                gold_aa = self.get_gold_aa(pairwise_gold_aln, full[seq_id], position)
                add_asterisk = False
                if res != "-" and res != " ":
                    if score[0] and score[1] == 1:
                        new_res = "<span class=featOK>{}</span>".format(res)
                        new_gold_res = "<span class=featOK>{}</span>".format(gold_aa)
                    else:
                        level = self.get_level_cmplx(score[1])
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                                level, res)
                        new_gold_res = "<span class=featWRONG{}>{}</span>".format(
                                level, gold_aa)
                        add_asterisk = True

                else:
                    if res == "-" and gold_aa == "-":
                        new_res = "<span class=noFeat>" + res + "</span>"
                        new_gold_res = "<span class=noFeat>" + gold_aa + "</span>"
                    elif res == "-":
                        new_res = "<span class=featMeh>" + res + "</span>"
                        new_gold_res = "<span class=featMeh>" + gold_aa + "</span>"
                        add_asterisk = True
                    else:
                        raise RuntimeError("res: {}; gold res: {}".format(res, gold_aa))
                if add_asterisk:
                    asterisk_line += "<span class=asteriskFull>*</span>"
                else:
                    asterisk_line += "<span class=asteriskBlank> </span>"
                # if res.upper() != gold_aa.upper() and new_res != "-":
                #     asterisk_line += "<span class=asteriskFull>*</span>"
                # else:
                #     asterisk_line += "<span class=asteriskBlank> </span>"

                html_sequence += new_res
                html_gold_sequence += new_gold_res
            html_out += asterisk_line + "\n"
            html_out += html_sequence + "\n"
            html_out += html_gold_sequence + "\n"
            html_out += "<br>"
            print seq_id, "OK"
        html_out += "</div>"
        _log.info("Finished creating pairwise html")
        return html_out

    @ staticmethod
    def get_gold_aa(pairwise_gold_aln, full_seq, pos):
        """
        Get residue in the gold aln that's aligned to residue 'pos'
        """
        try:
            aln_solutions = pairwise_gold_aln[str(pos)]
        except:
            for k in sorted(pairwise_gold_aln.keys()):
                print "R", k, pairwise_gold_aln[k]

            raise
        if "*" in aln_solutions:
            # this residue should not be aligned to anything
            gold_aa = "-"
        elif len(aln_solutions) == 1:
            # there is only one solution, take it
            aa_index = int(aln_solutions.keys()[0]) - 1
            gold_aa = full_seq[aa_index]
        else:
            # multiple solutions, go through them and select the best one
            max_score = 0
            index_of_max = -1
            for index, score in aln_solutions.iteritems():
                if score > max_score:
                    max_score = score
                    index_of_max = int(index)

            gold_aa = full_seq[index_of_max - 1]
        return gold_aa

    @staticmethod
    def int_or_gap(c):
        if c.isdigit():
            return int(c)
        else:
            return c

    @staticmethod
    def find_first_res(sequence):
        for i in sequence:
            if i != '-':
                return i

    def find_residues_neighbouring_insertions(self, gold_aln_cores, full_seqs):
        lowercase_res = {}
        for seq_id, seq_cores in gold_aln_cores.iteritems():
            #print seq_cores
            #seq_cores = map(self.int_or_gap, seq_cores)
            #lowercase_res[seq_id] = []
            lowercase_res[seq_id] = set()
            # if first residue is not one make it lowercase
            next_res = self.find_first_res(seq_cores)
            if next_res != 1:
                lowercase_res[seq_id].add(next_res)
            # if last residue is not the actual last make it lowercase
            next_res = self.find_first_res(seq_cores[::-1])
            if next_res != len(full_seqs[seq_id]):
                lowercase_res[seq_id].add(next_res)

            for i, curr_pos in enumerate(seq_cores):
                lower = False
                if curr_pos != "-":
                    last_res_pos = curr_pos
                else:
                    continue

                if i < len(seq_cores) - 1:
                    next_res = self.find_first_res(seq_cores[i + 1:])
                    # this is NOT the last residue
                    if next_res is not None and curr_pos + 1 != next_res:
                        print "after: comparing {} to {}, {}".format(
                            curr_pos + 1, seq_cores[i + 1], i + 1)
                        # there is insertion after this one, make lowercase
                        lowercase_res[seq_id].add(curr_pos)
                        lowercase_res[seq_id].add(next_res)
        return lowercase_res


    def aln_to_html_pairwise(self, aa_aln, gold_aln, full, wrong, order):
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"
        aln_length = len(aa_aln)

        _log.info("Creating pairiwse html")
        gold_lowercase_residues = self.find_residues_neighbouring_insertions(gold_aln["cores"], full)
        print gold_lowercase_residues
        for seq_id in order:
            if seq_id not in gold_aln["cores"]:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue

            seq = aa_aln[seq_id]
            html_sequence = "<b>TEST</b> {}    ".format(seq_id)
            html_gold_sequence = "<b>GOLD</b> {}    ".format(seq_id)
            asterisk_line = "         {}".format(" " * len(seq_id))
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
                    if gold_aa_index + 1 in gold_lowercase_residues[seq_id]:
                        print "change to lower"
                        gold_aa = gold_aa.lower()

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
                    if res == "-" and gold_aa == "-":
                        new_res = "<span class=noFeat>" + res + "</span>"
                        new_gold_res = "<span class=noFeat>" + gold_aa + "</span>"
                    elif res == "-":
                        new_res = "<span class=featMeh>" + res + "</span>"
                        new_gold_res = "<span class=featMeh>" + gold_aa + "</span>"
                    else:
                        raise RuntimeError("res: {}; gold res: {}".format(res, gold_aa))

                if res.upper() != gold_aa.upper() and new_res != "-":
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
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"

        for seq_id in order:
            seq = aa_aln[seq_id]
            html_sequence = "{}    ".format(seq_id)
            for r, res in enumerate(seq):
                if res != "-" and res != " ":
                    score = wrong[seq_id][r + 1]
                    if score[0] and score[1] == 1:
                        new_res = "<span class=featOK>{}</span>".format(res)
                    else:
                        level = self.get_level_cmplx(score[1])
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                            level, res)
                else:
                    new_res = "<span class=noFeat>" + res + "</span>"
                html_sequence += new_res
            html_out += html_sequence + "\n"
        return html_out

    def aln_to_html(self, aa_aln, wrong, order):
        # html_out = ""
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"
        aln_length = len(aa_aln)
        for seq_id in order:
            if seq_id not in wrong:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue

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
                    new_res = "<span class=noFeat>" + res + "</span>"
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
