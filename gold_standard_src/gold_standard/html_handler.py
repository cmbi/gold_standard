import copy
import logging
import os

from .aln_analyzer import get_score_mod_value

from .paths import TEMPLATE


_log = logging.getLogger("__main__")


class HtmlHandler(object):
    def __init__(self, long_len=20):
        self.long_len = long_len

    def write_html(self, quality_data, outname, mode="cores"):
        if mode in ["var", "var_short"]:
            outtxt = self.aln_to_html_var(quality_data, mode)
        elif mode in ["var_complex", "var_short_complex"]:
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
                    quality_data['aa_aln'], quality_data["num_aln"],
                    quality_data["gold_aln"], quality_data['wrong_cols'],
                    quality_data["order"], quality_data["target_id"],
                    quality_data["full"])
        else:
            return

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
        .featWRONG0{
            font-family: monospace;
            background: #FFE618;
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

    @staticmethod
    def get_insertion_positions(num_aln):
        insertion_positions = set()
        return insertion_positions

    def make_corvar_new(self, full_seqs, num_aln):
        insertion_positions = self.get_insertion_positions(num_aln)

    def aln_to_html_var(self, quality_data, mode):
        short_var = mode in ["var_short", "var_short_complex"]
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"
        aln_length = len(quality_data['aa_aln'])
        # quality_data["core_indexes"] = [0, 4]
        num_aln_c = self.split_cores(quality_data['num_aln'],
                                     quality_data['core_indexes'])
        num_aln_v = self.split_vars(num_aln_c)

        aa_aln_corvar = self.make_corvar(quality_data['full'], num_aln_v)
        var_lengths = self.get_max_var_lengths(num_aln_v, short_var)
        #aa_aln_corvar, num_aln_v = self.make_corvar_new(quality_data['full'], quality_data["num_aln"])
        #var_lengths = self.get_max_var_lengths_new(num_aln_v, short_var)

        for seq_id in quality_data["order"]:
            if seq_id not in quality_data["gold_ids"]:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue
            seq = aa_aln_corvar[seq_id]
            if mode.endswith("complex"):
                html_seq = self.make_html_var_seq_complex(
                        seq, quality_data['wrong_cols'][seq_id], var_lengths,
                        aln_length, short_var)
            else:
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

    @staticmethod
    def get_full_seq_pos(merged_corvar_seq, core_index):
        c = 0
        gaps = 0
        for pos, res in enumerate(merged_corvar_seq):
            if res == "-":
                gaps += 1
                c += 1
            if not res.islower() and res != "-":
                if c == core_index:
                    return pos - gaps
                c += 1
        raise RuntimeError("Did not find full seq pos for core index %d in seq: %s" % (core_index, merged_corvar_seq))

    def make_html_var_seq_complex(self, corvar_seq, wrong, max_lengths, aln_length, short_var=False):
        html_seq = ""
        r_index = 0
        corvar_merged = ""
        for c, core in enumerate(corvar_seq["cores"]):
            var = corvar_seq["var"][c].lower()
            full_var = var
            if short_var:
                var = self.make_short_var(var)
            corvar_merged += full_var + core
            var = " " + var + " " * (max_lengths[c] - len(var)) + " "
            html_seq += var
            for i, res in enumerate(core):
                if res == "-":
                    new_res = "<span>-</span>"
                else:
                    full_seq_pos = self.get_full_seq_pos(corvar_merged, r_index)
                    score = wrong[full_seq_pos + 1]
                    if score[0] and score[1] == 1:
                        new_res = "<span class=featOK>{}</span>".format(res)
                    else:
                        level = self.get_level_cmplx(score[1])
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                                level, res)
                r_index += 1
                html_seq += new_res
        var = corvar_seq["var"][-1].lower()
        if short_var:
            var = self.make_short_var(var)
        html_seq += " " + var + " " * (max_lengths[-1] - len(var))
        return html_seq

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

    @staticmethod
    def make_html_target_sequence(target_seq, target_id):
        header = "<span style='color:gray;'><b>    TARGET</b> {}   </span>".format(target_id)
        sequence = ""
        for i, res_i in enumerate(target_seq, start=1):
            if i % 10 == 0:
                res_i = "<b>" + res_i + "</b>"
                res_i = "<span class=noFeat style='color:black;'>{}</span>".format(res_i)
            else:
                res_i = "<span class=noFeat style='color:gray;'>{}</span>".format(res_i)
            sequence += res_i
        html_target_sequence = header + sequence
        return html_target_sequence

    @staticmethod
    def make_ruler(target_sequence, spaces_no=19):
        ruler = " " * spaces_no
        ticks = " " * spaces_no
        for i in range(1, len(target_sequence) + 1):
            if i % 10 == 0:
                new_number = " " * (10 - len(str(i))) + str(i)
                new_number = "<span class=noFeat style='color:gray;'>{}</span>".format(new_number)
                ruler += new_number
                new_tick = "<span class=noFeat style='color:gray;'>{}</span>".format(9 * " " + "|")
                ticks += new_tick

        return ruler, ticks

    def aln_to_html_pairwise_complex(self, quality_data):
        num_aln = quality_data["num_aln"]
        aa_aln = quality_data["aa_aln"]
        gold_aln = quality_data["gold_aln"]
        full = quality_data["full"]
        wrong = quality_data["wrong_cols"]
        order = quality_data["order"]
        target_id = quality_data["target_id"]
        gold_corvar = quality_data["gold_corvar"]["alns"]

        target_seq = full[target_id]
        gold_lowercase_residues = self.find_residues_neighbouring_insertions(gold_corvar["cores"], full)

        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"

        _log.info("Creating pairiwse html")
        master_num_seq = num_aln["cores"][target_id]
        html_target_sequence = self.make_html_target_sequence(target_seq, target_id)

        ruler, ticks = self.make_ruler(target_seq)
        html_out += ruler + "\n"
        html_out += ticks + "\n"

        for i, seq_id in enumerate(order, start=1):
            if seq_id not in gold_aln:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue
            number = " " * (3 - len(str(i))) + str(i)
            html_sequence = "<b>{} TEST   {}</b>   ".format(number, seq_id)
            html_gold_sequence = "<b>    GOLD   {}</b>   ".format(seq_id)
            asterisk_line = "<span class=asteriskBlank>         {}</span>".format(" " * len(seq_id))

            pairwise_gold_aln = gold_aln[seq_id]
            num_seq = num_aln["cores"][seq_id]
            for aln_pos, full_seq_pos in enumerate(num_seq):
                # full seq pos is 1-based, aln pos is 0-based, will change
                # full_seq pos to 0-based right below
                if full_seq_pos != "-":
                    full_seq_pos = int(full_seq_pos) - 1
                    # res = full[seq_id][full_seq_pos]
                    res = aa_aln[seq_id][aln_pos]
                    score = wrong[seq_id][full_seq_pos + 1]
                else:
                    res = "-"
                    score = 0

                master_index = int(master_num_seq[aln_pos])

                # gold_aa_full_seq_pos is also 0-based just like full_seq_pos
                gold_aa, gold_aa_full_seq_pos = self.get_gold_aa(pairwise_gold_aln, full[seq_id], master_index)

                if gold_aa != "-" and gold_aa_full_seq_pos + 1 in gold_lowercase_residues[seq_id]:
                    gold_aa = gold_aa.lower()

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
                        master_score = self.get_master_score(master_index, gold_aln[seq_id])
                        level = self.get_level_cmplx(master_score)
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                                level, res)
                        new_gold_res = "<span class=featWRONG{}>{}</span>".format(
                                level, gold_aa)
                        add_asterisk = True
                    else:
                        raise RuntimeError("res: {}; gold res: {}".format(res, gold_aa))
                if add_asterisk:
                    asterisk_line += "<span class=asteriskFull>*</span>"
                else:
                    asterisk_line += "<span class=asteriskBlank> </span>"

                html_sequence += new_res
                html_gold_sequence += new_gold_res
            # html_out += asterisk_line + "\n"
            html_out += html_target_sequence + "\n"
            html_out += html_sequence + "\n"
            html_out += html_gold_sequence + "\n"
            html_out += "<br>"
        html_out += "</div>"
        _log.info("Finished creating pairwise html")
        return html_out

    @staticmethod
    def get_master_score(master_index, gold_aln):
        """
        Find what's the score for alignment of this residue from the master sequence
        (this is for false negatives)
        :param master_index: 1-based position in the mnaster sequence
        :param gold_aln: pairwise alignment (like in the final_core.json file)
        :return:
        """
        # this master residue might have multiple alns allowed, find all then take the one
        # with the highest score
        highest_score = -1000
        found = False
        for res_i, alns in gold_aln.iteritems():
            if str(master_index) in alns:
                score = get_score_mod_value(alns[str(master_index)])
                if score > highest_score:
                    highest_score = score
                    found = True
        if not found:
            raise RuntimeError("could not find aln for master residue {}".format(highest_score))

        return 1 - highest_score

    @ staticmethod
    def get_gold_aa(pairwise_gold_aln, full_seq, pos):
        """
        Get residue in the gold aln that's aligned to residue 'pos'
        """
        aln_solutions = {}
        gold_aa_full_seq_pos = None
        for res, res_alns in pairwise_gold_aln.iteritems():
            if str(pos) in res_alns:
                aln_solutions[res] = res_alns[str(pos)]

        if not aln_solutions:
            # this residue should not be aligned to anything
            gold_aa = "-"
        elif len(aln_solutions) == 1:
            # there is only one solution, take it
            tmp_index = aln_solutions.keys()[0]
            if tmp_index != "-":
                aa_index = int(tmp_index) - 1
                gold_aa = full_seq[aa_index]
                gold_aa_full_seq_pos = aa_index
            else:
                gold_aa = "-"

        else:
            # multiple solutions, go through them and select the best one
            max_score = 0
            index_of_max = -1
            for index, score in aln_solutions.iteritems():
                if score > max_score:
                    max_score = score
                    index_of_max = int(index)

            gold_aa = full_seq[index_of_max - 1]
            gold_aa_full_seq_pos = index_of_max - 1
        return gold_aa, gold_aa_full_seq_pos

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
                        # there is insertion after this one, make lowercase
                        lowercase_res[seq_id].add(curr_pos)
                        lowercase_res[seq_id].add(next_res)
        return lowercase_res

    def aln_to_html_pairwise(self, aa_aln, gold_aln, full, wrong, order):
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"
        aln_length = len(aa_aln)

        _log.info("Creating pairiwse html")
        gold_lowercase_residues = self.find_residues_neighbouring_insertions(gold_aln["cores"], full)
        for seq_id in order:
            if seq_id not in gold_aln["cores"]:
                _log.warning("Sequence %s from the test aln is not present in the gold aln", seq_id)
                continue

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
                    if gold_aa_index + 1 in gold_lowercase_residues[seq_id]:
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

    def complex_aln_to_html(self, aa_aln, num_aln, gold_aln, wrong, order, target_id, full):
        html_out = "<div class=monospacediv style='font-family:monospace;'>\n<br>"

        ruler, ticks = self.make_ruler(aa_aln[target_id], spaces_no=9)
        html_out += ruler + "\n"
        html_out += ticks + "\n"

        master_num_seq = num_aln["cores"][target_id]
        for seq_id in order:
            seq = aa_aln[seq_id]
            html_sequence = "{}    ".format(seq_id)
            num_seq = num_aln["cores"][seq_id]
            pairwise_gold_aln = gold_aln[seq_id]
            for r, res in enumerate(seq):
                if res != "-" and res != " ":
                    real_pos = num_seq[r]

                    score = wrong[seq_id][real_pos]
                    if score[0] and score[1] == 1:
                        new_res = "<span class=featOK>{}</span>".format(res)
                    else:
                        level = self.get_level_cmplx(score[1])
                        new_res = "<span class=featWRONG{}>{}</span>".format(
                            level, res)
                elif res == "-":
                    master_residue = master_num_seq[r]
                    if master_residue == "-":
                        new_res = "<span class=noFeat>" + res + "</span>"
                    else:
                        master_index = int(master_residue)
                        gold_aa, gold_aa_full_seq_pos = self.get_gold_aa(pairwise_gold_aln, full[seq_id], master_index)
                        if gold_aa == "-":
                            new_res = "<span class=noFeat>" + res + "</span>"
                        else:
                            master_score = self.get_master_score(master_index, gold_aln[seq_id])
                            level = self.get_level_cmplx(master_score)
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

        ruler, ticks = self.make_ruler(aa_aln[order[0]], spaces_no=9)
        html_out += ruler + "\n"
        html_out += ticks + "\n"

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

        if n > 0.8:
            return 5
        elif n > 0.5:
            return 4
        elif n > 0.4:
            return 3
        elif n > 0.2:
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
