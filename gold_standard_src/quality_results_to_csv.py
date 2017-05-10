"""
Input format:

program, parameters, results paths
"""
import argparse
import re
from sys import argv


def parse_files(paths):
    results = {}
    uniq_tmpls = set()
    for batch, p in enumerate(paths):
        with open(p) as a:
            result_lines = a.read().replace(', ', ',').splitlines()[3:]
        # get unique aligned templates:
        for l in result_lines:
            if l.startswith("# "):
                uniq_tmpls.add(l.split()[1])
                uniq_tmpls.add(l.split()[2])

        results[batch] = {
            'sensitivity': result_lines[0].split(":")[1],
            'specificity': result_lines[4].split(":")[1],
            'mcc': result_lines[2].split(":")[1],
            'npv': result_lines[1].split(":")[1],
            'ppv': result_lines[3].split(":")[1],
            'aligned templates': result_lines[6].split(":")[1]
        }
    return results, uniq_tmpls


def get_string(s):
    int_float_reg = r'^(?=.*\d)\d*(?:\.\d+)?$'
    reg = re.compile(int_float_reg)
    if isinstance(s, float):
        return "%.2f" % s
    elif isinstance(s, int):
        return str(s)
    elif reg.match(s.replace(' ', '')):
        return "%.2f" % float(s)
    else:
        return s


def get_line(values, key):
    newline = ["", "", key]
    for i in values:
        newline.append(values[i][key])
    return newline


def convert_to_csv(inpath, outpath):
    with open(inpath) as a:
        filedata = []
        for i in a.read().splitlines():
            filedata.append((i.split(',')[0], i.replace(', ', ',').split(",")))

    header = "program, transition modifiers, results, batch 1, batch 2, batch 3, batch 4, averages, total aligned templates"
    csv_lines = [header]
    for program, info in filedata:
        paths = info[2:]
        batches_no = len(paths)
        program = info[0]
        params = info[1]

        r, uniq_tmpls = parse_files(paths)
        newlines = [
            get_line(r, "sensitivity"),
            get_line(r, "specificity"),
            get_line(r, "ppv"),
            get_line(r, "npv"),
            get_line(r, "mcc"),
            get_line(r, "aligned templates")
        ]
        for i in range(len(newlines)):
            av = sum(map(float, newlines[i][3:])) / batches_no
            newlines[i].append(av)
        newlines = map(lambda x: map(get_string, x), newlines)
        newlines = map(lambda x: ",".join(x), newlines)

        first_line = "{},{},,,,,,,{}".format(program, params, len(uniq_tmpls))
        new_lines = [first_line] + newlines
        csv_lines.extend(new_lines)

    with open(outpath, 'w') as o:
        o.write("\n".join(csv_lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("inpath")
    parser.add_argument("outpath")
    args = parser.parse_args()

    convert_to_csv(args.inpath, args.outpath)
