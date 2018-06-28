fs = frozenset


BLOSUM = {
    fs(['B', 'N']): 3, fs(['W', 'L']): -2, fs(['G', 'G']): 6, fs(['X', 'S']): 0,
    fs(['X', 'D']): -1, fs(['K', 'G']): -2, fs(['S', 'E']): 0, fs(['X', 'M']): -1,
    fs(['Y', 'E']): -2, fs(['W', 'R']): -3, fs(['I', 'R']): -3, fs(['X', 'Z']): -1,
    fs(['H', 'E']): 0, fs(['V', 'M']): 1, fs(['N', 'R']): 0, fs(['I', 'D']): -3,
    fs(['F', 'D']): -3, fs(['W', 'C']): -2, fs(['N', 'A']): -2, fs(['W', 'Q']): -2,
    fs(['L', 'Q']): -2, fs(['S', 'N']): 1, fs(['Z', 'K']): 1, fs(['V', 'N']): -3,
    fs(['Q', 'N']): 0, fs(['M', 'K']): -1, fs(['V', 'H']): -3, fs(['G', 'E']): -2,
    fs(['S', 'L']): -2, fs(['P', 'R']): -2, fs(['D', 'A']): -2, fs(['S', 'C']): -1,
    fs(['E', 'D']): 2, fs(['Y', 'G']): -3, fs(['W', 'P']): -4, fs(['X', 'X']): -1,
    fs(['Z', 'L']): -3, fs(['Q', 'A']): -1, fs(['V', 'Y']): -1, fs(['W', 'A']): -3,
    fs(['G', 'D']): -1, fs(['X', 'P']): -2, fs(['K', 'D']): -1, fs(['T', 'N']): 0,
    fs(['Y', 'F']): 3, fs(['W', 'W']): 11, fs(['Z', 'M']): -1, fs(['L', 'D']): -4,
    fs(['M', 'R']): -1, fs(['Y', 'K']): -2, fs(['F', 'E']): -3, fs(['M', 'E']): -2,
    fs(['S', 'S']): 4, fs(['X', 'C']): -2, fs(['Y', 'L']): -1, fs(['H', 'R']): 0,
    fs(['P', 'P']): 7, fs(['K', 'C']): -3, fs(['S', 'A']): 1, fs(['P', 'I']): -3,
    fs(['Q', 'Q']): 5, fs(['L', 'I']): 2, fs(['P', 'F']): -4, fs(['B', 'A']): -2,
    fs(['Z', 'N']): 0, fs(['M', 'Q']): 0, fs(['V', 'I']): 3, fs(['Q', 'C']): -3,
    fs(['I', 'H']): -3, fs(['Z', 'D']): 1, fs(['Z', 'P']): -1, fs(['Y', 'W']): 2,
    fs(['T', 'G']): -2, fs(['B', 'P']): -2, fs(['P', 'A']): -1, fs(['C', 'D']): -3,
    fs(['Y', 'H']): 2, fs(['X', 'V']): -1, fs(['B', 'B']): 4, fs(['Z', 'F']): -3,
    fs(['M', 'L']): 2, fs(['F', 'G']): -3, fs(['S', 'M']): -1, fs(['M', 'G']): -3,
    fs(['Z', 'Q']): 3, fs(['S', 'Q']): 0, fs(['X', 'A']): 0, fs(['V', 'T']): 0,
    fs(['W', 'F']): 1, fs(['S', 'H']): -1, fs(['X', 'N']): -1, fs(['B', 'Q']): 0,
    fs(['K', 'A']): -1, fs(['I', 'Q']): -3, fs(['X', 'W']): -2, fs(['N', 'N']): 6,
    fs(['W', 'T']): -2, fs(['P', 'D']): -1, fs(['B', 'C']): -3, fs(['I', 'C']): -1,
    fs(['V', 'K']): -2, fs(['X', 'Y']): -1, fs(['K', 'R']): 2, fs(['Z', 'R']): 0,
    fs(['W', 'E']): -3, fs(['T', 'E']): -1, fs(['B', 'R']): -1, fs(['L', 'R']): -2,
    fs(['Q', 'R']): 1, fs(['X', 'F']): -1, fs(['T', 'S']): 1, fs(['B', 'D']): 4,
    fs(['Z', 'A']): -1, fs(['M', 'N']): -2, fs(['V', 'D']): -3, fs(['F', 'A']): -2,
    fs(['X', 'E']): -1, fs(['F', 'H']): -1, fs(['M', 'A']): -1, fs(['K', 'Q']): 1,
    fs(['Z', 'S']): 0, fs(['X', 'G']): -1, fs(['V', 'V']): 4, fs(['W', 'D']): -4,
    fs(['X', 'H']): -1, fs(['S', 'F']): -2, fs(['X', 'L']): -1, fs(['B', 'S']): 0,
    fs(['S', 'G']): 0, fs(['P', 'M']): -2, fs(['Y', 'M']): -1, fs(['H', 'D']): -1,
    fs(['B', 'E']): 1, fs(['Z', 'B']): 1, fs(['I', 'E']): -3, fs(['V', 'E']): -2,
    fs(['X', 'T']): 0, fs(['X', 'R']): -1, fs(['R', 'R']): 5, fs(['Z', 'T']): -1,
    fs(['Y', 'D']): -3, fs(['V', 'W']): -3, fs(['F', 'L']): 0, fs(['T', 'C']): -1,
    fs(['X', 'Q']): -1, fs(['B', 'T']): -1, fs(['K', 'N']): 0, fs(['T', 'H']): -2,
    fs(['Y', 'I']): -1, fs(['F', 'Q']): -3, fs(['T', 'I']): -1, fs(['T', 'Q']): -1,
    fs(['P', 'L']): -3, fs(['R', 'A']): -1, fs(['B', 'F']): -3, fs(['Z', 'C']): -3,
    fs(['M', 'H']): -2, fs(['V', 'F']): -1, fs(['F', 'C']): -2, fs(['L', 'L']): 4,
    fs(['M', 'C']): -1, fs(['C', 'R']): -3, fs(['D', 'D']): 6, fs(['E', 'R']): 0,
    fs(['V', 'P']): -2, fs(['S', 'D']): 0, fs(['E', 'E']): 5, fs(['W', 'G']): -2,
    fs(['P', 'C']): -3, fs(['F', 'R']): -3, fs(['B', 'G']): -1, fs(['C', 'C']): 9,
    fs(['I', 'G']): -4, fs(['V', 'G']): -3, fs(['W', 'K']): -3, fs(['G', 'N']): 0,
    fs(['I', 'N']): -3, fs(['Z', 'V']): -2, fs(['A', 'A']): 4, fs(['V', 'Q']): -2,
    fs(['F', 'K']): -3, fs(['T', 'A']): 0, fs(['B', 'V']): -3, fs(['K', 'L']): -2,
    fs(['L', 'N']): -3, fs(['Y', 'N']): -2, fs(['F', 'F']): 6, fs(['L', 'G']): -4,
    fs(['B', 'H']): 0, fs(['Z', 'E']): 4, fs(['Q', 'D']): 0, fs(['X', 'B']): -1,
    fs(['Z', 'W']): -3, fs(['S', 'K']): 0, fs(['X', 'K']): -1, fs(['V', 'R']): -3,
    fs(['K', 'E']): 1, fs(['I', 'A']): -1, fs(['P', 'H']): -2, fs(['B', 'W']): -4,
    fs(['K', 'K']): 5, fs(['H', 'C']): -3, fs(['E', 'N']): 0, fs(['Y', 'Q']): -1,
    fs(['H', 'H']): 8, fs(['B', 'I']): -3, fs(['C', 'A']): 0, fs(['I', 'I']): 4,
    fs(['V', 'A']): 0, fs(['W', 'I']): -3, fs(['T', 'F']): -2, fs(['V', 'S']): -2,
    fs(['T', 'T']): 5, fs(['F', 'M']): 0, fs(['L', 'E']): -3, fs(['M', 'M']): 5,
    fs(['Z', 'G']): -2, fs(['D', 'R']): -2, fs(['M', 'D']): -3, fs(['W', 'H']): -2,
    fs(['G', 'C']): -3, fs(['S', 'R']): -1, fs(['S', 'I']): -2, fs(['P', 'Q']): -1,
    fs(['Y', 'A']): -2, fs(['X', 'I']): -1, fs(['E', 'A']): -1, fs(['B', 'Y']): -3,
    fs(['K', 'I']): -3, fs(['H', 'A']): -2, fs(['P', 'G']): -2, fs(['F', 'N']): -3,
    fs(['H', 'N']): 1, fs(['B', 'K']): 0, fs(['V', 'C']): -1, fs(['T', 'L']): -1,
    fs(['P', 'K']): -1, fs(['W', 'S']): -3, fs(['T', 'D']): -1, fs(['T', 'M']): -1,
    fs(['P', 'N']): -2, fs(['K', 'H']): -1, fs(['T', 'R']): -1, fs(['Y', 'R']): -2,
    fs(['L', 'C']): -1, fs(['B', 'L']): -4, fs(['Z', 'Y']): -2, fs(['W', 'N']): -4,
    fs(['G', 'A']): 0, fs(['S', 'P']): -1, fs(['E', 'Q']): 2, fs(['C', 'N']): -3,
    fs(['H', 'Q']): 0, fs(['D', 'N']): 1, fs(['Y', 'C']): -2, fs(['L', 'H']): -3,
    fs(['E', 'C']): -4, fs(['Z', 'H']): 0, fs(['H', 'G']): -2, fs(['P', 'E']): -1,
    fs(['Y', 'S']): -2, fs(['G', 'R']): -2, fs(['B', 'M']): -3, fs(['Z', 'Z']): 4,
    fs(['W', 'M']): -1, fs(['Y', 'T']): -2, fs(['Y', 'P']): -3, fs(['Y', 'Y']): 7,
    fs(['T', 'K']): -1, fs(['Z', 'I']): -3, fs(['T', 'P']): -1, fs(['V', 'L']): 1,
    fs(['F', 'I']): 0, fs(['G', 'Q']): -2, fs(['L', 'A']): -1, fs(['M', 'I']): 1}

# set of all residue codes that are present in the BLOSUM matrix
ALLOWED_CODES = set.union(*(map(set, BLOSUM)))


def get_blosum_score(res1, res2):
    if res1 not in ALLOWED_CODES:
        res1 = 'X'
    if res2 not in ALLOWED_CODES:
        res2 = 'X'

    return BLOSUM[fs([res1, res2])]
