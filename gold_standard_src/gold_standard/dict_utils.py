def merge_dicts(dict1, dict2):
    """
    Sum up values with the same keys
    """
    return {k: dict1.get(k, 0) + dict2.get(k, 0)
            for k in set(dict1) | set(dict2)}


def merge_nested_dicts(dict1, dict2):
    """
    Add elements nested in dicts of dicts
    e.g. :
        {"a": {"1": 0, "2": 1}}
        +
        {"a": {"3": 3}}

        =
        {"a": {"1": 0, "2": 1, "3": 3}}
    Will not change the values inside the nested dicts, keys in the
    nested dicts cannot be duplicate
    You expect that all of the outer keys in dict1 are also present in dict2
    (and the other way round)
    """
    for outer_key, inner_dict1 in dict1.iteritems():
        inner_dict2 = dict2[outer_key]

        for inner_key, val in inner_dict2.iteritems():
            inner_dict1[inner_key] = val

    return dict1
