def remove_positions(aln_dict, positions_to_remove):
    for i in sorted(positions_to_remove, reverse=True):
        for seq_id, seq in aln_dict.iteritems():
            if i < len(seq) - 1:
                new_seq = seq[:i-1] + seq[i-1].lower() + seq[i+1].lower() + seq[i+2:]
            else:
                new_seq = seq[:i-1] + seq[i-1].lower()
            aln_dict[seq_id] = new_seq


def find_positions_to_remove(master_aln_seq):
    """
    Find positions with gaps in the master sequence - these will be removed
    """
    positions_to_remove = []
    for i, res in enumerate(master_aln_seq):
        if res == "-":
            positions_to_remove.append(i)
    return positions_to_remove


def find_positions_to_fill_in(master_full_seq, master_aln_seq):
    """
    find on which positions do we need to add gaps
    """
    # position in full seq
    pos_full = 0
    # position in aligned seq
    pos_aln = 0

    # positions in the alignment on which gaps need to be added
    positions_to_fill_in = []
    while pos_full < len(master_full_seq):
        # if pos_aln < len(master_aln_seq) - 1 and master_full_seq[pos_full] != master_aln_seq[pos_aln].upper():
        if pos_aln < len(master_aln_seq) and master_full_seq[pos_full] != master_aln_seq[pos_aln].upper():
            # different residue on this position, we need to insert a gap here
            positions_to_fill_in.append(tuple([pos_aln, master_full_seq[pos_full]]))
            if master_aln_seq[pos_aln] == "-":
                pos_aln += 1
        # elif pos_aln >= len(master_aln_seq) - 1:
        elif pos_aln >= len(master_aln_seq):
            positions_to_fill_in.append(tuple([pos_aln, master_full_seq[pos_full]]))
        else:
            # same residue in the full and aligned sequence, can go to the next
            # position in both
            pos_aln += 1

        # moving along the full sequence with each iteration
        pos_full += 1
    return positions_to_fill_in


def fill_in_gaps(aln_dict, positions_to_fill_in, master_id, master_full_seq):
    # need to reverse positions order, to change the alignment starting
    # from C-terminus (then we don't need to update the alignment positions
    # to fill in)
    # rev_gap_positions = sorted(positions_to_fill_in, reverse=True, key=lambda x: x[0])
    rev_gap_positions = positions_to_fill_in[::-1]

    for pos_tuple in rev_gap_positions:
        # aln_position = pos_tuple[0] + 1
        aln_position = pos_tuple[0]
        master_res = pos_tuple[1]
        for seq_id, seq in aln_dict.iteritems():
            if seq_id == master_id:
                new_seq = seq[:aln_position] + master_res + seq[aln_position:]
            else:
                new_seq = seq[:aln_position] + "-" + seq[aln_position:]
            aln_dict[seq_id] = new_seq

    # there should be no residues missing in the master seq now,
    # we can make it all uppercase
    aln_dict[master_id] = aln_dict[master_id].upper()
    assert aln_dict[master_id].replace("-", "") == master_full_seq


def make_master_seq_full(aln_dict, full_seqs, gold_ids):
    """
    fill in the alignment with gaps so that the full master sequence is in
    the test alignment

    and remove positions on which master sequence has gaps
    """

    # retrieve data about the master sequence
    master_id = gold_ids[0]

    master_full_seq = full_seqs[master_id]
    master_aln_seq = aln_dict[master_id]

    # first check if it's needed
    if len(master_full_seq) == len(master_aln_seq) and "-" not in master_aln_seq:
        # alignment ok, master sequence is full and there are no gaps in it
        return aln_dict

    positions_to_remove = find_positions_to_remove(master_aln_seq)
    remove_positions(aln_dict, positions_to_remove)

    master_aln_seq = aln_dict[master_id]
    positions_to_fill_in = find_positions_to_fill_in(master_full_seq, master_aln_seq)

    fill_in_gaps(aln_dict, positions_to_fill_in, master_id, master_full_seq)

    return aln_dict
