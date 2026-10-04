# TAMA cigar_list / mapped_seq_length / trans_coordinates, copied so the
# oracle does not need Biopython. Behavior matches tama_collapse.py.
from __future__ import print_function

import re
import sys


def cigar_list(cigar):
    cigar = re.sub("=", "M", cigar)
    cigar = re.sub("X", "M", cigar)
    cig_char = re.sub(r"\d", " ", cigar)
    cig_char_list = cig_char.split()
    cig_digit = re.sub(r"[a-zA-Z]+", " ", cigar)
    cig_dig_list = cig_digit.split()
    return cig_dig_list, cig_char_list


def mapped_seq_length(cigar):
    cig_dig_list, cig_char_list = cigar_list(cigar)
    map_seq_length = 0
    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]
        if cig_flag == "H":
            continue
        elif cig_flag == "S":
            continue
        elif cig_flag == "M":
            map_seq_length = map_seq_length + int(cig_dig_list[i])
            continue
        elif cig_flag == "I":
            continue
        elif cig_flag == "D":
            map_seq_length = map_seq_length + int(cig_dig_list[i])
            continue
        elif cig_flag == "N":
            continue
    return map_seq_length


def trans_coordinates(start_pos, cigar):
    cig_dig_list, cig_char_list = cigar_list(cigar)
    end_pos = int(start_pos)
    exon_start_list = []
    exon_end_list = []
    exon_start_list.append(int(start_pos))
    for i in xrange(len(cig_dig_list)):
        cig_flag = cig_char_list[i]
        if cig_flag == "H":
            continue
        elif cig_flag == "S":
            continue
        elif cig_flag == "M":
            end_pos = end_pos + int(cig_dig_list[i])
            continue
        elif cig_flag == "I":
            continue
        elif cig_flag == "D":
            end_pos = end_pos + int(cig_dig_list[i])
            continue
        elif cig_flag == "N":
            exon_end_list.append(end_pos)
            end_pos = end_pos + int(cig_dig_list[i])
            exon_start_list.append(end_pos)
            continue
    exon_end_list.append(end_pos)
    return end_pos, exon_start_list, exon_end_list


CASES = [
    (1, "10M"),
    (100, "5S10M3N4M"),
    (50, "2H5M1I3M1D2M"),
    (10, "10=5X"),
    (20, "10M5S"),
    (1, "1M100N1M"),
    (8, "3S2M1I2D4N5M1H"),
    (1, "10P5M"),
]


def emit(start, cigar):
    length = mapped_seq_length(cigar)
    end_pos, starts, ends = trans_coordinates(start, cigar)
    print(
        "%d\t%s\t%d\t%d\t%s\t%s"
        % (
            start,
            cigar,
            length,
            end_pos,
            ",".join(str(x) for x in starts),
            ",".join(str(x) for x in ends),
        )
    )


if len(sys.argv) > 1:
    with open(sys.argv[1]) as handle:
        for raw in handle:
            if raw.startswith("@") or raw.strip() == "":
                continue
            fields = raw.rstrip("\n").split("\t")
            emit(int(fields[3]), fields[5])
else:
    for start, cigar in CASES:
        emit(start, cigar)
