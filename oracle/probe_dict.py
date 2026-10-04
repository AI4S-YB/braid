# Probe CPython 2.7 dict iteration versus insertion order.
# Run under PYTHONHASHSEED=0 and under another seed.
from __future__ import print_function

import os
import sys


def line(kind, *parts):
    print("\t".join([kind] + list(parts)))


def as_text(keys):
    return ",".join(repr(k) for k in keys)


def report(label, keys):
    mapping = {}
    for key in keys:
        mapping[key] = 1
    order = list(mapping)
    line(
        "CASE",
        label,
        "1" if order == list(keys) else "0",
        as_text(keys),
        as_text(order),
    )
    line("HASH", label, ",".join(str(hash(k)) for k in keys))


line("VERSION", sys.version.replace("\n", " "))
line("HASHSEED", os.environ.get("PYTHONHASHSEED", "<unset>"))

read_ids = [
    "rec_c99144/3/1053",
    "rec_c104542/4/1069",
    "2_c312362/1/1255",
    "G1.1",
    "G1.2",
    "G10.1",
    "11_c95610/24/952",
]
report("read_ids", read_ids)
report("priorities", [3, 1, 0, 2])
report("coords", [900, 867, 1200, 1053, 10, 20])
report("gene_nums", [1, 2, 10, 3, 11])
report("errors", ["0", "na", "10M>3.A", "5I_2M", "aa-bb"])
report("many", ["k%04d" % i for i in range(64)])
report("bigints", [34210, 867, 100000, 2 ** 31 - 1, 0, -1])

updated = {}
for key in read_ids:
    updated[key] = 1
updated[read_ids[0]] = 2
updated["new_id"] = 1
line("UPDATE", as_text(list(updated)))

popped = {}
for key in read_ids:
    popped[key] = 1
popped.pop(read_ids[1])
line("POP", as_text(list(popped)))

# Nested dict, same shape as priority_error_start_dict[priority][coord][error] = 1
nested = {}
for priority, coord, err in (
    (1, 100, "na"),
    (1, 100, "3.A"),
    (0, 90, "0"),
    (1, 110, "na"),
    (2, 80, "5I"),
):
    if priority not in nested:
        nested[priority] = {}
    if coord not in nested[priority]:
        nested[priority][coord] = {}
    nested[priority][coord][err] = 1
line("NEST_PRI", as_text(list(nested)))
line("NEST_COORD_1", as_text(list(nested[1])))
line("NEST_ERR", as_text(list(nested[1][100])))
