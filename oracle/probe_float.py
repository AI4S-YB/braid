# Python 2 str(round(x, 2)), the format TAMA writes for coverage, identity, and polyA.
from __future__ import print_function


def emit(kind, expr, value):
    rounded = round(value, 2)
    print("%s\t%s\t%s\t%s\t%s" % (kind, expr, repr(value), repr(rounded), str(rounded)))


literals = [
    0,
    0.0,
    1.0,
    100.0,
    99.0,
    99.5,
    99.995,
    99.994,
    99.996,
    85.0,
    85.125,
    85.135,
    85.126,
    85.124,
    70.0,
    50.0,
    2.5,
    3.5,
    1.005,
    1.015,
    1.225,
    8.195,
    0.005,
    0.015,
    10.0 / 3.0,
    1.0 / 3.0,
    -1.0,
    -1.5,
    -2.5,
]
for value in literals:
    emit("LIT", repr(value), value)

for num, den in (
    (100, 100),
    (99, 100),
    (999, 1000),
    (85, 100),
    (1, 3),
    (2, 3),
    (10, 7),
    (22, 7),
    (1000, 3),
    (1, 6),
    (5, 6),
    (1, 1),
    (0, 1),
    (123456, 1000),
):
    emit("DIV", "%d/%d*100" % (num, den), float(num) / float(den) * 100.0)

for i in range(0, 21):
    emit("HALF", repr(i + 0.5), i + 0.5)
    emit("MILLI", repr(i / 10.0 + 0.005), i / 10.0 + 0.005)
