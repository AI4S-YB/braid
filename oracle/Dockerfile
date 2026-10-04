# Frozen TAMA oracle. Digest recorded after the image was pulled:
# python:2.7.18@sha256:cfa62318c459b1fde9e0841c619906d15ada5910d625176e24bf692cf8a2601d
FROM python:2.7.18@sha256:cfa62318c459b1fde9e0841c619906d15ada5910d625176e24bf692cf8a2601d

# SeqIO only. Biopython 1.76 is the last release that supports Python 2.
# Debian Buster's apt repos are archived; install from PyPI and skip the system index.
RUN pip install --no-cache-dir 'biopython==1.76'

WORKDIR /work
ENV PYTHONHASHSEED=0
