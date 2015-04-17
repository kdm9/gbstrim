#!/usr/bin/env python

# Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import skbio
from screed.fastq import fastq_iter
from itertools import izip
import sys

SSW = skbio.alignment.StripedSmithWaterman

def pairitr(reads):
    pair = [None, None]
    for rnum, read in enumerate(reads):
        if rnum % 2 == 0:
            pair[0] = read
        else:
            pair[1] = read
            yield pair

def printrecs(r1, r2, stream=sys.stdout):
    for r in [r1, r2]:
        stream.write(
            "@{name} {annotations}\n{sequence}\n+\n{quality}\n".format(**r))

def is_dimer(r1, r2):
    """Excuse the magic numbers. """
    r1seq = r1["sequence"]
    r2seq = r2["sequence"]
    adapt = SSW("AGATCGGAAGAG")
    r1pos = adapt(r1seq[6:20]).target_begin
    r2pos = adapt(r2seq[6:20]).target_begin
    if r1pos-4 == r2pos-5 and r1pos-4 < 20:
        printrecs(r1, r2, stream=sys.stderr)
        return True
    return False

def main_ilfq(ilfq_name):
    pst1 = SSW("CTGCA")
    fh = open(ilfq_name)
    reads = fastq_iter(fh)
    for rpair in pairitr(reads):
        r1, r2 =  rpair
        r1aln = pst1(r1["sequence"][5:])
        r2aln = pst1(r2["sequence"][5:])
        if r1aln.query_begin != 0 or r1aln.query_end != 4 or \
                r2aln.query_begin != 0 or r2aln.query_end != 4:
            # Not read-through, but might be dimer.
            if not is_dimer(r1, r2):
                printrecs(r1, r2)
        if r1aln.target_begin == r2aln.target_begin:
            r1["sequence"] = r1["sequence"][:r1aln.target_begin]
            r1["quality"] = r1["quality"][:r1aln.target_begin]
            r2["sequence"] = "N"
            r2["quality"] = "#"
            printrecs(r1, r2)
    fh.close()

if __name__ == "__main__":
    main_ilfq(sys.argv[1])
