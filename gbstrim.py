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
import docopt
import gzip

SSW = skbio.alignment.StripedSmithWaterman

CLI_ARGS = """
USAGE: gbstrim -e RE_SITE <input_ILFQ>

OPTIONS:
    -e RE_SITE      Restriction site sequence
"""

LONGEST_BARCODE=8
ADAPTOR = "AGATCGGAAGAG"
ADAPTOR_LEN = len(ADAPTOR)

RE = {
    'PST1': {'site': "CTGCAG",
             'overhang_idx': 1},
}

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

def is_dimer(r1, r2, re_len):
    """
    Checks if the start of the adaptor is in the first few bases of both reads.
    """
    adapt = SSW(ADAPTOR)
    r1seq = r1["sequence"]
    r2seq = r2["sequence"]
    align_until = re_len + ADAPTOR_LEN + LONGEST_BARCODE
    r1aln = adapt(r1seq[re_len:align_until])
    r2aln = adapt(r2seq[re_len:align_until])
    r1len = len(r1aln.aligned_query_sequence)
    r2len = len(r2aln.aligned_query_sequence)
    if r1len != ADAPTOR_LEN or r2len != ADAPTOR_LEN:
        return False
    r1pos = r1aln.target_begin
    r2pos = r2aln.target_begin
    # If both alignments are within the alignment region, then we can only
    # assume there's dimer or a stupidly small fragment.
    if r1pos-re_len < align_until and r2pos-re_len < align_until:
        # printrecs(r1, r2, stream=sys.stderr)
        return True
    return False

def main_ilfq(ilfq_name, re_site):
    pst1 = SSW(re_site)
    re_len = len(re_site);
    if ilfq_name.endswith("gz"):
        fh = gzip.open(ilfq_name)
    else:
        fh = open(ilfq_name)
    reads = fastq_iter(fh)
    for rpair in pairitr(reads):
        r1, r2 =  rpair
        r1aln = pst1(r1["sequence"][re_len:])
        r2aln = pst1(r2["sequence"][re_len:])
        if r1aln.query_begin != 0 or r1aln.query_end != 4 or \
                r2aln.query_begin != 0 or r2aln.query_end != 4:
            # Not read-through, but might be dimer.
            if not is_dimer(r1, r2, re_len):
                printrecs(r1, r2)
        if r1aln.target_begin == r2aln.target_begin:
            r1["sequence"] = r1["sequence"][:r1aln.target_begin]
            r1["quality"] = r1["quality"][:r1aln.target_begin]
            r2["sequence"] = "N"
            r2["quality"] = "#"
            printrecs(r1, r2)
    fh.close()

if __name__ == "__main__":
    opts = docopt.docopt(CLI_ARGS)
    main_ilfq(opts['<input_ILFQ>'], opts['-e'])
