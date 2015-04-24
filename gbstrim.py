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

from __future__ import print_function
import skbio
from screed.fastq import fastq_iter
from itertools import izip
import sys
import docopt
import gzip

SSW = skbio.alignment.StripedSmithWaterman

CLI_ARGS = """
USAGE: gbstrim -e RE_SITE [-d DIMER_FILE]  <input_ILFQ>

OPTIONS:
    -e RE_SITE      Restriction site sequence
    -d DIMER_FILE   Keep dimers in file
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
        print("@{name} {annotations}".format(**r), file=stream)
        print("{sequence}".format(**r), file=stream)
        print("+", file=stream)
        print("{quality}".format(**r), file=stream)

def is_dimer(r1, r2, re_len, dimer_file=None):
    """
    Checks if the start of the adaptor is in the first few bases of both reads.
    """
    adapt = SSW(ADAPTOR)
    r1seq = r1["sequence"]
    r2seq = r2["sequence"]
    align_until = re_len + ADAPTOR_LEN + LONGEST_BARCODE
    aln = adapt(r1seq[re_len:align_until])
    aln_len = len(aln.aligned_query_sequence)
    if aln_len / float(ADAPTOR_LEN) < 0.9:
        # adaptor alignment isn't complete
        return False
    r1pos = aln.target_begin
    if r1pos - re_len < align_until:
        if dimer_file is not None:
            printrecs(r1, r2, stream=dimer_file)
        return True
    return False

def main_ilfq(ilfq_name, re_site, dimer_file=None):
    pst1 = SSW(re_site)
    if dimer_file is not None:
        dimer_file = open(dimer_file, 'w')
    n_reads = 0
    n_adapt = 0
    n_trimmed = 0
    re_len = len(re_site);
    if ilfq_name.endswith("gz"):
        fh = gzip.open(ilfq_name)
    else:
        fh = open(ilfq_name)
    reads = fastq_iter(fh)
    for rpair in pairitr(reads):
        n_reads += 1
        r1, r2 =  rpair
        r1aln = pst1(r1["sequence"][re_len:])
        r2aln = pst1(r2["sequence"][re_len:])
        if r1aln.query_begin != 0 or r1aln.query_end != 4 or \
                r2aln.query_begin != 0 or r2aln.query_end != 4:
            # Not read-through, but might be dimer.
            if not is_dimer(r1, r2, re_len, dimer_file):
                printrecs(r1, r2)
            else:
                n_adapt += 1
        if r1aln.target_begin == r2aln.target_begin:
            r1["sequence"] = r1["sequence"][:r1aln.target_begin]
            r1["quality"] = r1["quality"][:r1aln.target_begin]
            r2["sequence"] = "N"
            r2["quality"] = "#"
            n_trimmed += 1
            printrecs(r1, r2)
        if n_reads % 1000 == 0:
            print("Processed {:0.0f}K read pairs".format(n_reads /1000.0),
                  file=sys.stderr, end='\r')
    fh.close()
    if dimer_file is not None:
        dimer_file.close();
    print("Processed", n_reads, "read pairs", file=sys.stderr)
    print("Trimmed", n_trimmed, file=sys.stderr)
    print("Adaptor in", n_adapt, file=sys.stderr)

if __name__ == "__main__":
    opts = docopt.docopt(CLI_ARGS)
    main_ilfq(opts['<input_ILFQ>'], opts['-e'], opts['-d'])
