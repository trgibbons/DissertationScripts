#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 18:17:11 2014

@author: tgibbons
"""

#TODO: Split this into 3 separate programs:
#      -polyA.py
#      -internalDinoSL.py
#      -serialDinoSL.py

import sys
import argparse
import re
import os.path as op
from bisect import bisect
from glob import glob

from Bio import SeqIO
from fuzzywuzzy import fuzz

def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.

    I typically try to be pretty good about following PEP-8 guidelines,
    wrapping my code and what not, and not hard-coding information into a
    program without having a runtime option to override it. I work on a pretty
    large screens though, and this program contains many long strings that I
    just found to be much more readable when I let them stretch out.

    The hard-coded information in this program is current as of Aug 7, 2014.

    Args:
        None

    Returns:
        An exit status (hopefully 0)
    """
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()

    if args.resume:
        import shutil
        shutil.move("polyadenylated_seqs.txt", "polyadenylated_seqs.old")
        shutil.move("trans-spliced_seqs.txt", "trans-spliced_seqs.old")
        shutil.move("serially_trans-spliced_seqs.txt",
                    "serially_trans-spliced_seqs.old")
        out = open("mature_transcript_evidence_summary_stats.txt", 'w')
        paout = open("polyadenylated_seqs.txt", 'w')
        dslout = open("trans-spliced_seqs.txt", 'w')
        sdslout = open("serially_trans-spliced_seqs.txt", 'w')

        seq_stats, final_seqs = initialize_vars()
        read_in_paout(paout, seq_stats, final_seqs)

    else:
        out = open("mature_transcript_evidence_summary_stats.txt", 'w')
        paout = open("polyadenylated_seqs.txt", 'w')
        dslout = open("trans-spliced_seqs.txt", 'w')
        sdslout = open("serially_trans-spliced_seqs.txt", 'w')

        paout.write("Assembly\tSeqID\tOrientation\tTailLength\tSequence\n")
        dslout.write("Assembly\tSeqID\tOrientation\tPosition\tSequence\n")
        sdslout.write("Assembly\tSeqID\tOrientation\tPosition\t" + \
                      "CopyNumber\tPercID\tSubSeq\tSequence\n")

        polyA_seqs = set()
        dsl_seqs = set()
        sdsl_seqs = set()

    dinosl = "CCGTAGCCATTTTGGCTCAAG"
    slsuf = "CCATTTTGGCTCAAG"

    ass_rgx = "/Volumes/Research/2014/DAToL/dinosl/A*.fasta"
    for ass_file in glob(ass_rgx):
        assembly = op.basename(ass_file)[:-6]

        out.write("{0}\n".format(assembly))

        for rec in SeqIO.parse(ass_file, "fasta"):
            fws = str(rec.seq)  # ForWard Sequence
            rcs = str(rec.seq.reverse_complement())  # RC Sequence

            polya_tail(assembly, rec.id, polyA_seqs, fws, rcs, paout)

            sliding_window(assembly, rec.id, fws, rcs,
                           dinosl, dsl_seqs, dslout,
                           slsuf, sdsl_seqs, sdslout)

        out.write("PolyA:\t{0}\n".format(len(polyA_seqs)))
        out.write("ExactDinoSL:\t{0}\n".format(len(dsl_seqs)))
        out.write("FuzzySerialDinoSL:\t{0}\n\n".format(len(sdsl_seqs)))
        polyA_seqs = set()
        dsl_seqs = set()
        sdsl_seqs = set()

    out.close()
    paout.close()
    dslout.close()
    sdslout.close()


def get_parsed_args():
    """Parse the command line arguments

    Parses command line arguments using the argparse package, which is a
    standard Python module starting with version 2.7.

    Args:
        None, argparse fetches them from user input

    Returns:
        args: An argparse.Namespace object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="This program aggregates statistics about the relative " +
            "abundance of DinoSL suffixes within a set of precomputed " +
            "startmer distributions. Everything is hard coded. There are no " +
            "arguments or options. The output is printed to a file called " +
            "'exact_dinosl_suffix_contig_startmer_abundances.Rtab'.")

    parser.add_argument('-r', '--resume', dest='resume',
                        default=False, action='store_true',
                        help="Resume previous run [def=False]")

    args = parser.parse_args()

    return args


def initialize_vars():
    ass_rgx = "/Volumes/Research/2014/DAToL/dinosl/P*.fasta"
    for ass_file in glob(ass_rgx):
        assembly = op.basename(ass_file)[:-6]


def read_in_paout(paout):
    polyA_seqs = set()
    last_assembly = ''
    last_seq = ''

    for line in open("polyadenylated_seqs.old"):
        temp = line.strip().split()
        if not temp:
            continue
        elif temp[0][0] == '#':
            continue
        elif len(temp) != 5:
            break

        if not temp[0] == last_assembly:
            for entry in polyA_seqs:
                paout.write('\t'.join(entry)+'\n')

        last_assembly = temp[0]
        last_seq = temp[1]

        polyA_seqs.append(temp)

    return polyA_seqs, last_assembly, last_seq


def polya_tail(assembly, sid, polyA_seqs, fws, rcs, paout):
    fw_tail = re.search('A{5,}$', fws)
    rc_tail = re.search('A{5,}$', rcs)

    if fw_tail:
        tail_len = str(len(fw_tail.group()))
        paout.write('\t'.join([assembly, sid, "+", tail_len, fws])+'\n')

    if rc_tail:
        tail_len = str(len(rc_tail.group()))
        paout.write('\t'.join([assembly, sid, "-", tail_len, rcs])+'\n')


def sliding_window(assembly, sid, fws, rcs,
                   dinosl, dsl_seqs, dslout,
                   slsuf, sdsl_seqs, sdslout):
    """Sliding window subsequence comparisons, both exact and fuzzy"""
    fw_stack = False
    rc_stack = False

    # Indexing past the end of a string doesn't throw an error
    for i in range(len(fws)):
        exact_dinosl(assembly, sid, i,
                     dinosl, dsl_seqs, fws, rcs, dslout)

        for j in range(5, 0, -1):
            sub_seq = fws[i:i+15*j]
            sub_rcs = rcs[i:i+15*j]
            fw_ratio = fuzz.ratio(sub_seq, slsuf*j)
            rc_ratio = fuzz.ratio(sub_rcs, slsuf*j)
        
            if fw_ratio > 75:
                sdsl_seqs.add(sid)
                fw_list = [assembly, sid, '+', i+1, j, fw_ratio, sub_seq, fws]
                if not fw_stack:
                    fw_stack = [fw_list]
                elif any([(x[3]-1)+15*x[4] >= i for x in fw_stack]):
                    fw_stack.append(fw_list)
                else:
                    process_stack(fw_stack, sdslout)
                    fw_stack = [fw_list]
        
            if rc_ratio > 75:
                sdsl_seqs.add(sid)
                rc_list = [assembly, sid, '-', i+1, j, rc_ratio, sub_rcs, rcs]
                if not rc_stack:
                    rc_stack = [rc_list]
                elif any([(x[3]-1)+15*x[4] >= i for x in rc_stack]):
                    rc_stack.append(rc_list)
                else:
                    process_stack(rc_stack, sdslout)
                    rc_stack = [rc_list]

    # EOF: process final entries
    if fw_stack:
        process_stack(fw_stack, sdslout)
    if rc_stack:
        process_stack(rc_stack, sdslout)


def exact_dinosl(assembly, sid, i, dinosl, dsl_seqs, fws, rcs, dslout):
    """Print exact matches to the canonical 21bp DinoSL sequence"""
    if fws[i:i+21] == dinosl:
        dsl_seqs.add(sid)
        fw_list = [str(x) for x in [assembly, sid, '+', i+1, fws]]
        dslout.write('\t'.join(fw_list)+'\n')
    if rcs[i:i+21] == dinosl:
        dsl_seqs.add(sid)
        rc_list = [str(x) for x in [assembly, sid, '-', i+1, rcs]]
        dslout.write('\t'.join(rc_list)+'\n')


def process_stack(stack, sdslout):
    """Avoid printing duplicate matches
    
    An exact match will allow up to 25% of the sequence to be missing from
    either end, creating many hits for each instance of highly-conserved dinosl
    sequences. This fuction takes a stack of adjacent matches, discards hits
    that do not contain the largest number of serial-SL multiples observed in
    the stack, then discards remaining hits whos sequence identities do not
    match the maximum observed score. All remaining hits are then printed.
    """
    better_stack = []
    best_stack = []

    max_rep = max([x[4] for x in stack])
    for i in range(len(stack)):
        if stack[i][4] == max_rep:
            better_stack.append(stack[i])

    max_pid = max([x[5] for x in better_stack])
    for i in range(len(better_stack)):
        if better_stack[i][5] == max_pid:
            best_stack.append([str(x) for x in better_stack[i]])

    for i in range(len(best_stack)):
        sdslout.write('\t'.join(best_stack[i])+'\n')


if __name__ == "__main__":
    sys.exit(main())
