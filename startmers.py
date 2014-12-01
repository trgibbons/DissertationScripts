#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 12:46:23 2014

@author: tgibbons
"""

import sys
import argparse
from Bio import SeqIO


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

    w = len(str(args.max))  # width of smer-length field in output file names

    for k in range(args.min, args.max+1):

        if args.split:
            wmers, cmers, smers = count_startmers(args.fasta, s=k, split=True)
        else:
            smers = count_startmers(args.fasta, s=k, split=False)

        args.fasta.seek(0)
    
        pref = "{0}_s{1}".format(args.opref, str(k).zfill(w))
        print_startmers(smers, out=pref+".smers", lex=args.lex)
    
        if args.split:
            print_startmers(wmers, out=pref+".wmers", lex=args.lex)
            print_startmers(cmers, out=pref+".cmers", lex=args.lex)


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
        description="Count 'startmers' in a multi-FASTA file")

    parser.add_argument('-m', '--min',
                        dest='min',
                        type=int,
                        action='store',
                        default=4,
                        help="Minimum startmer length [def=4]")

    parser.add_argument('-M', '--max',
                        dest='max',
                        type=int,
                        action='store',
                        default=22,
                        help="Maximum startmer length [def=22]")

    parser.add_argument('--lex',
                        dest='lex',
                        action='store_true',
                        default=False,
                        help="Sort output lexicographicaly instead of by " +
                             "abundance [def=False]")

    parser.add_argument('--split',
                        dest='split',
                        action='store_true',
                        default=False,
                        help="Print individual output files corresponding " +
                             "to the Watson and Crick strands, in addition " +
                             "to the default output file that combines both " +
                             "[def=False]")

    parser.add_argument('fasta',
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="Multi-FASTA format input file")

    parser.add_argument('opref',
                        action='store',
                        help="Prefix for the output files: w(atson)mers, " +
                             "c(rick)mers, and combined s(tart)mers")

    args = parser.parse_args()

    return args


def count_startmers(fasta, s, split=False):
    """Count startmers in a multi-FASTA file
    """
    smers = dict()  # Combined startmer counts
    if split:
        wmers = dict()  # Watson-mers - kmers at the 5' end
        cmers = dict()  # Crick-mers - reverse-complemented kmers at the 3' end

    for rec in SeqIO.parse(fasta, "fasta"):
        wmer = str(rec.seq[:s])
        cmer = str(rec.seq.reverse_complement()[:s])

        increment_smer(smers, wmer)
        increment_smer(smers, cmer)

        if split:
            increment_smer(wmers, wmer)
            increment_smer(cmers, cmer)

    if split:
        return wmers, cmers, smers
    else:
        return smers


def increment_smer(smer_dict, smer):
    try:
        smer_dict[smer] += 1
    except KeyError:
        smer_dict[smer] = 1


def print_startmers(smers, out, lex=False):
    handle = open(out, 'w')
    if lex:
        for smer, count in sorted(smers.iteritems()):
            handle.write("{0}\t{1}\n".format(smer, count))
    else:
        for smer, count in sorted(smers.items(),
                                  key=lambda x: x[1],
                                  reverse=True):
            handle.write("{0}\t{1}\n".format(smer, count))
    handle.close()


if __name__ == "__main__":
    sys.exit(main())
