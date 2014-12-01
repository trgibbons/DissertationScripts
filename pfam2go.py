#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 17:43:13 2014

@author: tgibbons
"""

import sys


def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.

    Args:
        None

    Returns:
        An exit status (hopefully 0)
    """
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()
    pf2go = import_pfam2go(args.pfam2go)
    convert_pfam_dat(args.dat, pf2go, args.out, args.NA)


def get_parsed_args():
    """Parse the command line arguments

    Parses command line arguments using the argparse package, which is a
    standard Python module starting with version 2.7.

    Args:
        None, argparse fetches them from user input

    Returns:
        args: An argparse.Namespace object containing the parsed arguments
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Convert a HMMer domtblout '.dat' file into a "
                    "tab-delimited file containing sequence IDs, Pfam IDs " +
                    "and descriptions, and GO IDs and descriptions.")

    parser.add_argument('--NA',
                        dest='NA',
                        default=False,
                        action='store_true',
                        help="Print lines corresponding to Pfam domains for" +
                             "which pfam2go offers no GO mapping, listing " +
                             "instead 'NA' for both the GO ID and " +
                             "GO description [def=False]")

    parser.add_argument('pfam2go',
                        type=argparse.FileType('r'),
                        help="pfam2go file published by the Gene Ontology " +
                             "Consortium")

    parser.add_argument('dat',
                        type=argparse.FileType('r'),
                        help="Tab-delimited Pfam '.dat' file")

    parser.add_argument('out',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="Output file name [def=stdout]")

    args = parser.parse_args()

    return args


def import_pfam2go(pfam2go_handle):
    pf2go = dict()
    for line in pfam2go_handle:
        if line[0] == '!':
            continue
        pf, go = line.strip().split(' > ')
        pf_id, pf_name = parse_pfam2go_pf_line(pf)
        go_id, go_desc = parse_pfam2go_go_line(go)
        try:
            pf2go[pf_id]['go'][go_id] = go_desc
        except KeyError:
            pf2go[pf_id] = dict(name=pf_name,
                                go={go_id:go_desc})

    return pf2go


def parse_pfam2go_pf_line(pf):
    pf = pf[5:]
    pf_id, pf_desc = pf.split()
    return pf_id, pf_desc


def parse_pfam2go_go_line(go):
    go_desc, go_id = go.split(' ; ')
    go_desc = go_desc[3:]
    return go_id, go_desc


def convert_pfam_dat(pfam_dat_handle, pf2go, out_handle, NA):
    """Reduce a HMMer domtblout file to a few columns + GO annotations
    
    The function iterates through a HMMer 3.0 domtblout '.dat' file generated
    by annotating protein sequences with hmmscan and Pfam. Pfam domains not
    listed in the pfam2go file
    """
    for line in pfam_dat_handle:
        temp = line.strip().split()
        seq_id = str(temp[3])
        pf_id = str(temp[1]).split('.')[0]
        try:
            pf_name = pf2go[pf_id]['name']
            for go_id, go_desc in pf2go[pf_id]['go'].iteritems():
                out_handle.write("\t".join([seq_id, pf_id, pf_name, 
                                            go_id, go_desc])+"\n")
        except KeyError:
            if NA:
#                pf_desc = ' '.join(temp[22:])
                pf_name = temp[0]
                out_handle.write("\t".join([seq_id, pf_id, pf_name,
                                            "NA", "NA"])+"\n")


if __name__ == "__main__":
    sys.exit(main())
