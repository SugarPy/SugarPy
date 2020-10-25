#!/usr/bin/env python3
# encoding: utf-8

import sugarpy
import pprint

def main():
    '''
    For a given list of peptides and a set of monosaccharides,
    this script builds all theoretical combinations of glycopeptides.

    Usage::
    
        ./build_glycopeptide_combination.py
    '''
    sp_run = sugarpy.run.Run(
        monosaccharides={
            "dHex": 'C6H10O4',
            "Hex": 'C6H10O5',
            "HexNAc": 'C8H13NO5',
            "NeuAc": 'C11H17NO8',
        },
    )

    pep_unimod_list = [
        'GLYCANSTPEPTIDE',
    ]

    glycopeptide_combinations = sp_run.add_glycans2peptide(
        peptide_list=pep_unimod_list,
        max_tree_length=5,
    )

    pprint.pprint(glycopeptide_combinations)

if __name__ == '__main__':
    main()