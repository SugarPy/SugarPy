#!/usr/bin/env python3
# encoding: utf-8

import sugarpy
import sys
import pprint

def main(input_path=None):
    '''
    This script generates a peptide lookup that can be used as an input
    for further SugarPy functions. Here, an Ursgal result file is parsed
    and pepide sequences as well as corresponding informations are stored
    in a dictionary. However, This dictionary can also be generated in 
    different ways, but should follow this structure::

        peptide_lookup = {
            'sequence#unimod:pos': {    
                'rt': set(),
                'spec_id': set(),
                'accuracy': [],
                'protein': 'description',
            }
        }

    Usage::

        ./parse_peptide_sequence_input_file.py <inpu_result_file.csv>
    '''
    sp_run = sugarpy.run.Run()

    peptide_lookup = sp_run.parse_ident_file(
        ident_file=input_path,
        unimod_glycans_incl_in_search=[
            'HexNAc',
            'HexNAc(2)',
        ],
    )
    pprint.pprint(peptide_lookup)

if __name__ == '__main__':
    main(input_path=sys.argv[1])