#!/usr/bin/env python3.5
# encoding: utf-8

# '''
# SugarPy - discovery-driven analysis of glycan compositions from IS-CID
# of intact glycopeptides.

# Copyright (c) 2016 S. Schulze, J. Kraegenbring, A. Oltmann, C. Fufezan,
# M. Hippler

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# '''

from __future__ import absolute_import
from itertools import combinations_with_replacement
from itertools import combinations
from collections import namedtuple
import ursgal
import sugarpy
import sys
import pymzml
import csv
import pyqms
import os
import pickle
import math
import re


# SPS_key = namedtuple(
#     'SPS_key',
#     []
# )

class Run(object):
    '''
    The SugarPy run class comprises the main functionality for glycopeptide matching.
    A typical workflow contain the following steps:
    1. parse_ident_file:
        Ursgal result files are parsed and peptide sequences, as well as their modifications 
        (except monosaccharides that would be part of the glycan) and retention times (RTs), are extracted.
    2. build_combinations and add_glycans2peptide:
        For a set of monosaccharides (param: monosaccharides) and maximal glycan length (param: max_tree_length), 
        all possible combinations of monosaccharides are calculated and the chemical compositions of the resulting 
        theoretical glycans (taking into account glycans with the same mass) are added to the chemical compositions
        of the extracted set of peptides.
    3. quantify:
        pyQms is used to build isotope envelope libraries for the theoretical glycopeptides 
        and to match them against all MS1 spectra within the given RT windows. 
        It should be noted that isotope envelopes consist of the theoretical m/z and relative intensity 
        for all isotopic peaks. Therefore, the quality of the resulting matches is indicated by an mScore, 
        which comprises the accuracies of the measured m/z and intensity.
    4. sort_results and validate_results:
        For each matched molecule, a score (VL) is calculated as the length of a vector for the mScore 
        (ranging from 0 to 1) and intensity (normalized by the maximum intensity of matched glycopeptides 
        within the run, therefore also ranging from 0 to 1).
        For each spectrum, all matched molecules are sorted by the glycan length (number of monosaccharides).
        Subsequently, starting with the longest glycan, for each glycan length, all glycan compositions
        are checked if they are part of any glycan composition of the previous level (longer glycan).
        Glycan compositions that are true subsets of larger, matched glycans (subtrees of those) 
        are considered fragment ions and are therefore merged with the larger, final glycopeptides.
        It should be noted that glycan compositions can be subtrees of multiple final glycans.
        Furthermore, fragmentation pathways are not taken into account, however, if Y1-ions are matched 
        (peptide harboring one monosaccharide), the corresponding monosaccharide is noted as the reducing end.
        For all final glycopeptides within one spectrum, the subtree coverage is calculated.
        Finally, the SugarPy score is calculated for each glycopeptide as the sum of vector lengths 
        from all corresponding subtrees (fragment ions Y0 to Yn) multiplied by the subtree coverage.
    '''

    def __init__(
        self,
        monosaccharides={}
    ):
        self.validated_results = []
        self.monosaccharides = monosaccharides
        # globals()['SPS_key'] = namedtuple(
        #     'SPS_key',
        #     sorted(monosaccharides.keys()) + ['End']
        # )

    def parse_ident_file(
        self,
        ident_file=None,
        unimod_glycans_incl_in_search=[],
    ):
        '''
        Parses an Ursgal results .csv file and extracts identified peptides
        together with their retention times.
        Glycans that were included in the search as modifications are removed.

        Keyword Arguments:
            ident_file (str): Path to the Ursgal result .csv file.
                This file should only include (potential) glycopeptides,
                i.e. it should be filtered.
            unimod_glycans_incl_in_search (list): List of Unimod PSI-MS names
                corresponding to glycans that were included in the database search
                as modification (will be removed from the peptide).

        Returns:
            dict: Lookup containing retention times and accuracies of all
                PSMs for each identified peptidoform (Peptide#Unimod:Pos)
        '''
        print('[ SugarPy  ] Parsing ident file ... ', ident_file)
        peptide_lookup = {}
        with open(ident_file, 'r') as in_file:
            csv_input = csv.DictReader(in_file)
            for line_dict in csv_input:
                peptide = line_dict['Sequence']
                rt = float(line_dict['Retention Time (s)']) / 60
                spec_id = int(line_dict['Spectrum ID'])

                # Organize unimods, remove monosaccharides from unimods, add
                # glycans to peptides
                mod_pattern = re.compile( r''':(?P<pos>[0-9]*$)''' )
                new_unimods = []
                for unimod in line_dict['Modifications'].split(';'):
                    match = mod_pattern.search(unimod)
                    if match is not None:
                        pos = int(match.group('pos'))
                        mod = unimod[:match.start()]
                        if mod in unimod_glycans_incl_in_search:
                            continue
                        new_unimods.append(unimod)
                if new_unimods != []:
                    peptide_unimod = '{0}#{1}'.format(
                        peptide, ';'.join(new_unimods))
                else:
                    peptide_unimod = peptide #+ '#'
                if peptide_unimod not in peptide_lookup.keys():
                    peptide_lookup[peptide_unimod] = {
                        'rt': set(),
                        'spec_id': set(),
                        'accuracy': [],
                        'protein': '',
                    }
                peptide_lookup[peptide_unimod]['rt'].add(float(rt))
                peptide_lookup[peptide_unimod]['spec_id'].add(spec_id)
                peptide_lookup[peptide_unimod]['accuracy'].append(
                    float(line_dict['Accuracy (ppm)'])
                )
                peptide_lookup[peptide_unimod][
                    'protein'] = line_dict['Protein ID']

        return peptide_lookup

    def add_glycans2peptide(
        self,
        peptide_list=[],
        max_tree_length=None,
        monosaccharides=None
    ):
        '''
        Adds chemical composition of glycans to a given list of peptides.
        Peptides need to be in unimod style (Peptide#Modifications).
        The chemical composition of the original peptidoform is returned as well.

        Keyword Arguments:
            peptide_list (list): List of peptides in unimod style
            max_tree_length (int): maximum number of monosaccharides in one combination
            monosaccharides(dict): dictionary containing name and chemical composition of monosaccharides

        Returns:
            dict: { 'Sequence#Modifications : {glycan_hill_notation': ['Name']}}
        '''
        if monosaccharides is None:
            monosaccharides == self.monosaccharides

        monosacch_combinations = self.build_combinations(
            max_tree_length=max_tree_length,
            monosaccharides=monosaccharides
        )
        peptides_with_glycans = {}
        for peptide in peptide_list:
            if peptide in peptides_with_glycans.keys():
                continue
            peptides_with_glycans[peptide] = {}
            # pep_with_glycans = {}
            cc = ursgal.ChemicalComposition()
            cc.use(peptide)
            hill_notation = cc.hill_notation_unimod()
            peptides_with_glycans[peptide][hill_notation] = [
                peptide]  # ['{0}'.format(peptide)]

            for composition in monosacch_combinations.keys():
                cc.add_chemical_formula(composition)
                hill_notation = cc.hill_notation_unimod()
                for combo in monosacch_combinations[composition]:
                    combo_name = ''
                    combo_dict = {}
                    for monosacch in combo:
                        if monosacch not in combo_dict.keys():
                            combo_dict[monosacch] = 0
                        combo_dict[monosacch] += 1
                    for monosacch in sorted(combo_dict.keys()):
                        combo_name += '{0}({1})'.format(
                            monosacch,
                            combo_dict[monosacch]
                        )
                    if hill_notation not in peptides_with_glycans[peptide].keys():
                        peptides_with_glycans[peptide][hill_notation] = []
                    peptides_with_glycans[peptide][hill_notation].append(
                        '{0}|{1}'.format(peptide, combo_name)
                    )
                cc.subtract_chemical_formula(composition)

        print('[ SugarPy  ] Added glycans to peptides.')
        return peptides_with_glycans

    def build_combinations(
        self,
        max_tree_length=None,
        monosaccharides=None,
        mode='replacement',
    ):
        '''
        Builds and returns a dictionary containing chemical compositions 
        of all combinations (with replacement, not ordered)
        of a given dict of monosaccharides and a maximal length of the tree.

        Keyword arguments:
            max_tree_length (int): Maximum number of monosaccharides in one combination
            monosaccharides(dict): Dictionary containing name and chemical composition of monosaccharides

        Returns:
            dict: keys: chemical compositions of all combinations (with replacement, not ordered), 
                values: combination(s) monosaccharide names corresponding to the chemical composition

        ToDo: change monosaccharides to list and get compositions from ursgal.ChemicalComposition(),
              keyword argument for calculate_formula?
        '''
        if monosaccharides is None:
            monosaccharides = self.monosaccharides

        if mode == 'replacement':
            print('[ SugarPy  ] Building combinations for:')
            print(
                '[ SugarPy  ]',
                len(monosaccharides),
                'given monosaccharides and a max tree length of',
                max_tree_length
            )
        # mode = 'sugarqb'
        # sugarqb_glycan_db = open('sugarqb_glycan_db.txt', 'w')
        glycan_combinations = {}
        for nr_repeats in range(1, max_tree_length + 1):
            if mode == 'combinations':
                glycan_combinations[nr_repeats] = []
                tmp_combinations = set()
                for combo in combinations(monosaccharides, nr_repeats):
                    tmp_combinations.add(combo)
                for tmp_combo in tmp_combinations:
                    glycan_dict = {}
                    for monosacch in set(tmp_combo):
                        count = tmp_combo.count(monosacch)
                        glycan_dict[monosacch] = count
                    glycan_combinations[nr_repeats].append(glycan_dict)
            elif mode == 'replacement':
                for combo in combinations_with_replacement(monosaccharides, nr_repeats):
                    cc = ursgal.ChemicalComposition()
                    for monosacch in combo:
                        cc.add_chemical_formula(monosaccharides[monosacch])
                    hill_notation = cc.hill_notation_unimod()
                    if hill_notation not in glycan_combinations:
                        glycan_combinations[hill_notation] = set()
                    glycan_combinations[hill_notation].add(combo)
            elif mode == 'sugarqb':
                for combo in combinations_with_replacement(monosaccharides, nr_repeats):
                    cc = ursgal.ChemicalComposition()
                    glycan_dict = {}
                    print(combo)
                    for monosacch in combo:
                        cc.add_glycan(monosacch)
                        count = combo.count(monosacch)
                        glycan_dict[monosacch] = count
                    sugarqb_list = []
                    for monosacch, count in glycan_dict.items():
                        sugarqb_list.append('{0}{1}'.format(monosacch, count))    
                    sugarqb_str = '{0}_N-Glycan {1}'.format(','.join(sorted(sugarqb_list)), round(cc._mass(), 6))
                    print(sugarqb_str, file=sugarqb_glycan_db)
            else:
                print('''
                    [ SugarPy  ] ERROR! the mode for build_combinations
                    [ SugarPy  ] is not available: {0}
                '''.format(mode))

        if mode == 'replacement':
            print('[ SugarPy  ] built', len(glycan_combinations), 'combinations')

        return glycan_combinations

    def build_rt_lookup(self, mzml_file, ms_level):
        '''
        Builds and returns a dictionary equivalent to the ursgal_lookup.pkl 
        It contains a dictionary (key is scan number, value is rt) for every
        mzml file name.

        Arguments:
            mzML_file: Path to the mzML file.
            ms_level: MS level for which the lookup should be built

        Returns:
            dict
        '''
        pymzml_run = pymzml.run.Reader(
            mzml_file,
            # extraAccessions = [
            #     ('MS:1000016', ['value', 'unitName'])
            # ],
            # obo_version = '1.1.0'
        )
        lookup = {}
        mzml_basename = os.path.basename(
            mzml_file
        ).replace('.mzML', '')
        lookup[mzml_basename] = {
            'scan_2_rt' : {}
        }
        for spectrum in pymzml_run:
            # scan_time, unit = spectrum.get('MS:1000016', (None,None))
            if spectrum.ms_level == ms_level:
                # if unit == 'second':
                #     rt_corr_factor = 60
                # else:
                    # rt_corr_factor = 1
                rt = spectrum.scan_time_in_minutes()  # / rt_corr_factor
                lookup[mzml_basename]['scan_2_rt'][spectrum.ID] = rt

        return lookup

    def quantify(
        self,
        molecule_name_dict=None,
        rt_window=None,
        ms_level=1,
        charges=None,
        params=None,
        pkl_name='',
        mzml_file=None,
        spectra=None,
        return_all=False,
        collect_precursor=False,
        force=False,
    ):
        '''
        Quantify a list of molecules in a given mzML file using pyQms.
        Quantification is done by default on MS1 level and can be specified
        for a retention time window.

        Keyword Arguments:
            molecule_name_dict (dict): contains for the molecules that should be quantified
                as hill notations (keys) a list of corresponding trivial names (values)
            rt_window (dict): optional argument to define a retention time window 
                in which the molecules are quantified (use 'min' and 'max' as keys in the dict)
            ms_level: MS level for which quantification should be performed
            charges (list): list of charge states that are quantified
            params (dict): pyQms parameters (see pyQms manual for further information)
            pkl_name (str): name of the result pickle containing the pyQms results
            mzml_file (str): path to the mzML file used for the quantification
            spectra (list): optional list of spectrum IDs that should be quantified
            return_all (bool): if True, in addition to the results pkl, the IsotopologueLibrary
                as well as the spectrum peaks are returned. This should only be used for 
                a single spectrum.

        Returns:
            str: path to the results pickle
        '''
        # quantify the shizznit
        print('[ SugarPy  ] Quantification for {0} molecules.'.format(
            len(molecule_name_dict)
        ))
        if os.path.exists(pkl_name) is False \
                or force or return_all or collect_precursor:

            molecules = []
            trivial_names = {}
            for formula in molecule_name_dict.keys():
                molecules.append('+{0}'.format(formula))
                trivial_names[
                    '+{0}'.format(formula)] = molecule_name_dict[formula]
            lib = pyqms.IsotopologueLibrary(
                molecules=molecules,
                charges=charges,
                metabolic_labels=None,
                fixed_labels=None,
                verbose=False,
                trivial_names=trivial_names,
                params=params
            )
            run = pymzml.run.Reader(
                mzml_file,
                # extraAccessions = [
                #     ('MS:1000016', ['value', 'unitName'])
                # ],
                # obo_version = '1.1.0'
            )
            results = None
            peaks = []
            precursor_to_rt_id = {}
            for n, spectrum in enumerate(run):
                if n % 100 == 0:
                    print(
                        '[ SugarPy  ] Processing spectrum number: {0}'.format(
                            n,
                        ),
                        end='\r'
                    )
                if collect_precursor and spectrum.ms_level >= 2:
                    rt = float(spectrum.scan_time_in_minutes())
                    selected_precursors = spectrum.selected_precursors
                    if selected_precursors is not None:
                        for precursor_dict in selected_precursors:
                            precursor_mz = precursor_dict['mz']
                            rounded_precursor_mz = round(precursor_mz, 3)
                            if rounded_precursor_mz not in precursor_to_rt_id.keys():
                                precursor_to_rt_id[rounded_precursor_mz] = []
                            precursor_to_rt_id[rounded_precursor_mz].append(
                                (rt, spectrum.ID))

                if spectrum.ms_level == ms_level:
                    if spectra != None:
                        if spectrum.ID not in spectra:
                            continue
                        # if spectrum['id'] != 7313: # 3136
                        #     continue
                    rt = float(spectrum.scan_time_in_minutes())
                    if rt_window != None:
                        rt_min, rt_max = rt_window
                        if rt < rt_min:
                            continue
                        elif rt > rt_max:
                            break
                    if return_all == True:
                        peaks.append(spectrum.peaks('centroided'))
                    results = lib.match_all(
                        mz_i_list=spectrum.peaks('centroided'),
                        file_name=mzml_file,
                        spec_id=spectrum.ID,
                        spec_rt=rt,
                        results=results
                    )

            results.lookup['formula to evidences'] = {}
            for molecule, formula in results.lookup['molecule to formula'].items():
                if formula not in results.lookup['formula to evidences'].keys():
                    results.lookup['formula to evidences'][formula] = {}
                if molecule not in results.lookup['formula to evidences'][formula].keys():
                    results.lookup['formula to evidences'][formula][molecule] = {
                        'trivial_names' : []
                    }
                if len(results.lookup['formula to trivial name'][formula]) >= 2:
                    print('this should never happen')
                results.lookup['formula to evidences'][ formula ][ molecule ]['trivial_names'] += \
                    results.lookup['formula to trivial name'][formula][0]

            pickle.dump(
                results,
                open(pkl_name, 'wb')
            )
        if return_all:
            # potential memory overkill
            results = pickle.load(
                open(
                    pkl_name,
                    'rb'
                )
            )
            return results, lib, peaks
        elif collect_precursor:
            return pkl_name, precursor_to_rt_id, lib
        else:
            return pkl_name

    def validate_results(
        self,
        pyqms_results_dict=None,
        min_spec_number=0,
        min_tree_length=0,
        monosaccharides=None
    ):
        '''
        Parse through pyQms results list and validate the results which includes the following:
        * sort_results: determines vector length for each matched spectrum with vectors beeing defined 
            by the mScore and the normalized intensity (normalized scaling factor). Also filters for
            a minimum number of spectra (for each molecule) in the results
        * build_match_dict: extracts information about glycan compositions, sorts glycans by their length
        * starting with the longest glycans, for each level (glycan length) the corresponding glycans
            are determined (glycans that are subtrees of longer glycans) and merged
        * the quality of glycan assignments is assessed by calculating the SugarPy_score ( (Sum of vector lengths)*subtree coverage),
            the subtree_coverage (Number of unique matched subtree lengths/Total number of unique subtree lengths) and the number
            of matched subtrees

        Keyword Arguments:
            pyqms_results_dict (dict): dictionary containing the Peptides#Unimod (key) and corresponding pyqms result pkl (value)
            min_spec_number (int): defines the minimum number of spectra for one matched formula.
            min_tree_length (int): minimum number of monosaccharides per glycan
            monosaccharides (dict): dictionary containing name and chemical composition of monosaccharides

        Returns
            results class object (dict): class (dict) containing all scored_glycans as well as the spec_collector 
                for every peptide_unimod
        '''
        if monosaccharides is None:
            monosaccharides = self.monosaccharides

        print(
            '[ SugarPy  ] Validating results, matching subtrees, calculating SugarPy_scores')
        # globals()['SPS_key'] = namedtuple(
        #     'SPS_key',
        #     sorted(monosaccharides.keys()) + ['End']
        # )
        sp_results = sugarpy.results.Results(
            monosaccharides=monosaccharides
        )
        n = 0
        for peptide, pkl_name in pyqms_results_dict.items():
            n += 1
            print('[ SugarPy  ] for peptide {0}: {1} out of {2}'.format(
                n, peptide, len(pyqms_results_dict)
            ))
            with open(pkl_name, 'rb') as open_pkl:
                results = pickle.load(open_pkl)
            spec_collector = self.sort_results(
                results=results, min_spec_number=min_spec_number
            )
            scored_glycans = {}
            for spec in spec_collector.keys():
                scored_glycans[spec] = {}
                final_glycans = {}
                match_dict, rt = self.build_match_dict(
                    spectrum_dict=spec_collector[spec])
                end_glycans2delete = set()
                for lvl, tree_length in enumerate(sorted(match_dict.keys(), reverse=True)):
                    for glyc in match_dict[tree_length]:
                        new_end = True
                        sp_glyc_dict = {}
                        for monosacch in sorted(monosaccharides.keys()):
                            sp_glyc_dict[monosacch] = glyc[
                                'glycan_comp'][monosacch]
                        if lvl != 0:
                            corresponding_trees = []
                            for end_glycan in final_glycans.keys():
                                end_glycan_dict = dict(end_glycan)
                                corresponding_tree = True

                                if tree_length != 0:  # not just peptide
                                    for mosa in sorted(glyc['glycan_comp'].keys()):
                                        if end_glycan_dict[mosa] < glyc['glycan_comp'][mosa]:
                                            corresponding_tree = False
                                else:  # just peptide
                                    final_glycans[end_glycan][
                                        'SugarPy_score'] += glyc['vector_length']
                                    final_glycans[end_glycan][
                                        'num_subtrees'] += 1
                                    final_glycans[end_glycan][
                                        'subtrees'].append(glyc['formula'])
                                    new_end = False
                                    continue

                                # either glycan (>1) or just peptide (0)
                                if corresponding_tree == True and tree_length != 1:
                                    final_glycans[end_glycan][
                                        'SugarPy_score'] += glyc['vector_length']
                                    final_glycans[end_glycan][
                                        'num_subtrees'] += 1
                                    final_glycans[end_glycan][
                                        'suc0r'].add(tree_length)
                                    final_glycans[end_glycan][
                                        'subtrees'].append(glyc['formula'])
                                    new_end = False

                                # make new entry/glycan for different endings of
                                # glycans
                                elif corresponding_tree == True and tree_length == 1:
                                    corresponding_trees.append(end_glycan)

                        # new_tree, either glycan (>1) or just peptide (0)
                        if new_end == True and tree_length != 1:
                            sp_glyc_dict['End'] = None
                            sp_tuple = tuple(sorted(sp_glyc_dict.items()))
                            final_glycans[sp_tuple] = {
                                'tree_length': tree_length,
                                'SugarPy_score': glyc['vector_length'],
                                'num_subtrees': 1,
                                'suc0r': set([tree_length]),
                                'formula': glyc['formula'],
                                'subtrees': [glyc['formula']],
                            }

                        # new_tree for different endings of glycans
                        elif new_end == True and tree_length == 1:
                            if lvl == 0:
                                sp_glyc_dict['End'] = list(glyc['glycan_comp'].keys())[
                                    list(glyc['glycan_comp'].values()).index(1)]
                                sp_tuple = tuple(sorted(sp_glyc_dict.items()))
                                final_glycans[sp_tuple] = {
                                    'tree_length': tree_length,
                                    'SugarPy_score': glyc['vector_length'],
                                    'num_subtrees': 1,
                                    'suc0r': set([tree_length]),
                                    'formula': glyc['formula'],
                                    'subtrees': [glyc['formula']],
                                }
                            else:
                                for end_glycan in corresponding_trees:
                                    end_glycan_dict = dict(end_glycan)
                                    # those have an end monosacch already
                                    if end_glycan_dict['End'] != None:
                                        continue
                                    ext_sp_dict = {}
                                    for monosacch, count in end_glycan:
                                        if monosacch == 'End':
                                            continue
                                        ext_sp_dict[monosacch] = count
                                    ext_sp_dict['End'] = list(glyc['glycan_comp'].keys())[
                                        list(glyc['glycan_comp'].values()).index(1)]
                                    ext_sp_tuple = tuple(
                                        sorted(ext_sp_dict.items())
                                    )
                                    if ext_sp_tuple not in final_glycans.keys():
                                        final_glycans[ext_sp_tuple] = {
                                            'tree_length': final_glycans[end_glycan]['tree_length'],
                                            'SugarPy_score': final_glycans[end_glycan]['SugarPy_score'],
                                            'num_subtrees': final_glycans[end_glycan]['num_subtrees'],
                                            'suc0r': final_glycans[end_glycan]['suc0r'] | {(tree_length)},
                                            'formula': final_glycans[end_glycan]['formula'],
                                            'subtrees': final_glycans[end_glycan]['subtrees'],
                                        }
                                    final_glycans[ext_sp_tuple][
                                        'SugarPy_score'] += glyc['vector_length']
                                    final_glycans[ext_sp_tuple][
                                        'num_subtrees'] += 1
                                    final_glycans[ext_sp_tuple][
                                        'subtrees'] += [glyc['formula']]
                                    end_glycans2delete.add(end_glycan)

                for end_glycan2delete in end_glycans2delete:
                    del final_glycans[end_glycan2delete]
                for final_glycan in final_glycans.keys():
                    if final_glycans[final_glycan]['tree_length'] < min_tree_length:
                        continue
                    # only peptide, no glycans matched
                    elif final_glycans[final_glycan]['tree_length'] == 0:
                        final_glycans[final_glycan]['sub_cov'] = 0
                    else:
                        final_glycans[final_glycan]['sub_cov'] = len(final_glycans[final_glycan]['suc0r']) /\
                            final_glycans[final_glycan]['tree_length']
                        final_glycans[final_glycan]['SugarPy_score'] = final_glycans[
                            final_glycan]['SugarPy_score'] * final_glycans[final_glycan]['sub_cov']
                    scored_glycans[spec][
                        final_glycan] = final_glycans[final_glycan]
            sp_results.add_results(
                peptide_unimod=peptide,
                scored_glycans=scored_glycans,
                spec_collector=spec_collector,
            )
        return sp_results

    def sort_results(self, results=None, min_spec_number=1):
        '''
        Parse through pyQms results, determines vector length for each matched spectrum 
        with vectors beeing defined by the mScore and the normalized intensity (normalized scaling factor).

        Keyword Arguments:
            min_spec_number (int): defines the minimum number of spectra for one matched formula.
            results (dict): pyQms results dictionary

        Returns
            dict: spec_collector = { matched_spectrum : { formula : {
                'vector' : [],
                'charge' : [],
                'trivial_name' : [],
                'glycan_comp' : [],
                'glycan_trees' : [],
            }
        '''
        # formula_collector = {}
        spec_collector = {}
        for m_key, results_value in results.items():
            ''' m_key contains named tuples: file_name, formula, charge, label_percentiles'''
            if results_value['len_data'] < min_spec_number:
                continue
            formula = m_key.formula
            trivial_name = results.lookup[
                'formula to trivial name'][formula][0]
            max_scaling_factor = 0
            vector_collector = {}
            for matched_spectrum in results_value['data']:
                if matched_spectrum.spec_id not in vector_collector.keys():
                    vector_collector[matched_spectrum.spec_id] = {
                        'vector_length' : 0,
                        'scaling_factor' : 0,
                        'mscore' : 0,
                        'rt' : 0,
                        'peaks' : [],
                    }
                vector_collector[matched_spectrum.spec_id][
                    'scaling_factor'] = matched_spectrum.scaling_factor
                if max_scaling_factor < matched_spectrum.scaling_factor:
                    max_scaling_factor = matched_spectrum.scaling_factor
                vector_collector[matched_spectrum.spec_id][
                    'mscore'] = matched_spectrum.score
                vector_collector[matched_spectrum.spec_id][
                    'rt'] = matched_spectrum.rt
                vector_collector[matched_spectrum.spec_id][
                    'peaks'] = matched_spectrum.peaks
            for matched_spectrum in vector_collector.keys():
                norm_scaling_factor = vector_collector[matched_spectrum][
                    'scaling_factor'] / max_scaling_factor
                vector_collector[matched_spectrum]['vector_length'] =  \
                    math.sqrt((norm_scaling_factor**2) +
                              (vector_collector[matched_spectrum]['mscore']**2))
                if matched_spectrum not in spec_collector.keys():
                    spec_collector[matched_spectrum] = {}
                if formula not in spec_collector[matched_spectrum].keys():
                    spec_collector[matched_spectrum][formula] = {
                        'vector' : [],
                        'charge' : [],
                        'trivial_name' : [],
                        'glycan_comp' : [],
                        'glycan_trees' : [],
                    }
                spec_collector[matched_spectrum][formula][
                    'vector'].append(vector_collector[matched_spectrum])
                spec_collector[matched_spectrum][formula][
                    'charge'].append(m_key.charge)
                spec_collector[matched_spectrum][formula][
                    'trivial_name'].append(trivial_name)

        return spec_collector

    def build_match_dict(self, spectrum_dict=None):
        '''
        Uses a spec_collector (see sort_results) spectrum to extract the peptide
        and glycan composition and to sort glycans according to their length.

        Keyword Arguments:
            spectrum_dict (dict): dictionary for a spectrum from spec_collector
                (see sort_results)

        Returns:
            dict: { n (tree length) : [{
                                        formula: ,
                                        glycan_comp: ,
                                        vector_length: ,
                                        }, ... ]
        }
        '''
        match_dict = {}
        for formula in spectrum_dict.keys():
            for i, name_list in enumerate(spectrum_dict[formula]['trivial_name']):
                for name in name_list:
                    pep, glycan_comp = self.extract_pep_and_glycan_comp(
                        trivial_name=name
                    )
                    tree_length = 0
                    for monosacch in glycan_comp.keys():
                        tree_length += glycan_comp[monosacch]
                    if tree_length not in match_dict.keys():
                        match_dict[tree_length] = []
                    match_dict[tree_length].append({
                        'formula' : formula,
                        'glycan_comp' : glycan_comp,
                        'vector_length' : spectrum_dict[formula]['vector'][i]['vector_length'],
                        'charge' : spectrum_dict[formula]['charge'][i]
                    })
        rt = spectrum_dict[formula]['vector'][0]['rt']
        return match_dict, rt

    def extract_pep_and_glycan_comp(self, trivial_name=None):
        '''
        Use the trivial name to extract information about the peptide and glycan composition.
        This also works for just extracting the glycan composition from a glycan string.

        Keyword Argumens:
            trivial_name('str'): trivial name of the glycan ('HexNAc(2)Hex(5)') or 
                glycopeptide ('PEPTIDE|HexNAc(2)Hex(5)')

        Returns:
            peptide(str): peptide sequence
            glycan_comp(dict): glycan composition as dictionary with monosaccharides as keys
                and their number as value
        '''
        glycan_comp = {}
        for monosacch in self.monosaccharides.keys():
            glycan_comp[monosacch] = 0
        if '|' in trivial_name:
            pep, glycan = trivial_name.split('|')
        elif '(' in trivial_name:
            pep = '',
            glycan = trivial_name
        else:
            pep = trivial_name
            glycan = ''
        pattern = re.compile(r'(?P<element>[A-z,0-9]*)[(](?P<count>[0-9]*)[)]')
        for match in pattern.finditer(glycan):
            if match.group('count') == '':
                count = 1
            else:
                count = int(match.group('count'))
            glycan_comp[match.group('element')] += count

        return pep, glycan_comp

if __name__ == '__main__':
    print(__doc__)
