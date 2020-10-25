.. code-block:: python

	#!/usr/bin/env python
	# encoding: utf-8
	
	import sugarpy
	import ursgal
	import sys
	import os
	import statistics
	import pickle
	from collections import namedtuple
	
	
	def main(
	    ident_file=None,
	    unimod_glycans_incl_in_search=[
	        'HexNAc',
	        'HexNAc(2)',
	    ],
	    max_tree_length=10,
	    monosaccharides={
	        "dHex": 'C6H10O4',
	        "Hex": 'C6H10O5',
	        "HexNAc": 'C8H13NO5',
	        "NeuAc": 'C11H17NO8',
	    },
	    mzml_file=None,
	    scan_rt_lookup=None,
	    rt_border_tolerance=1,
	    charges=[1, 2],
	    output_file=None,
	    pyqms_params={},
	    min_spec_number=1,
	    min_tree_length=1,
	    max_trees_per_spec=5,
	    min_sugarpy_score=0,
	    min_sub_cov=0.5,
	    ms_level=1,
	    force=True,
	):
	    '''
	    This example script takes the following two files as input and
	    than performs a search for intact glycopeptides:
	
	    * Ursgal result file: this file contains the identified glycopeptide sequences in .csv format
	    * mzML file: the raw MS data in .mzML format
	
	    Usage::
	
	        ./sugarpy_run_complete.py <ursgal_result_file.csv> <ms_data_file.mzML>
	    '''
	
	    # initialize Sugary run
	    sp_run = sugarpy.run.Run(
	        monosaccharides=monosaccharides
	    )
	
	    # read the Ursgal result file and store pepides in lookup
	    peptide_lookup = sp_run.parse_ident_file(
	        ident_file=ident_file,
	        unimod_glycans_incl_in_search=unimod_glycans_incl_in_search,
	    )
	
	    # build all theoretical glycopeptide combinations
	    pep_unimod_list = list(peptide_lookup.keys())
	    peps_with_glycans = sp_run.add_glycans2peptide(
	        peptide_list=pep_unimod_list,
	        max_tree_length=max_tree_length,
	    )
	
	    lookup_dict = sp_run.build_rt_lookup(mzml_file, ms_level)
	    mzml_basename = os.path.basename(mzml_file).replace('.mzML', '')
	    scan2rt = lookup_dict[mzml_basename]['scan_2_rt']
	    output_folder = os.path.dirname(output_file)
	
	    # run pyqms to identify isotope envelopes for all theoretical glycopeptides
	    pyqms_results = {}
	    for n, peptide_unimod in enumerate(sorted(peps_with_glycans)):
	        pkl_name = '-'.join(peptide_unimod.split(';'))
	        pkl_name = '_'.join(pkl_name.split(':'))
	        rt_min = min(peptide_lookup[peptide_unimod]
	                     ['rt']) - rt_border_tolerance
	        rt_max = max(peptide_lookup[peptide_unimod]
	                     ['rt']) + rt_border_tolerance
	        print(
	            '[ SugarPy  ] Quantification for peptide ',
	            peptide_unimod,
	            '#{0} out of {1}'.format(
	                n + 1,
	                len(peps_with_glycans)
	            )
	        )
	        results_pkl = sp_run.quantify(
	            molecule_name_dict=peps_with_glycans[peptide_unimod],
	            rt_window=(rt_min, rt_max),
	            ms_level=ms_level,
	            charges=charges,
	            params=pyqms_params,
	            pkl_name=os.path.join(
	                output_folder,
	                '{0}_{1}_pyQms_results.pkl'.format(
	                    mzml_basename,
	                    pkl_name
	                )
	            ),
	            mzml_file=mzml_file,
	        )
	        pyqms_results[peptide_unimod] = results_pkl
	
	    # combine matched glycopeptide fragments to inteact glycopeptides
	    validated_results = sp_run.validate_results(
	        pyqms_results_dict=pyqms_results,
	        min_spec_number=min_spec_number,
	        min_tree_length=min_tree_length,
	    )
	
	    # initiate results class and save Sugary results as csv file
	    sp_results = sugarpy.results.Results(
	        validated_results=validated_results,
	        monosaccharides=monosaccharides
	    )
	    sp_results.write_results2csv(
	        output_file=output_file,
	        max_trees_per_spec=max_trees_per_spec,
	        min_sugarpy_score=min_sugarpy_score,
	        min_sub_cov=min_sub_cov,
	        peptide_lookup=peptide_lookup,
	        scan_rt_lookup=scan2rt,
	        mzml_basename=mzml_basename,
	    )
	    output_pkl = output_file.replace('.csv', '.pkl')
	    with open(output_pkl, 'wb') as out_pkl:
	        pickle.dump(sp_results, out_pkl)
	    print('Done.')
	
	if __name__ == '__main__':
	    if len(sys.argv) < 3:
	        print(main.__doc__)
	        sys.exit(1)
	    mzml_file = sys.argv[1]
	    ident_file = sys.argv[2]
	    output_file = ident_file.replace('.csv', '_sugarpy.csv')
	    main(
	        ident_file=ident_file,
	        mzml_file=mzml_file,
	        output_file=output_file,
	    )
