#!/usr/bin/env python
# encoding: utf-8

import sugarpy
import ursgal
import sys
import os
import pickle

def main(
    mzml_file=None,
    validated_results_pkl=None,
    ursgal_result_file=None,
    ms_level=1,
    output_file=None,
    min_sugarpy_score=1,
    min_sub_cov=0.5,
    min_spec_number=1,
    charges=[1, 2, 3, 4, 5],
    monosaccharides={
        "dHex": 'C6H10O4',
        "Hex": 'C6H10O5',
        "HexNAc": 'C8H13NO5',
        "NeuAc": 'C11H17NO8',
    },
    rt_border_tolerance=1,
    max_tree_length=10,
    unimod_glycans_incl_in_search=[
        'HexNAc',
        'HexNAc(2)',
    ],
    max_trees_per_spec=5,
):
    '''
    This example script takes the following three files as input:

    * Ursgal result file: this file contains the identified glycopeptide sequences in .csv format
    * mzML file: the raw MS data in .mzML format
    * SugarPy results pkl: pkl file that contains SugarPy results (see sugarpy_run_exampe.py on how to generate it)

    From those files, a SugarPy results csv is generated, as well as a glycopeptide elution profile.

    Usage::

        ./sugarpy_results_example.py <ursgal_result_file.csv> <ms_data_file.mzML> <sugarpy_results.pkl>
    '''

    # load results from pkl
    with open(validated_results_pkl, 'rb') as results_pkl:
        validated_results = pickle.load(results_pkl)

    mzml_basename = os.path.basename(mzml_file)
    sp_run = sugarpy.run.Run()
    lookup_dict = sp_run.build_rt_lookup(mzml_file, ms_level)
    scan2rt = lookup_dict[mzml_basename]['scan_2_rt']

    # load results class
    sp_results = sugarpy.results.Results(
        monosaccharides=monosaccharides,
        scan_rt_lookup=scan2rt,
        validated_results=validated_results,
    )

    # write results csv
    sp_run = sugarpy.run.Run()
    peptide_lookup = sp_run.parse_ident_file(
        ident_file=ursgal_ident_file,
        unimod_glycans_incl_in_search=unimod_glycans_incl_in_search
    )
    sp_result_file = sp_results.write_results2csv(
        output_file=validated_results_pkl.replace('.pkl', '.csv'),
        max_trees_per_spec=max_trees_per_spec,
        min_sugarpy_score=min_sugarpy_score,
        min_sub_cov=min_sub_cov,
        peptide_lookup=peptide_lookup,
        monosaccharides=monosaccharides,
        scan_rt_lookup=scan2rt,
        mzml_basename=mzml_basename
    )


    # plot glycopeptide elution profile
    plot_molecule_dict = sp_results.parse_result_file(sp_result_file)
    elution_profile_file = sp_results.plot_glycan_elution_profile(
        peptide_list=sorted(plot_molecule_dict.keys()),
        min_sugarpy_score=min_sugarpy_score,
        min_sub_cov=min_sub_cov,
        x_axis_type='retention_time',
        score_type='top_scores',
        output_file=output_file.replace('.txt', '_glycan_profile.html'),
        plotly_layout={
            'font' : {
                'family':'Arial',
                'size': 26,
                'color' : 'rgb(0, 0, 0)',
            },
            'width':1700,
            'height':1200,
            'margin':{
                'l':120,
                'r':50,
                't':50,
                'b':170,
            },
            'xaxis':{
                'type':'linear',
                'color':'rgb(0, 0, 0)',
                'title': 'Retention time [min]',
                'title_font':{
                    'family':'Arial',
                    'size':26,
                    'color':'rgb(0, 0, 0)',
                },
                'autorange':True,
                'fixedrange':False,
                'tickmode':'auto',
                'showticklabels':True,
                'tickfont':{
                    'family':'Arial',
                    'size':22,
                    'color':'rgb(0, 0, 0)',
                },
                'tickangle':0,
                'ticklen':5,
                'tickwidth':1,
                'tickcolor':'rgb(0, 0, 0)',
                'ticks':'outside',
                'showline':True,
                'linecolor':'rgb(0, 0, 0)',
                'mirror':False,
                'showgrid':False,
                'anchor':'y',
                'side':'bottom',
            },
            'yaxis':{
                'type':'linear',
                'color':'rgb(0, 0, 0)',
                'title':'Max. SugarPy score',
                'title_font':{
                    'family':'Arial',
                    'size':26,
                    'color':'rgb(0, 0, 0)',
                },
                'autorange':True,
                'fixedrange':False,
                'tickmode':'auto',
                'showticklabels':True,
                'tickfont':{
                    'family':'Arial',
                    'size':22,
                    'color':'rgb(0, 0, 0)',
                },
                'tickangle':0,
                'ticklen':5,
                'tickwidth':1,
                'tickcolor':'rgb(0, 0, 0)',
                'ticks':'outside',
                'showline':True,
                'linecolor':'rgb(0, 0, 0)',
                'mirror':False,
                'showgrid':False,
                'zeroline':False,
                'anchor':'x',
                'side':'left',
            },
            'legend':{
                'font':{
                    'family':'Arial',
                    'size':5,
                    'color':'rgb(0, 0, 0)',
                },
                'orientation':'v',
                'traceorder':'normal',
            },
            'showlegend':True,
            'paper_bgcolor':'rgba(0,0,0,0)',
            'plot_bgcolor':'rgba(0,0,0,0)',
        },
    )

    print('Done.')
    return output_file


if __name__ == '__main__':
    ursgal_ident_file = sys.argv[1]
    mzml_file = sys.argv[2]
    validated_results_pkl = sys.argv[3]
    output_file = result_file.replace('.csv', '.html')
    main(
        mzml_file=mzml_file,
        validated_results_pkl=validated_results_pkl,
        ursgal_result_file=ursgal_result_file,
        output_file=output_file,
    )