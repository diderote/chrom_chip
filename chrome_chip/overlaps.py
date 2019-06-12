#!/usr/bin/env python3

import os
from datetime import datetime

import pandas as pd
import rpy2.robjects as ro
import rpy2.rinterface as ri
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import gseapy

from chrome_chip.plot import plot_venn2, plot_venn2_set, plot_venn3_counts, plot_venn3_set
from chrome_chip.common import out_result, rout_write, output, make_folder, load_bedtool, read_pd, bed2df


def overlaps(exp):
    '''
    Performs overlaps of two or more de_sig lists.
    '''

    out_dir = make_folder(f'{exp.scratch}/Overlaps/')

    for comparison, overlap_list in exp.overlaps.items():
        comp_dir = make_folder(f'{out_dir}{comparison}_Overlap/')
        peakset = 'overlap_peak' if 'none' in [exp.sample_files[condition]['idr_optimal_peak'] for condition in overlap_list] else 'idr_optimal_peak'

        if (peakset == 'overlap_peak') and ('none' in [exp.sample_files[condition]['overlap_peak'] for condition in overlap_list]):
            output(f'ENCODE processing did not finish for at least one of the samples in the comparison {comparison}.  Skipping overlap...', log_file=exp.log_file, run_main=exp.run_main)
            with open(f'{comp_dir}SKIPPING_OVERLAP.txt', 'w') as file:
                file.write('Cannot find the peaks for at least one sample.  Skipping overlap...')
            continue

        for condition in overlap_list:
            exp.sample_files[condition]['peaktype'] = peakset

        bed_dict = {condition: load_bedtool(exp.sample_files[condition][peakset]) for condition in overlap_list}

        genome_list = exp.IPs.loc[exp.IPs['Condition'].isin(overlap_list), 'Genome'].unique().tolist()
        if len(genome_list) > 1:
            output(f'Cannot overlap peaks from different genomes for {condition}.', log_file=exp.log_file, run_main=exp.run_main)
            with open(f'{comp_dir}SKIPPING_OVERLAP.txt', 'w') as file:
                file.write('Cannot overlap peaks from different genomes for this condition.')
            continue
        else:
            genome = genome_list[0]

        if len(overlap_list) == 2:
            exp.overlap_results[comparison] = overlap_two(bed_dict, comparison, comp_dir, exp.log_file, genome=genome, run_main=exp.run_main)
        elif len(overlap_list) == 3:
            exp.overlap_results[comparison] = overlap_three(bed_dict, comparison, comp_dir, exp.log_file, genome=genome, run_main=exp.run_main)
        else:
            output(f'Cannot overlap more than three samples for {condition}.', log_file=exp.log_file, run_main=exp.run_main)
            with open(f'{comp_dir}SKIPPING_OVERLAP.txt', 'w') as file:
                file.write('Cannot overlap more than three samples.')
            continue

    output(f'Overlap analysis complete: {datetime.now():%Y-%m-%d %H:%M:%S}\n', log_file=exp.log_file, run_main=exp.run_main)

    return exp


def annotation(exp):

    from requests.exceptions import RetryError
    from time import sleep

    out_dir = make_folder(f'{exp.scratch}/Annotated/')

    condition_list = exp.IPs['Condition'].unique().tolist()
    for condition in condition_list:
        if 'peaktype' not in list(exp.sample_files[condition].keys()):
            peakset = 'overlap_peak' if exp.sample_files[condition]['idr_optimal_peak'] == 'none' else 'idr_optimal_peak'
            exp.sample_files[condition]['peaktype'] = peakset

    peakfiles = {condition: read_pd(exp.sample_files[condition][exp.sample_files[condition]['peaktype']])
                 for condition in exp.IPs['Condition'].unique().tolist()
                 if exp.sample_files[condition][exp.sample_files[condition]['peaktype']] != 'none'
                 }

    for condition, file in peakfiles.items():
        genome = exp.IPs.loc[exp.IPs.Condition == condition, 'Genome'].unique().tolist()[0]

        cond_dir = make_folder(f'{out_dir}{condition}/')
        anno_results = annotate_peaks({condition: file}, cond_dir, genome, db='UCSC', check=False, log_file=exp.log_file, run_main=exp.run_main)[f'{condition}_annotated']
        anno_list = anno_results.SYMBOL.unique().tolist()

        try:
            sleep(1)
            enrichr(anno_list, f'enrichr_{condition}', cond_dir, scan=None, max_terms=10, figsize=(12, 6), run_main=exp.run_main)
        except RetryError:
            output(f'No stable enrichr connection.  Skipping enrichr for {condition}.', log_file=exp.log_file, run_main=exp.run_main)

        exp.anno_results = {**exp.anno_results, **anno_results}

    return exp


def annotate_peaks(dict_of_dfs, folder, genome, log_file, db='UCSC', check=False, run_main=False):
    '''
    Annotate a dictionary of dataframes from bed files to the genome using ChIPseeker and Ensembl annotations.

    Inputs
    ------
    dict_of_beds: dictionary of bed files
    folder: output folder
    genome: hg38, hg19, mm10
    db: default UCSC, but can also accept Ensembl
    check: bool. checks whether annotation file already exists

    Returns
    -------
    dictionary of annotated bed files as dataframe

    '''
    pandas2ri.activate()

    ri.set_writeconsole_regular(rout_write)
    ri.set_writeconsole_warnerror(rout_write)

    folder = make_folder(folder)

    chipseeker = importr('ChIPseeker')
    genomicFeatures = importr('GenomicFeatures')
    makeGR = ro.r("makeGRangesFromDataFrame")

    check_df = {key: os.path.isfile(f'{folder}{key.replace(" ", "_")}_annotated.xlsx') for key in dict_of_dfs.keys()}
    return_bool = False not in set(check_df.values())
    if return_bool & check:
        return {f'{key}_annotated': pd.from_excel(f'{folder}{key.replace(" ", "_")}_annotated.xlsx') for key in dict_of_dfs.keys()}

    if db.lower() == 'ucsc':
        species = ('Mmusculus' if genome.lower() == 'mm10' else 'Hsapiens')
        TxDb = importr(f'TxDb.{species}.UCSC.{genome.lower()}.knownGene')
        txdb = ro.r(f'txdb <- TxDb.{species}.UCSC.{genome.lower()}.knownGene')
    elif db.lower() == 'ensembl':
        pwd = 'todo'
        loadDb = ro.r('loadDb')
        txdb = loadDb(pwd.format(genome.lower()))
    else:
        raise ValueError('UCSC or Ensembl only.')

    if genome.lower() == 'mm10':
        annoDb = importr('org.Mm.eg.db')
        anno = 'org.Mm.eg.db'
    elif genome.lower() == 'hg38' or genome.lower() == 'hg19':
        annoDb = importr('org.Hs.eg.db')
        anno = 'org.Hs.eg.db'

    return_dict = {}

    output('Annotating Peaks...', log_file=log_file, run_main=run_main)
    for key, df in dict_of_dfs.items():
        if check & check_df[key]:
            return_dict[f'{key}_annotated'] = pd.from_excel(f'{folder}{key.replace(" ", "_")}_annotated.xlsx')
        else:
            col_len = len(df.columns)
            df.columns = ["chr", "start", "end"] + list(range(col_len - 3))
            GR = makeGR(df)
            GR_anno = chipseeker.annotatePeak(GR, overlap='all', TxDb=txdb, annoDb=anno)
            return_dict[f'{key}_annotated'] = ro.pandas2ri.ri2py(chipseeker.as_data_frame_csAnno(GR_anno))
            return_dict[f'{key}_annotated'].to_excel(f'{folder}{key.replace(" ", "_")}_annotated.xlsx', index=None)

    return return_dict


def overlap_two(bed_dict, overlap_name, out_folder, log_file, genome=None, run_main=False):
    '''
    Takes a dictionary of two bed-like format files.
    Merges all overlapping peaks for each bed into a master file.
    Intersects beds to merged master file.
    Performs annotations with ChIPseeker if genome is specified.
    Plots venn diagrams of peak overlaps
    If genome is specified, also plots venn diagrams of annotated gene sets.

    Inputs
    ------
    bed_dict:  dictionary of BedTool files
    genome: 'hg38','hg19','mm10'

    Returns
    -------
    Returns a dictionary of dataframes from unique and overlap peaks.
    If genome is specified, includes a dictionary of annotated peaks.
    '''

    names = list(bed_dict.keys())

    out_folder = make_folder(out_folder)

    output(f'Output files for {overlap_name} are found in {out_folder}', log_file=log_file, run_main=run_main)

    masterfile = bed_dict[names[0]].cat(bed_dict[names[1]]).sort().merge()
    sorted_dict = {key: bed.sort().merge() for key, bed in bed_dict.items()}
    overlap_dict = {'overlap': masterfile.intersect(sorted_dict[names[0]]).intersect(sorted_dict[names[1]])}
    for key, bed in sorted_dict.items():
        other = {other_key: other_bed for other_key, other_bed in sorted_dict.items() if other_key != key}
        overlap_dict[f'{key}_unique_peak'] = masterfile.intersect(sorted_dict[key]).intersect(list(other.values())[0], v=True)

    for key, bed in overlap_dict.items():
        if len(bed) == 0:
            open(f'{out_folder}{key.replace(" ", "_")}-unique-peaks-from-mergedPeaks.bed', 'w').close()  # Can't convert empty bed file to dataframe
        else:
            bed2df(bed).to_csv(f'{out_folder}{key.replace(" ", "_")}-unique-peaks-from-mergedPeaks.bed', header=None, index=None, sep="\t")

    overlap_numbers = pd.Series({names[0]: len(overlap_dict[f'{names[0]}_unique_peak']),
                                 names[1]: len(overlap_dict[f'{names[1]}_unique_peak']),
                                 'overlap': len(overlap_dict['overlap'])
                                 },
                                index=[names[0], names[1], 'overlap']
                                )

    # Venn
    plot_venn2(overlap_numbers, overlap_name.replace('_', ' '), out_folder)
    out_result(f'{out_folder}venn_plot/{overlap_name.replace(" ","_")}-overlap.png', f"{overlap_name.replace('_',' ')} Peak Venn Overlap", run_main=run_main)

    if bool(genome):
        # output(f'Annotating overlaping peaks for {overlap_name.replace("_"," ")}...', log_file)
        # Annotate with ChIPseeker
        unikey = '{}_unique'
        unianno = '{}_unique_annotated'
        return_dict = annotate_peaks({unikey.format(key): bed2df(bed) for key, bed in overlap_dict.items() if len(bed) > 0}, out_folder, genome=genome, log_file=log_file, run_main=run_main)
        for key, bed in overlap_dict.items():
            if len(bed) == 0:
                return_dict[unianno.format(key)] = None

        Set1_unique = set() if return_dict[unianno.format(f'{names[0]}_unique_peak')] is None else set(return_dict[unianno.format(f'{names[0]}_unique_peak')].SYMBOL.unique().tolist())
        Set2_unique = set() if return_dict[unianno.format(f'{names[1]}_unique_peak')] is None else set(return_dict[unianno.format(f'{names[1]}_unique_peak')].SYMBOL.unique().tolist())
        Overlap_Set = set() if return_dict[unianno.format('overlap')] is None else set(return_dict[unianno.format('overlap')].SYMBOL.unique().tolist())

        venn2_dict = {names[0]: (Set1_unique | Overlap_Set),
                      names[1]: (Set2_unique | Overlap_Set)
                      }

        plot_name = f'{overlap_name.replace("_"," ")} Annotated Gene'
        plot_venn2_set(venn2_dict, plot_name, out_folder)
        out_result(f'{out_folder}venn_plot/{plot_name.replace(" ","_")}-overlap.png', f"{overlap_name.replace('_',' ')} Venn Annotated Gene Overlap", run_main=run_main)

        gene_overlaps = {}
        gene_overlaps[f'{names[0]}_unique_genes'] = Set1_unique - (Set2_unique | Overlap_Set)
        gene_overlaps[f'{names[1]}_unique_genes'] = Set2_unique - (Set1_unique | Overlap_Set)
        gene_overlaps['Overlap_Gene_Set'] = (Set1_unique & Set2_unique) | Overlap_Set

        return_dict = {key: bed2df(bed) for key, bed in overlap_dict.items()}

        for key, item in gene_overlaps.items():
            return_dict[key] = item

    else:
        return_dict = {key: bed2df(bed) for key, bed in overlap_dict.items()}

    return return_dict


def overlap_three(bed_dict, overlap_name, out_folder, log_file, genome=None, run_main=False):
    '''
    Takes a dictionary of three bed-like format files.
    Merges all overlapping peaks for each bed into a master file.
    Intersects beds to merged master file.
    Performs annotations with ChIPseeker if genome is specified.
    Plots venn diagrams of peak overlaps
    If genome is specified, also plots venn diagrams of annotated gene sets.

    Inputs
    ------
    bed_dict:  dictionary of BedTool files
    genome: 'hg38','hg19','mm10'

    Returns
    -------
    Returns a dictionary of dataframes from unique and overlap peaks.
    If genome is specified, includes a dictionary of annotated peaks.
    '''
    from collections import OrderedDict

    names = list(bed_dict.keys())

    out = make_folder(out_folder)

    output(f'Output files are found in {out}', log_file=log_file, run_main=run_main)
    output(f'A: {names[0]}, B: {names[1]}, C: {names[2]}', log_file=log_file, run_main=run_main)
    with open(f'{out}README.txt', 'w') as file:
        file.write('All peaks are unique, meaning that each peak is in only one group.\n')
        file.write('Capital letter means this sample peak is included in the overlap.\n')
        file.write('Lowercase letter means the sample is excluded in the overlap.\n\n')
        file.write(f'A: {names[0]}\nB: {names[1]}\nC: {names[2]}')

    master = bed_dict[names[0]].cat(bed_dict[names[1]]).cat(bed_dict[names[2]]).sort().merge()

    A = bed_dict[names[0]].sort().merge()
    B = bed_dict[names[1]].sort().merge()
    C = bed_dict[names[2]].sort().merge()

    sorted_dict = OrderedDict({'master': master, 'A': A, 'B': B, 'C': C})
    sorted_dict['Abc'] = master.intersect(A).intersect(B, v=True).intersect(C, v=True)
    sorted_dict['aBc'] = master.intersect(B).intersect(A, v=True).intersect(C, v=True)
    sorted_dict['ABc'] = master.intersect(A).intersect(B).intersect(C, v=True)
    sorted_dict['abC'] = master.intersect(C).intersect(A, v=True).intersect(B, v=True)
    sorted_dict['AbC'] = master.intersect(A).intersect(C).intersect(B, v=True)
    sorted_dict['aBC'] = master.intersect(B).intersect(C).intersect(A, v=True)
    sorted_dict['ABC'] = master.intersect(A).intersect(B).intersect(C)

    labTup = tuple(key for key in sorted_dict.keys())
    lenTup = tuple(len(bed) for bed in sorted_dict.values())

    output(f'{labTup}\n{lenTup}', log_file=log_file, run_main=run_main)

    plot_venn3_counts(lenTup[4:], names, f'{overlap_name} Peak', out)
    out_result(f'{out}venn_plot/{overlap_name}_Peak-overlap.png', f"{overlap_name} Peak Venn Overlap", run_main=run_main)

    for key, bed in sorted_dict.items():
        if len(bed) == 0:
            open(f'{out}{key.replace(" ", "_")}-peaks-from-mergedPeaks.bed', 'w').close()  # Can't convert empty bed file to dataframe
        else:
            bed2df(bed).to_csv(f"{out}{key.replace(' ', '_')}-peaks-from-mergedPeaks.bed", header=None, index=None, sep="\t")

    if bool(genome):
        output('Annotating ovelapped peaks...', log_file=log_file)
        unikey = '{}_unique'
        unianno = '{}_unique_annotated'
        return_dict = annotate_peaks({unikey.format(key): bed2df(bed) for key, bed in sorted_dict.items() if len(bed) > 0}, out, genome=genome, log_file=log_file, run_main=run_main)
        for key, bed in sorted_dict.items():
            if len(bed) == 0:
                return_dict[unianno.format(key)] = None

        Set1 = set() if return_dict[unianno.format('A')] is None else set(return_dict[unianno.format('A')].SYMBOL.unique().tolist())
        Set2 = set() if return_dict[unianno.format('B')] is None else set(return_dict[unianno.format('B')].SYMBOL.unique().tolist())
        Set3 = set() if return_dict[unianno.format('C')] is None else set(return_dict[unianno.format('C')].SYMBOL.unique().tolist())

        plot_venn3_set({names[0]: Set1, names[1]: Set2, names[2]: Set3}, f'{overlap_name}_annotated_genes', out)
        out_result(f'{out}venn_plot/{overlap_name}_annotated_genes-overlap.png', f"{overlap_name.replace('_',' ')} Gene Venn Overlap", run_main=run_main)

    return_sorted_dict = {key: bed2df(bed) for key, bed in sorted_dict.items()}

    return return_sorted_dict if genome is None else {**return_sorted_dict, **return_dict}


def enrichr(gene_list, description, out_dir, scan=None, max_terms=10, figsize=(12, 6), run_main=False):
    '''
    Performs GO Molecular Function, GO Biological Process and KEGG enrichment on a gene list.
    Uses enrichr.

    Inputs
    ------
    gene_list: list of genes to perform enrichment on
    description: string description for title
    out_dir: output director
    scan: dictionary with additional enrichr dbs to scan (http://amp.pharm.mssm.edu/Enrichr/#stats)
    max_terms: limit return plot to this max
    load: load results
    figsize: change fig size

    Returns
    -------

    None

    '''

    out_dir = make_folder(out_dir)

    testscan = {'KEGG': 'KEGG_2016',
                'GO_biological_process': 'GO_Biological_Process_2017b',
                'ChIP-X_Consensus_TFs': 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
                'ChEA': 'ChEA_2016',
                'OMIM_Disease': 'OMIM_Disease'
                }

    if isinstance(scan, dict):
        testscan = {**testscan, **scan}

    for nick, name in testscan.items():
        try:
            gseapy.enrichr(gene_list=gene_list,
                           figsize=figsize,
                           top_term=max_terms,
                           description=f'{description}_{nick}',
                           gene_sets=name,
                           outdir=out_dir,
                           format='png'
                           )

            out_result(f'{out_dir}{nick}.{name}.enrichr.reports.png', f'Enrichr: {nick} for {description}', run_main=run_main)

        except:
            output(f'Error in enrichr submission for {description} {nick}. Gene list is {len(gene_list)}.  Enrichr has problems > 6000 genes.', run_main=False)

    out_list = pd.DataFrame({'Gene Name': gene_list}, index=range(len(gene_list)))
    out_list.to_excel(f'{out_dir}{description}_genes.xlsx', index=None)


def heatmaps(exp):
    return exp
