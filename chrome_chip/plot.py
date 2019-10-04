#!/usr/bin/env python3

import glob
import re

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
import seaborn as sns
from scipy import stats

from chrome_chip.common import val_folder, out_result, output, make_folder, submission_prepend


def plot_col(df, title, ylabel, out='', xy=(None, None), xticks=[''], plot_type=['violin', 'swarm'], pvalue=False, compare_tags=None, log_file=None, run_main=False):
    '''
    One or two column boxplot from dataframe.  Titles x axis based on column names.

    Inputs
    ------
    df: dataframe (uses first two columns)
    title: string of title
    ylabel: string of y label
    xy: If specified, will x is the label column and y is the data column. (default: (None,None): Data separated into two columns).
    xticks: list of xtick names (default is none)
    pvalue: bool to perform ttest (default is False).  Will only work if xy=(None,None) or ther are only two labels in x.
    plot_type: list of one or more: violin, box, swarm (default=violin)
    compare_tags:  if xy and pvalue is specified and there are more than two tags in x, specify the tags to compare. eg. ['a','b']
    out: out parent directory.  if none returns into colplot/
    log_file: log_file

    Returns
    ------
    None
    '''

    out = make_folder(f'{val_folder(out)}plots/')

    plt.clf()
    sns.set(context='paper', font='Arial', font_scale=2, style='white', rc={'figure.dpi': 300, 'figure.figsize': (5, 6)})

    if type(plot_type) != list:
        plot_type = plot_type.split()
    lower_plot_type = [x.lower() for x in plot_type]

    if len(lower_plot_type) == 0:
        raise IOError('Input a plot type.')
    elif True not in {x in lower_plot_type for x in ['violin', 'box', 'swarm']}:
        raise IOError('Did not recognize plot type.')

    if 'swarm' in lower_plot_type:
        if xy == (None, None):
            fig = sns.swarmplot(data=df, color='black', s=4)
        else:
            fig = sns.swarmplot(data=df, x=xy[0], y=xy[1], color='black', s=4)
    if 'violin' in lower_plot_type:
        if xy == (None, None):
            fig = sns.violinplot(data=df)
        else:
            fig = sns.violinplot(data=df, x=xy[0], y=xy[1])
    if 'box' in lower_plot_type:
        if xy == (None, None):
            fig = sns.boxplot(data=df)
        else:
            fig = sns.boxplot(data=df, x=xy[0], y=xy[1])

    fig.yaxis.set_label_text(ylabel)
    fig.set_title(title)
    if xticks:
        fig.xaxis.set_ticklabels(xticks)
        fig.xaxis.set_label_text('')
        for tick in fig.xaxis.get_ticklabels():
            tick.set_fontsize(12)

    if pvalue:
        if xy == (None, None):
            _, pvalue = stats.ttest_ind(a=df.iloc[:, 0], b=df.iloc[:, 1])
            compare_tags = df.columns
        else:
            _, pvalue = stats.ttest_ind(a=df[df[xy[0]] == compare_tags[0]][xy[1]], b=df[df[xy[0]] == compare_tags[1]][xy[1]])
        fig.text(s=f'p-value = {pvalue:.03g}, {compare_tags[0]} v {compare_tags[1]}', x=0, y=-.12, transform=fig.axes.transAxes, fontsize=12)

    sns.despine()
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.17, top=0.9)
    plt.savefig(f"{out}{title.replace(' ', '_')}.png", dpi=300)
    if run_main:
        plt.close()

    out_result(f"{out}{title.replace(' ', '_')}.png", f'{title} Plot', run_main=run_main)
    output(f"{title.replace(' ', '_')}.png found in {out}", log_file=log_file, run_main=run_main)


def deeptools(regions, signals, matrix_name, out_name, pegasus_folder, title='', bps=(1500, 1500, 4000), type='center', scaled_names=('TSS', 'TES'), make=('matrix', 'heatmap', 'heatmap_group', 'profile', 'profile_group')):
    '''
    Inputs
    ------
    regions: dictionary {'region_name':'/path/to/ssh/bedfile'}
    signals: dictionary {'signal_name':'/path/to/ssh/bigwigfile'}
    matrix_name: string of matrix name or matrix to be named (before .matrix.gz)
    out_name: name for output file
    tite: plot title (optional)
    bps: tuple of region width on either side of center or scaled.  center ignores last number.  default is (1500,1500,4000)
    type: 'center' or 'scaled'
    scaled_names: optional names for scaled start and end (default ('TSS','TES'))
    make: tuple of deeptool commands.  options: matrix, heatmap, heatmap_group, profile, profile_group
    copy: bool.  Copy region and signal files to peagasus
    copy_folder: folder to copy into

    Returns
    -------
    string of commands for ssh_job

    '''
    pegasus_folder = make_folder(pegasus_folder)

    make_lower = [x.lower() for x in make]

    if type.lower() == 'center':
        deepMat = 'reference-point --referencePoint center'
        deepHeat = "--refPointLabel 'Peak Center'"
        deepProf = "--refPointLabel 'Peak Center'"
    else:
        deepMat = f'scale-regions --regionBodyLength {str(bps[2])}'
        deepHeat = f'--startLabel {scaled_names[0]} --endLabel {scaled_names[1]}'
        deepProf = f'--startLabel {scaled_names[0]} --endLabel {scaled_names[1]}'

    cmd_list = [submission_prepend()]

    pegasus_region_path = ' '.join([f"{pegasus_folder}{region_path.split('/')[-1]}" for region_path in regions.values()])
    pegasus_signal_path = ' '.join([f"{pegasus_folder}{signal_path.split('/')[-1]}" for signal_path in signals.values()])

    if 'matrix' in make_lower:
        signal_name = ' '.join([signal_name for signal_name in signals.keys()])
        computeMatrix = f"computeMatrix {deepMat} -a {str(bps[0])} -b {str(bps[1])} -p 4 -R {pegasus_region_path} -S {pegasus_signal_path} --samplesLabel {signal_name} -o {matrix_name}.matrix.gz"
        cmd_list.append(computeMatrix)

    if 'heatmap' in make_lower or 'heatmap_group' in make_lower:
        region_name = ' '.join([region_name for region_name in regions.keys()])
        plotHeatmap_base = f"plotHeatmap -m {matrix_name}.matrix.gz --dpi 300 {deepHeat} --regionsLabel {region_name} --plotTitle '{title.replace('_', ' ')}' --whatToShow 'heatmap and colorbar' --colorMap Reds -out {out_name}_heatmap"
        if 'heatmap' in make_lower:
            cmd_list.append(f"{plotHeatmap_base}.png")
        if 'heatmap_group' in make_lower:
            cmd_list.append(f"{plotHeatmap_base}_perGroup.png --perGroup")

    if 'profile' in make_lower or 'profile_group' in make_lower:
        region_name = ' '.join([region_name for region_name in regions.keys()])
        plotProfile_base = f"plotProfile -m {matrix_name}.matrix.gz --dpi 300 {deepProf} --plotTitle '{title.replace('_', ' ')}' --regionsLabel {region_name} -out {out_name}_profile"
        if 'profile' in make_lower:
            cmd_list.append(f"{plotProfile_base}.png")
        if 'profile_group' in make_lower:
            cmd_list.append(f"{plotProfile_base}_perGroup.png --perGroup")

    return cmd_list


def plot_venn2(Series, overlap_name, folder):
    '''
    Series with with overlaps 10,01,11
    Plots a 2 way venn.
    Saves to file.
    '''

    folder = make_folder(f"{val_folder(folder)}venn_plot")

    plt.clf()
    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    # make venn
    venn_plot = venn2(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]), set_labels=[name.replace('_', ' ') for name in Series.index.tolist()])
    patch = ['10', '01', '11']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')

    c = venn2_circles(subsets=(Series.iloc[0], Series.iloc[1], Series.iloc[2]))
    colors = ['green', 'blue']
    for circle, color in zip(c, colors):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(2)

    plt.title(overlap_name.replace('_', ' ') + " overlaps")
    plt.tight_layout()
    plt.savefig(f'{folder}{overlap_name.replace(" ", "_")}-overlap.svg')
    plt.savefig(f'{folder}{overlap_name.replace(" ", "_")}-overlap.png', dpi=300)


def plot_venn2_set(dict_of_sets, overlap_name, folder):
    '''
    Plots a 2 way venn from a dictionary of sets
    Saves to file.

    Inputs
    ------
    dict_of_sets: dictionary of sets to overlap
    overlap_name: string with name of overlap
    folder: output folder

    Returns
    -------
    None

    '''
    folder = make_folder(f"{val_folder(folder)}venn_plot")

    plt.clf()
    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    set_list = []
    set_names = []
    for name, setlist in dict_of_sets.items():
        set_list.append(setlist)
        set_names.append(name.replace('_', ' '))

    # make venn
    venn_plot = venn2(subsets=set_list, set_labels=set_names)
    patch = ['10', '01', '11']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')

    c = venn2_circles(subsets=set_list)
    colors = ['green', 'blue']
    for circle, color in zip(c, colors):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(2)

    plt.title(f"{overlap_name.replace('_', ' ')} overlaps")
    plt.tight_layout()
    plt.savefig(f'{folder}{overlap_name.replace(" ", "_")}-overlap.svg')
    plt.savefig(f'{folder}{overlap_name.replace(" ", "_")}-overlap.png', dpi=300)


def plot_venn3_set(dict_of_sets, overlap_name, folder):
    '''
    Makes 3 way venn from 3 sets.
    Saves to file.

    Inputs
    ------
    dict_of_sets: dictionary of sets to overlap
    overlap_name: string with name of overlap
    folder: output folder

    Returns
    -------
    None

    '''
    folder = make_folder(f"{val_folder(folder)}venn_plot")

    plt.clf()
    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    set_list = []
    set_names = []
    for name, setlist in dict_of_sets.items():
        set_list.append(setlist)
        set_names.append(name.replace('_', ' '))

    # make venn
    venn_plot = venn3(subsets=set_list, set_labels=set_names)
    patch = ['100', '110', '101', '010', '011', '001', '111']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')

    # make
    c = venn3_circles(subsets=set_list)
    colors_list = ['green', 'blue', 'grey']
    for circle, color in zip(c, colors_list):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(4)

    plt.title(f"{overlap_name.replace('_', ' ')} Overlaps")
    plt.tight_layout()
    plt.savefig(f"{folder}{overlap_name.replace(' ', '_')}-overlap.svg")
    plt.savefig(f"{folder}{overlap_name.replace(' ', '_')}-overlap.png", dpi=300)
    plt.close()


def plot_venn3_counts(element_list, set_labels, overlap_name, folder):
    '''
    Plot three way venn based on counts of specific overlaping numbers.
    Saves to file.

    Inputs
    ------
    element_list: tuple with counts of the the overlaps from (Abc,aBc,ABc,abC,AbC,ABC)
    set_labels: list or tuple with names of the overlaps ('A','B','C')
    overlap_name: string with name of overlap
    folder: output folder

    Returns
    -------
    None

    '''
    folder = make_folder(f"{val_folder(folder)}venn_plot")

    plt.clf()
    plt.figure(figsize=(7, 7))

    font = {'family': 'sans-serif',
            'weight': 'normal',
            'size': 16,
            }

    plt.rc('font', **font)

    # make venn
    venn_plot = venn3(subsets=element_list, set_labels=[name.replace('_', ' ') for name in set_labels])
    patch = ['100', '110', '101', '010', '011', '001', '111']
    for p in patch:
        if venn_plot.get_patch_by_id(p):
            venn_plot.get_patch_by_id(p).set_color('none')
            venn_plot.get_patch_by_id(p).set_alpha(.4)
            venn_plot.get_patch_by_id(p).set_edgecolor('none')

    # make
    c = venn3_circles(subsets=element_list)
    colors_list = ['green', 'blue', 'grey']
    for circle, color in zip(c, colors_list):
        circle.set_edgecolor(color)
        circle.set_alpha(0.8)
        circle.set_linewidth(4)

    plt.title(f"{overlap_name.replace('_', ' ')} Overlaps")
    plt.tight_layout()
    plt.savefig(f"{folder}{overlap_name.replace(' ', '_')}-overlap.svg")
    plt.savefig(f"{folder}{overlap_name.replace(' ', '_')}-overlap.png", dpi=300)


def extract_AQUAS_report_data(base_folder, out_folder='', histone=False, replicate=False):
    '''
    Inputs
    -----
    base_folder:  AQUAS results folder.  Will use subfolders for sample name and look for report in those subfolders.
    replicate: Whether the ChIPseq was performed as a repliate or not.

    Returns
    -----
    DataFrame of results
    '''

    reports = glob.glob(f'{base_folder}/*/*report.html')
    out_folder = val_folder(out_folder)
    base_folder = val_folder(base_folder)

    if replicate is True:
        raise AssertionError('Not set up for replicates yet.')

    results_df = pd.DataFrame(index=['Percent_mapped', 'Mapped_Reads', 'Fraction_Duplicated', 'S_JS_Distance', 'PBC1', 'RSC', 'Raw_Peak_Number', 'N_optimal_overlap_peaks', 'FrIP_IDR', 'N_IDR_peaks'])
    for file in reports:
        name = re.findall(r'.*/(.*)_report.html', file)[0]
        report = pd.read_html(file)
        series = pd.Series()
        series['Percent_mapped'] = report[1].iloc[7, 1]
        series['Mapped_Reads'] = report[2].iloc[5, 1]
        series['Fraction_Duplicated'] = report[3].iloc[7, 1]
        series['S_JS_Distance'] = report[4].iloc[7, 1]
        series['PBC1'] = report[5].iloc[6, 1]
        series['RSC'] = report[6].iloc[8, 1]
        series['Raw_Peak_Number'] = report[7].iloc[0, 1]
        series['N_optimal_overlap_peaks'] = report[10].iloc[4, 1]
        if histone is False:
            series['FrIP_IDR'] = report[11].iloc[0, 1]
            series['N_IDR_peaks'] = report[12].iloc[4, 1]
        results_df[name] = series

    for index in results_df.index.tolist():
        plot_col(results_df.loc[index], out=out_folder, title=f'{index}', ylabel=index.replace('_', ' '), plot_type=['violin', 'swarm'])

    return results_df


def scale_factor(spike1_tags, spike2_tags, IP1_tags, IP2_tags):
    '''
    rep1 should be the lowest of number of tags in all samples in the comparison.
    returns scaled ratio of IP2.
    '''
    correction_factor = int(spike2_tags) / int(spike1_tags)
    scaled_IP2 = int(IP2_tags) / correction_factor
    IP2_ratio = scaled_IP2 / int(IP1_tags)

    return IP2_ratio


def spike_in_plot(spike_df, sample_list, description, out_dir):
    '''
    df: spike-in df with spike-in reads counts in column 'spike_reads' and genome counts in 'genome_reads'.
    sample_list:  list of samples to compare

    '''
    import seaborn as sns

    df = spike_df.loc[sample_list].copy()
    replicates = True if len(df.Replicate.unique().tolist()) > 1 else False

    base_idx = df.genome_reads.astype(int).idxmin()
    base_spike = df.loc[base_idx, 'spike_reads']
    base_genome = df.loc[base_idx, 'genome_reads']

    df['Tag Ratio'] = df[['spike_reads', 'genome_reads']].apply(lambda x: scale_factor(base_spike, x[0], base_genome, x[1]), axis=1)

    sns.set(context='notebook', style='white', font='Arial')
    if replicates:
        sns.boxplot(x='Condition', y='Tag Ratio', data=df, palette='pastel')
        sns.swarmplot(x='Condition', y='Tag Ratio', hue='Replicate', data=df, size=8)
    else:
        sns.barplot(x='Condition', y='Tag Ratio', data=df)

    sns.mpl.pyplot.ylabel('Normalized Genome Read Ratio')
    sns.mpl.pyplot.title(description.replace('_', ' '))
    sns.despine()
    sns.mpl.pyplot.savefig(f'{out_dir}{description.replace(" ", "_")}.png')
    sns.mpl.pyplot.savefig(f'{out_dir}{description.replace(" ", "_")}.svg')

    return f'{out_dir}{description.replace(" ", "_")}.png'
