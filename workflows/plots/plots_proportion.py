import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import seaborn as sns
sns.set_style("darkgrid")

from sklearn import metrics

# calculate mean accuracies from accuracy tables and save as new files
sstar_1src_accuracy_proportion = pd.read_csv(snakemake.input.sstar_1src_accuracy_proportion, sep="\t").dropna()
sprime_1src_accuracy_proportion = pd.read_csv(snakemake.input.sprime_1src_accuracy_proportion, sep="\t").dropna()
skovhmm_1src_accuracy_proportion = pd.read_csv(snakemake.input.skovhmm_1src_accuracy_proportion, sep="\t").dropna()

sstar_1src_accuracy_proportion_grouped = sstar_1src_accuracy_proportion.groupby(['demography', 'scenario', 'sample', 'cutoff'], as_index=False)
sprime_1src_accuracy_proportion_grouped = sprime_1src_accuracy_proportion.groupby(['demography', 'sample', 'cutoff'], as_index=False)
skovhmm_1src_accuracy_proportion_grouped = skovhmm_1src_accuracy_proportion.groupby(['demography', 'sample', 'cutoff'], as_index=False)

sstar_1src_accuracy_proportion_mean = sstar_1src_accuracy_proportion_grouped.mean()
sstar_1src_accuracy_proportion_mean['precision_std'] = sstar_1src_accuracy_proportion_grouped.std()['precision']
sstar_1src_accuracy_proportion_mean['recall_std'] = sstar_1src_accuracy_proportion_grouped.std()['recall']

sprime_1src_accuracy_proportion_mean = sprime_1src_accuracy_proportion_grouped.mean()
sprime_1src_accuracy_proportion_mean['precision_std'] = sprime_1src_accuracy_proportion_grouped.std()['precision']
sprime_1src_accuracy_proportion_mean['recall_std'] = sprime_1src_accuracy_proportion_grouped.std()['recall']

skovhmm_1src_accuracy_proportion_mean = skovhmm_1src_accuracy_proportion_grouped.mean()
skovhmm_1src_accuracy_proportion_mean['precision_std'] = skovhmm_1src_accuracy_proportion_grouped.std()['precision']
skovhmm_1src_accuracy_proportion_mean['recall_std'] = skovhmm_1src_accuracy_proportion_grouped.std()['recall']

sstar_1src_accuracy_proportion_mean.to_csv(snakemake.output.sstar_1src_accuracy_proportion_mean, sep="\t", index=False)
sprime_1src_accuracy_proportion_mean.to_csv(snakemake.output.sprime_1src_accuracy_proportion_mean, sep="\t", index=False)
skovhmm_1src_accuracy_proportion_mean.to_csv(snakemake.output.skovhmm_1src_accuracy_proportion_mean, sep="\t", index=False)

sprime_1src_accuracy_proportion_mean['scenario'] = ['true'] * len(sprime_1src_accuracy_proportion_mean)
skovhmm_1src_accuracy_proportion_mean['scenario'] = ['true'] * len(skovhmm_1src_accuracy_proportion_mean)
sstar_1src_accuracy_proportion_mean['scenario'] = ['true'] * len(sstar_1src_accuracy_proportion_mean)

# define variables and layout of the plots
methods1 = ['sprime', 'skovhmm', 'sstar']
demography1 = ['HumanNeanderthal']
samples = ['nref_10_ntgt_1', 'nref_50_ntgt_1', 'nref_100_ntgt_1', 'nref_10_ntgt_10', 'nref_50_ntgt_10', 'nref_100_ntgt_10', 'nref_10_ntgt_50', 'nref_50_ntgt_50', 'nref_100_ntgt_50', 'nref_10_ntgt_100', 'nref_50_ntgt_100', 'nref_100_ntgt_100']
scenarios = ['true']
accuracy1 = {
    'sprime': sprime_1src_accuracy_proportion_mean,
    'skovhmm': skovhmm_1src_accuracy_proportion_mean,    
    'sstar': sstar_1src_accuracy_proportion_mean,
}


fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(8,6.5), dpi=300)
gridspec = axs[0,0].get_subplotspec().get_gridspec()
for a in axs[:, 1]:
    a.remove()
for c in axs[1, 0:1]:
    c.remove()


markers = {
    'nref_10_ntgt_1': {'symbol':'.', 'size': 6},
    'nref_50_ntgt_1': {'symbol':'*', 'size': 6},
    'nref_100_ntgt_1': {'symbol':'P', 'size': 6},
    'nref_10_ntgt_10': {'symbol':'p', 'size': 4},
    'nref_50_ntgt_10': {'symbol':'d', 'size': 4},
    'nref_100_ntgt_10': {'symbol':'D', 'size': 4},
    'nref_10_ntgt_50': {'symbol':'H', 'size': 4},
    'nref_50_ntgt_50': {'symbol':'+', 'size': 4},
    'nref_100_ntgt_50': {'symbol':'s', 'size': 4},
    'nref_10_ntgt_100': {'symbol':'X', 'size': 4},
    'nref_50_ntgt_100': {'symbol':'x', 'size': 4},
    'nref_100_ntgt_100': {'symbol':'o', 'size': 4},
}

colors = {
    'sprime': 'blue',
    'skovhmm': 'green',
    'sstar': {'true': 'orange'},
}


linestyles = {
    'true': 'solid',
}


titles = {
    'HumanNeanderthal': 'Human-Neanderthal model\n(Gower et al. 2021)',
}


# create file with calculated AUC (area under the precision-recall curve) score values
auc = open(snakemake.output.auc1, 'w')
auc.write('method\tdemography\tscenario\tsample\tAUC\n')

# create precision-recall plot
j = 0
for d in demography1:
    for s in samples:
        for sc in scenarios:
            for m in methods1:
                if (m == 'sstar') or (m == 'skovhmm'):
                    if (s == 'nref_10_ntgt_10') or (s == 'nref_50_ntgt_10') or (s == 'nref_100_ntgt_10') or (s == 'nref_10_ntgt_50') or (s == 'nref_50_ntgt_50') or (s == 'nref_100_ntgt_50') or (s == 'nref_10_ntgt_100') or (s == 'nref_50_ntgt_100') or (s == 'nref_100_ntgt_100'): continue
                if m == 'sstar': color = colors[m][sc]
                else: color = colors[m]
                df = accuracy1[m][
                        (accuracy1[m]['demography'] == d) &
                        (accuracy1[m]['sample'] == s) &
                        (accuracy1[m]['scenario'] == sc)
                    ].sort_values(by='recall', ascending=False)
                recall = df['recall']
                precision = df['precision']
                if (m == 'sprime') or (m == 'skovhmm'):
                    if sc != 'true': continue
                # calculate AUC scores
                auc_score = metrics.auc(recall/100, precision/100)
                auc.write(f'{m}\t{d}\t{sc}\t{s}\t{auc_score}\n')
                axs[j,0].plot(recall, precision,
                    marker=markers[s]['symbol'], ms=markers[s]['size'],
                    c=color)

    axs[0,j].set_xlabel('Recall (%)', fontsize=8)
    axs[0,j].set_ylabel('Precision (%)', fontsize=8)
    axs[0,j].set_xlim([-5, 105])
    axs[0,j].set_ylim([-5, 105])
    axs[0,j].set_title(titles[d], fontsize=8, weight='bold')

auc.close()
    

# create figure legend
subfig = fig.add_subfigure(gridspec[:,1])
handles, labels = subfig.gca().get_legend_handles_labels()
sprime_line = plt.Line2D([0], [0], label='SPrime', color=colors['sprime'])
skovhmm_line = plt.Line2D([0], [0], label='SkovHMM', color=colors['skovhmm'])
sstar_line = plt.Line2D([0], [0], label='sstar (true model)', color=colors['sstar']['true'])
nref_10_ntgt_1 = plt.Line2D([0], [0], marker=markers['nref_10_ntgt_1']['symbol'],
                            ms=5, label='Nref=10, Ntgt=1', color='black', linewidth=0)
nref_50_ntgt_1 = plt.Line2D([0], [0], marker=markers['nref_50_ntgt_1']['symbol'],
                            ms=5, label='Nref=50, Ntgt=1', color='black', linewidth=0)
nref_100_ntgt_1 = plt.Line2D([0], [0], marker=markers['nref_100_ntgt_1']['symbol'],
                             ms=4, label='Nref=100, Ntgt=1', color='black', linewidth=0)
nref_10_ntgt_10 = plt.Line2D([0], [0], marker=markers['nref_10_ntgt_10']['symbol'],
                             ms=5, label='Nref=10, Ntgt=10', color='black', linewidth=0)
nref_50_ntgt_10 = plt.Line2D([0], [0], marker=markers['nref_50_ntgt_10']['symbol'],
                             ms=4, label='Nref=50, Ntgt=10', color='black', linewidth=0)
nref_100_ntgt_10 = plt.Line2D([0], [0], marker=markers['nref_100_ntgt_10']['symbol'],
                             ms=4, label='Nref=100, Ntgt=10', color='black', linewidth=0)
nref_10_ntgt_50 = plt.Line2D([0], [0], marker=markers['nref_10_ntgt_50']['symbol'],
                             ms=5, label='Nref=10, Ntgt=50', color='black', linewidth=0)
nref_50_ntgt_50 = plt.Line2D([0], [0], marker=markers['nref_50_ntgt_50']['symbol'],
                             ms=4, label='Nref=50, Ntgt=50', color='black', linewidth=0)
nref_100_ntgt_50 = plt.Line2D([0], [0], marker=markers['nref_100_ntgt_50']['symbol'],
                             ms=4, label='Nref=100, Ntgt=50', color='black', linewidth=0)
nref_10_ntgt_100 = plt.Line2D([0], [0], marker=markers['nref_10_ntgt_100']['symbol'],
                             ms=5, label='Nref=10, Ntgt=100', color='black', linewidth=0)
nref_50_ntgt_100 = plt.Line2D([0], [0], marker=markers['nref_50_ntgt_100']['symbol'],
                             ms=4, label='Nref=50, Ntgt=100', color='black', linewidth=0)
nref_100_ntgt_100 = plt.Line2D([0], [0], marker=markers['nref_100_ntgt_100']['symbol'],
                             ms=4, label='Nref=100, Ntgt=100', color='black', linewidth=0)
handles.extend([sprime_line, skovhmm_line, sstar_line,
                nref_10_ntgt_1, nref_50_ntgt_1, nref_100_ntgt_1, nref_10_ntgt_10, nref_50_ntgt_10, nref_100_ntgt_10, nref_10_ntgt_50, nref_50_ntgt_50, nref_100_ntgt_50, nref_10_ntgt_100, nref_50_ntgt_100, nref_100_ntgt_100,])
subfig.legend(handles=handles, fontsize=8, handlelength=2)

plt.savefig(snakemake.output.accuracy, bbox_inches='tight')
