import itertools
import json
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mtick
from cycler import cycler
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from collections import defaultdict
from pathlib import Path
from rich.status import Status
from sklearn import metrics as sklearn_metric

from .analysis_function import ClassificationResults, precision_recall_f1, precision_recall_f1_threshold, calculate_l1_distance

mpl.rcParams["mathtext.default"] = 'regular'
plt.style.use('default')
plt.rcParams['axes.axisbelow'] = True
plt.rcParams['axes.grid'] = True
custom_cycler = cycler(color=['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6'])
plt.rcParams['axes.prop_cycle'] = custom_cycler

list_uc = ['no_hit',
           'ambigious',
           'unclassified',
           'missing_genus'
           ]

logger = logging.getLogger('classification_agg')

update_classifier_names = {"bracken": "Bracken",
                           "ccmetagen": "CCMetagen",
                           "centrifuge": "Centrifuge",
                           "kaiju": "Kaiju",
                           "kma": "KMA",
                           "kraken2": "Kraken2",
                           'metaphlan': 'MetaPhlAn3',
                           "mmseqs2": "MMSeqs2",
                           'motus': 'mOTUs2'
                           }
list_abundance_classifiers = ['MetaPhlAn3', 'mOTUs2']

# Determines order of the classifiers in figures
classifier_order = ['Kraken2', 'Bracken', 'KMA', 'CCMetagen', 'Centrifuge', 'Kaiju', 'MMSeqs2', 'MetaPhlAn3', 'mOTUs2']

ALPH = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

DATASET_ID = {'bei_resources_HM-276D_com_B_even': 'BeiRes_276',
              'bei_resources_HM-277D_com_B_staggered': 'BeiRes_277D',
              'strain_madness_1': 'StrainMad_1',
              'strain_madness_2': 'StrainMad_2',
              'strain_madness_3': 'StrainMad_3',
              'zymo_d6300_micro_com_standard': 'Zymo_D6300',
              'zymo_d6310_micro_com_standard_ii': 'Zymo_D6310',
              'zymo_d6322_hmw_dna_standard': 'Zymo_D6322',
              'zymo_d6331_gut_micro_standard': 'Zymo_D6331'}


def main(input_path, output_plots):
    type_sequence = 'reads'
    status = Status('')
    # status.start()
    for organism in ['genus', 'species']:
        output_plots_organism = output_plots / organism
        if not output_plots_organism.exists():
            output_plots_organism.mkdir(parents=True)
        status.update(f'Calculating metrics for level {organism.upper()}')

        ####################
        # READ IN DATASETS #
        ####################

        # Read in all datasets it can find in the given folder.
        # Collect in a dictionary with structure dict[classifier][datasets] = {ground_truth, metric object}

        classifier_collection = defaultdict(dict)
        for classifier_file in input_path.glob('**/' + type_sequence + '/' + organism + '/output_*.tsv'):
            classifier_name = Path(classifier_file).stem.split('_')[-1]
            classifier_name = update_classifier_names[classifier_name]
            dataset_name = DATASET_ID[classifier_file.parts[5]]
            classifier_collection[classifier_name][dataset_name] = {}
            if classifier_name in list_abundance_classifiers:
                ground_truth_file = classifier_file.parents[3] / 'ground_truth_tax_updated.yml'
            else:
                ground_truth_file = classifier_file.parents[3] / 'ground_truth_updated.yml'
            with open(ground_truth_file, 'r') as f:
                theoretical_abundance = {'species': yaml.load(f, Loader=yaml.FullLoader)}
            theoretical_abundance['genus'] = {}
            for k, v in theoretical_abundance['species'].items():
                genus_name = k.split()[0]
                if genus_name in theoretical_abundance['genus']:
                    theoretical_abundance['genus'][genus_name] += v
                else:
                    theoretical_abundance['genus'][genus_name] = v
            classifier_collection[classifier_name][dataset_name].update({'ground_truth': pd.Series(theoretical_abundance[organism],
                                                                                                   name='ground_truth')})
            classifier_data = pd.read_table(
                classifier_file,
                usecols=[0, 1],
                index_col=0,
                names=[None, classifier_name],
                header=None
            ).squeeze("columns")

            classifier_collection[classifier_name][dataset_name].update(
                {'metric_object': ClassificationResults(classifier_name, classifier_data, list_uc)})

        with open(output_plots / 'used_classifiers.json', 'w') as f:
            json_out = {classifier_name: list(classifier_datasets) for classifier_name, classifier_datasets in
                        classifier_collection.items()}
            json.dump(json_out, f)
        with open(output_plots / 'used_datasets.txt', 'w') as f:
            txt_out = []
            for classifier_name, classifier_datasets in classifier_collection.items():
                for classifier_dataset_name in list(classifier_datasets):
                    txt_out.append(classifier_dataset_name)
            for item in np.unique(txt_out):
                f.write("%s\n" % item)

        logger.info(f'Read in files for {organism.upper()}')

        ##################
        # GATHER METRICS #
        ##################

        # Go through the collected dictionary and use the metric object to calculate metrics

        agg_metrics = defaultdict(dict)

        for classifier_name, datasets in classifier_collection.items():
            for dataset, data in datasets.items():
                classifier_result = data['metric_object']
                ground_truth = data['ground_truth']

                precision, recall, f1 = precision_recall_f1(classifier_result.only_classified(), ground_truth)
                agg_metrics[('precision', classifier_name)][dataset] = precision
                agg_metrics[('recall', classifier_name)][dataset] = recall
                agg_metrics[('F1', classifier_name)][dataset] = f1
                # for auc, the precision and recall over different thresholds is needed
                precision_w_thresh, recall_w_thresh, _, _ = precision_recall_f1_threshold(classifier_result.only_classified(True),
                                                                                          ground_truth)

                agg_metrics[('auprc', classifier_name)][dataset] = sklearn_metric.auc(recall_w_thresh,
                                                                                      precision_w_thresh)

                if ground_truth.sum() == 0:
                    l1 = np.nan
                else:
                    l1 = calculate_l1_distance(classifier_result.only_classified(True), ground_truth)
                agg_metrics[('L1', classifier_name)][dataset] = l1

        # df is multindex (metric, classifier) and datasets as columns
        agg_metrics_df = pd.DataFrame.from_dict(agg_metrics, orient='index').sort_index(level=0)

        # Give the multiindex a custom order
        custom_order = [(metric, classifier) for metric in ['precision', 'recall', 'F1', 'L1', 'auprc']
                        for classifier in classifier_order]
        agg_metrics_df = agg_metrics_df.loc[custom_order]
        # Add mean and median to dataframe, but only in output.
        # Otherwise, downstream analysis will throw an error
        agg_metrics_df.assign(mean=agg_metrics_df.mean(axis=1), median=agg_metrics_df.median(axis=1)). \
            to_csv(output_plots_organism / 'metric_table.tsv', sep='\t', na_rep="nan")

        logger.info('Gathered metrics')

        ###################
        # BOXPLOT METRICS #
        ###################

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        fig = plt.figure(figsize=(12, 15))
        ax = [fig.add_subplot(3, 2, 1), fig.add_subplot(3, 2, 2), fig.add_subplot(3, 2, 3), fig.add_subplot(3, 2, 4),
              fig.add_subplot(3, 2, 5), fig.add_subplot(3, 2, 6)]

        metric_names = agg_metrics_df.index.get_level_values(0).unique().to_list()

        for i, metric_name in enumerate(metric_names):
            # Create dataframe with dataset as index and classifier as column
            df = agg_metrics_df.loc[metric_name].T
            # boxplot for each column (so for each classifier)
            sns.boxplot(data=df,
                        color='white',
                        linewidth=0.5,
                        ax=ax[i]
                        )
            for patch in ax[i].patches:
                r, g1, b, a = patch.get_facecolor()
                patch.set_facecolor((r, g1, b, 0))
                patch.set_edgecolor('black')
            for line in ax[i].lines:
                line.set_color('black')

            for j, dataset in enumerate(df.to_numpy()):
                if np.isnan(dataset).all():
                    continue
                ax[i].scatter(x=[k for k in range(len(dataset))],
                              y=dataset,
                              label=df.index[j],
                              c=colors[j],
                              edgecolors='face',
                              linewidths=2,
                              s=5,
                              zorder=2
                              )

            ax[i].set_title(metric_name)
            # Coordinates that have the last x value do not get plotted because their
            # value is np.nan. This line adds the last x back
            # if metric_name == 'L1':
            #     ax[i].set_xlim(right=len(df.columns) - 1)

            for xtick in ax[i].get_xticklabels():
                xtick.set(rotation=40, ha='right', rotation_mode='anchor')
            # ax[i].set_xticklabels(df.columns, rotation=40)
            if metric_name == 'L1':
                ax[i].set_ylim(-0.05, 2.05)
            else:
                ax[i].set_ylim(-0.05, 1.05)

        handles, legends = ax[0].get_legend_handles_labels()

        ax[-1].legend(handles, legends, loc="center", borderaxespad=0,
                      fontsize='medium', numpoints=1, markerscale=1.5)
        ax[-1].set_axis_off()

        labels = [ALPH[i] + ')' for i in range(0, len(ax))]
        for i, _ in enumerate(ax[:-1]):
            # label physical distance in and down:
            trans = mtransforms.ScaledTranslation(5 / 72, 15 / 72, fig.dpi_scale_trans)
            ax[i].text(0.0, 1.0, labels[i], transform=ax[i].transAxes + trans,
                       fontsize='large', verticalalignment='top', fontfamily='serif', fontweight='heavy')

        plt.tight_layout()
        # plt.subplots_adjust(right=0.85)
        # plt.show()
        plt.savefig(output_plots_organism / 'boxplot_metrics.svg')

        logger.info('Plotted boxplot')

        ##################
        # CUTOFF METRICS #
        ##################

        metric_cutoff = {}
        METRIC = 'median'
        cutoff_range = range(int(0), int(0.03 * 10000), int(0.0005 * 10000))

        def aggregate_metric(array_of_metrics, METRIC):
            if METRIC == 'median':
                return np.nanmedian(array_of_metrics, axis=0)
            elif METRIC == 'mean':
                return np.nanmean(array_of_metrics, axis=0)

        for classifier_name, datasets in classifier_collection.items():
            precision_cutoff, recall_cutoff, F1_cutoff = [], [], []
            for dataset, data in datasets.items():
                classifier_result = data['metric_object']
                ground_truth = data['ground_truth']
                temp_precision, temp_recall, temp_f1 = [], [], []
                for cutoff in cutoff_range:
                    cutoff = cutoff / 10000
                    precision, recall, f1 = precision_recall_f1(classifier_result.only_classified(True), ground_truth, threshold=cutoff)
                    temp_precision.append(precision)
                    temp_recall.append(recall)
                    temp_f1.append(f1)
                precision_cutoff.append(temp_precision)
                recall_cutoff.append(temp_recall)
                F1_cutoff.append(temp_f1)

            metric_cutoff[classifier_name] = {'precision': precision_cutoff,
                                              'recall': recall_cutoff,
                                              'f1': F1_cutoff,
                                              'threshold': [i / 10000 for i in cutoff_range]
                                              }

        logger.info('Calculated cutoff metrics')

        #######################
        # METRICS CUTOFF PLOT #
        #######################

        y_label_dict = {'precision': r'Precision',
                        'recall': r'Recall',
                        'f1': 'F1'
                        }

        fig, axs = plt.subplots(1, 3, figsize=(15, 9))
        ax = axs.flatten()
        fig.suptitle(f'Aggregated {METRIC} of all datasets per classifier')
        for i, metric in enumerate(['precision', 'recall', 'f1']):
            for classifier in classifier_order:
                ax[i].plot(metric_cutoff[classifier]['threshold'],
                           aggregate_metric(np.array(metric_cutoff[classifier][metric]), METRIC),
                           label=classifier
                           )
                ax[i].set_ylabel(y_label_dict[metric])
                ax[i].xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))

        labels = [ALPH[i] + ')' for i in range(0, len(ax))]
        for i, _ in enumerate(ax):
            # label physical distance in and down:
            trans = mtransforms.ScaledTranslation(5 / 72, 15 / 72, fig.dpi_scale_trans)
            ax[i].text(0.0, 1.0, labels[i], transform=ax[i].transAxes + trans,
                       fontsize='large', verticalalignment='top', fontfamily='serif', fontweight='heavy')

        ax[2].legend()
        fig.supxlabel('Relative abundance')
        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'cutoff_metrics_plot.svg')

        ##########################
        # METRICS CUTOFF WITH Qs #
        ##########################

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        fig, axs = plt.subplots(3, len(classifier_collection), figsize=(22, 10), sharey='all', sharex='all')
        fig.suptitle(f'Aggregated {METRIC} of all datasets per classifier')
        # plt.xticks(rotation=55)
        for i, metric in enumerate(['precision', 'recall', 'f1']):
            for j, classifier in enumerate(classifier_order):
                axs[i, j].plot(metric_cutoff[classifier]['threshold'],
                               aggregate_metric(np.array(metric_cutoff[classifier][metric]), METRIC),
                               label=METRIC,
                               color=colors[j]
                               )

                axs[i, j].fill_between(x=metric_cutoff[classifier]['threshold'],
                                       y1=np.nanpercentile(metric_cutoff[classifier][metric], 25, axis=0),
                                       y2=np.nanpercentile(metric_cutoff[classifier][metric], 75, axis=0),
                                       alpha=0.5,
                                       color=colors[j],
                                       label='Q1 to Q3'
                                       )
                axs[i, j].fill_between(x=metric_cutoff[classifier]['threshold'],
                                       y1=np.nanpercentile(metric_cutoff[classifier][metric], 75, axis=0),
                                       y2=np.nanpercentile(metric_cutoff[classifier][metric], 100, axis=0),
                                       alpha=0.2,
                                       color=colors[j],
                                       label='MIN to MAX'
                                       )
                axs[i, j].fill_between(x=metric_cutoff[classifier]['threshold'],
                                       y1=np.nanpercentile(metric_cutoff[classifier][metric], 0, axis=0),
                                       y2=np.nanpercentile(metric_cutoff[classifier][metric], 25, axis=0),
                                       alpha=0.2,
                                       color=colors[j],

                                       )
                if j == 0:
                    axs[i, j].set_ylabel(metric)
                if i == 0:
                    axs[i, j].set_title(classifier)
                    axs[i, j].legend(loc='lower right')
                axs[i, j].xaxis.get_offset_text().set_visible(False)

        for ax in axs[-1, :]:
            ax.set_xticks([0, 0.001, 0.005, 0.008, 0.012])
            ax.set_xticklabels([0, 0.001, 0.005, 0.008, 0.012], rotation=65, rotation_mode='anchor', ha='right')
            ax.xaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))

        plt.xlim(left=0, right=0.0125)

        fig.supxlabel('Relative abundance threshold')
        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'cutoff_metrics_plot_Q.svg')

        # colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # fig, axs = plt.subplots(3 * len(classifier_collection), figsize=(14, 50), sharey='all', sharex='all')
        # fig.suptitle(f'Aggregated {METRIC} of all datasets per classifier')
        # for j, classifier in enumerate(metric_cutoff):
        #     for i, metric in enumerate(['precision', 'recall', 'f1']):
        #         axs[j * 3 + i].plot(metric_cutoff[classifier]['threshold'],
        #                             aggregate_metric(np.array(metric_cutoff[classifier][metric]), METRIC),
        #                             label=METRIC,
        #                             color=colors[j]
        #                             )
        #
        #         axs[j * 3 + i].fill_between(x=metric_cutoff[classifier]['threshold'],
        #                                     y1=np.nanpercentile(metric_cutoff[classifier][metric], 25, axis=0),
        #                                     y2=np.nanpercentile(metric_cutoff[classifier][metric], 75, axis=0),
        #                                     alpha=0.5,
        #                                     color=colors[j],
        #                                     label='Q2 to Q3'
        #                                     )
        #         axs[j * 3 + i].fill_between(x=metric_cutoff[classifier]['threshold'],
        #                                     y1=np.nanpercentile(metric_cutoff[classifier][metric], 75, axis=0),
        #                                     y2=np.nanpercentile(metric_cutoff[classifier][metric], 100, axis=0),
        #                                     alpha=0.2,
        #                                     color=colors[j],
        #                                     label='Q0 to Q4'
        #                                     )
        #         axs[j * 3 + i].fill_between(x=metric_cutoff[classifier]['threshold'],
        #                                     y1=np.nanpercentile(metric_cutoff[classifier][metric], 0, axis=0),
        #                                     y2=np.nanpercentile(metric_cutoff[classifier][metric], 25, axis=0),
        #                                     alpha=0.2,
        #                                     color=colors[j],
        #
        #                                     )
        #         if j % 3 == 0:
        #             axs[j * 3 + i].set_ylabel(metric)
        #         if i == 0:
        #             axs[j * 3 + i].set_title(classifier)
        #             # axs[i, j].legend(loc='lower right')
        # plt.tight_layout()
        # # plt.show()
        # plt.savefig(output_plots_organism / 'cutoff_metrics_plot_Q.svg')

        ##############
        # 2D PLOT PR #
        ##############

        fig, axs = plt.subplots(1, 3, figsize=(15, 8))
        ax = axs.flatten()

        error_bar_threshs = [0.0005, 0.001]

        for classifier in metric_cutoff:
            # IQR PLOT #

            # first value of each dataset is the metric when abundance cutoff is 0
            first_precisions = np.array(metric_cutoff[classifier]['precision'])[:, 0]
            first_recalls = np.array(metric_cutoff[classifier]['recall'])[:, 0]
            x_precision = aggregate_metric(np.array(metric_cutoff[classifier]['precision'])[:, 0], METRIC)
            y_recall = aggregate_metric(np.array(metric_cutoff[classifier]['recall'])[:, 0], METRIC)
            xerrs_iqr = np.array([[abs(x_precision - np.nanpercentile(first_precisions, 25))],
                                  [abs(x_precision - np.nanpercentile(first_precisions, 75))]])
            yerrs_iqr = np.array(
                [[abs(y_recall - np.nanpercentile(first_recalls, 25))],
                 [abs(y_recall - np.nanpercentile(first_recalls, 75))]])

            ax[0].errorbar(x_precision,
                           y_recall,
                           xerr=xerrs_iqr,
                           yerr=yerrs_iqr,
                           label=classifier,
                           fmt='o',
                           capsize=3, capthick=3
                           )

            # CUTOFF PLOT #

            for i, error_bar_thresh in enumerate(error_bar_threshs):
                # find index of threshold value
                index_thresh = np.where(np.array(metric_cutoff[classifier]['threshold']) == error_bar_thresh)[0][0]
                # take mean of precision and recall of every dataset at the above threshold
                thresh_mean_precision = aggregate_metric(np.array(metric_cutoff[classifier]['precision'])[:, index_thresh],
                                                         'median')
                thresh_mean_recall = aggregate_metric(np.array(metric_cutoff[classifier]['recall'])[:, index_thresh],
                                                      'median')

                if x_precision - thresh_mean_precision < 0:
                    xerrs_5 = np.array([[0], [abs(x_precision - thresh_mean_precision)]])
                else:
                    xerrs_5 = np.array([[abs(x_precision - thresh_mean_precision)], [0]])
                if y_recall - thresh_mean_recall < 0:
                    yerrs_5 = np.array([[0], [abs(y_recall - thresh_mean_recall)]])
                else:
                    yerrs_5 = np.array([[abs(y_recall - thresh_mean_recall)], [0]])

                ax[i + 1].errorbar(x_precision,
                                   y_recall,
                                   xerr=xerrs_5,
                                   yerr=yerrs_5,
                                   label=classifier,
                                   fmt='o',
                                   capsize=3, capthick=3
                                   )

        fig.suptitle(f'{METRIC} value of PR for all datasets per classifier')
        fig.supxlabel(r'Precision')
        fig.supylabel(r'Recall')
        ax[0].set_title('Error bars represent IQR values', fontdict={'fontsize': 9})
        ax[1].set_title(f'Error bars represent PR values at {error_bar_threshs[0]*100}% relative abundance cutoff',
                        fontdict={'fontsize': 8})
        ax[2].set_title(f'Error bars represent PR values at {error_bar_threshs[1]*100}% relative abundance cutoff',
                        fontdict={'fontsize': 8})

        for ax in axs:
            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)

        labels = [ALPH[i] + ')' for i in range(0, len(axs))]
        for i, _ in enumerate(axs):
            # label physical distance in and down:
            trans = mtransforms.ScaledTranslation(5 / 72, 15 / 72, fig.dpi_scale_trans)
            axs[i].text(0.0, 1.0, labels[i], transform=axs[i].transAxes + trans,
                        fontsize='large', verticalalignment='top', fontfamily='serif', fontweight='heavy')

        plt.legend(loc='lower right')
        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'PR_errorbars.svg')

        logger.info('Plotted the rest')

    status.stop()


if __name__ == "__main__":
    main(Path('/scratch/alvanuffelen/hackaton/samples/SRR9198630/'), Path('/home/alvanuffelen/Downloads/test'))
