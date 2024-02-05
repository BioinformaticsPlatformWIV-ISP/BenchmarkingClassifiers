from .analysis_function import ClassificationResults, plot_top15, plot_classification_bar, \
    plot_ecdf, plot_ecdf_zoomed, plot_tpr_fps, plot_precision_recall, generate_l1_matrix, plot_l1_heatmap, tp_fp_fn, \
    precision_recall_f1

import json
import logging
import math
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from rich.status import Status
import yaml

logger = logging.getLogger('amr_report')

mpl.rcParams["mathtext.default"] = 'regular'
plt.style.use('bmh')
plt.rcParams['axes.axisbelow'] = True

# Number of rows if there are multiple subplots
ROWS = 2


# Number of cols based on the total number of classifiers and rows
def n_rows(total, rows):
    return int(math.ceil(total / rows))


list_uc = ['no_hit',
           'ambigious',
           'unclassified',
           'missing_genus'
           ]

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

def main(data, output_plots, allowed_classifiers):
    ########################
    # READ IN GROUND TRUTH #
    ########################

    theoretical_abundances = {}
    for ground_truth_file in data.parents[1].glob('ground_truth*_updated.yml'):
        with open(ground_truth_file, 'r') as f:
            theoretical_abundance = {'species': yaml.load(f, Loader=yaml.FullLoader)}
        theoretical_abundance['genus'] = {}

        for k, v in theoretical_abundance['species'].items():
            genus_name = k.split()[0]
            if genus_name in theoretical_abundance['genus']:
                theoretical_abundance['genus'][genus_name] += v
            else:
                theoretical_abundance['genus'][genus_name] = v

        if ground_truth_file.stem.split('_')[-2] == 'tax':
            theoretical_abundances['tax'] = theoretical_abundance
        else:
            theoretical_abundances['sequence'] = theoretical_abundance

    logger.info('Read in ground truth')

    status = Status('')
    # status.start()
    for organism in ['genus', 'species']:
        status.update(f'Processing {organism}')
        ground_truth = pd.Series(theoretical_abundances['sequence'][organism], name='ground_truth')
        ground_truth_tax = pd.Series(theoretical_abundances['tax'][organism], name='ground_truth_tax')
        # Input folder
        data_organism = data / organism
        # Output path for plots
        output_plots_organism = output_plots / organism
        if not output_plots_organism.exists():
            output_plots_organism.mkdir(parents=True)

        #################################
        # READ IN DATA FROM CLASSIFIERS #
        #################################

        logger.info(f'Reading in {organism} files')

        classifiers_input = []

        for input_classifier in sorted(data.glob(organism + '/' + '*.tsv')):
            classifier_name = input_classifier.stem.split(sep='_')[-1]
            classifier_name = update_classifier_names[classifier_name]
            if allowed_classifiers == 'all' or classifier_name in allowed_classifiers:
                classifier_table = pd.read_table(
                    data_organism / input_classifier.name,
                    usecols=[0, 1],
                    index_col=0,
                    names=[None, classifier_name],
                    header=None
                ).squeeze("columns")
                classifiers_input.append(ClassificationResults(classifier_name, classifier_table, list_uc))
            logger.debug(f'Read in {input_classifier}')

        ################
        # METRIC TABLE #
        ################
        # BOTH READ AND ABUNDANCE CLASSIFIERS #
        collected_metrics = {}
        for classifier_result in classifiers_input:
            collected_metrics[classifier_result.classifier] = {**tp_fp_fn(classifier_result.only_classified(), ground_truth),
                                                               **precision_recall_f1(classifier_result.only_classified(), ground_truth)}
        top_table = pd.DataFrame.from_dict(collected_metrics, orient='index').T
        top_table.to_csv(output_plots_organism / 'metric_table.tsv', sep='\t')

        ########
        # FNs? #
        ########

        FNs = {}
        for classifier_result in classifiers_input:
            FNs[classifier_result.classifier] = sorted([i for i in ground_truth.index.to_list()
                                                        if i not in classifier_result.only_classified().index.to_list()])

        json.dump(FNs,
                  (output_plots_organism / "FNs.json").open('w'),
                  indent=4)

        ###############################
        # TABLE (UN)CLASSFIED NUMBERS #
        ###############################
        # ONLY READ CLASSIFIERS #

        un_classified = {}
        for classifier_result in classifiers_input:
            if classifier_result.classifier not in list_abundance_classifiers:
                un_classified[classifier_result.classifier] = classifier_result.bundled_classified(False)
        un_classified_values = pd.DataFrame.from_dict(un_classified)
        un_classified_values.to_csv(output_plots_organism / 'number_table.tsv', sep='\t')

        #############################
        # BAR TOP 15 PER CLASSIFIER #
        #############################
        # ONLY READ CLASSIFIERS #
        gs = GridSpec(3, 3)
        fig = plt.figure(figsize=(18, 9))

        fig.suptitle(f'Top 15 of each classifier')
        for i, classifier_results in enumerate(classifiers_input):
            if classifier_results.classifier not in list_abundance_classifiers:
                ax = fig.add_subplot(gs[i])
                plot_top15(classifier_results.results, ground_truth, classifier_results.classifier, ax)

        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'top15.svg')

        #################################
        # BAR CLASSIFICATION PROPORTION #
        #################################
        # BOTH READ AND ABUNDANCE CLASSIFIERS #

        fig, ax = plt.subplots(figsize=(9, 6))

        if len(ground_truth.index) <= 20:
            classifier_bundled_misclassified_others = []
            for classifier in classifiers_input:
                classifier_bundled_misclassified_others.append(
                    classifier.bundled_misclassified_others(ground_truth, True))

            plot_classification_bar(classifier_bundled_misclassified_others, ground_truth, ground_truth_tax, ax)
        else:
            ax.text(0.5, 0.5, s=f'Too much organisms in the ground truth ({len(ground_truth.index)})', weight='bold', va='center',
                    ha='center', color='r')
        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'classification_bar_w_uc.svg')

        #################################################
        # BAR CLASSIFICATION PROPORTION W/O UNCLASSIFIED #
        #################################################
        # BOTH READ AND ABUNDANCE CLASSIFIERS #

        fig, ax = plt.subplots(figsize=(9, 6))

        classifier_bundled_misclassified_excl_others = []
        for classifier in classifiers_input:
            classifier_bundled_misclassified_excl_others.append(
                classifier.bundled_misclassified_excl_others(ground_truth, True))
        pd.DataFrame(classifier_bundled_misclassified_excl_others).to_csv(output_plots_organism / 'classification_bar_numbers.tsv',
                                                                          sep='\t')
        if len(ground_truth.index) <= 20:
            plot_classification_bar(classifier_bundled_misclassified_excl_others, ground_truth, ground_truth_tax, ax)
        else:
            ax.text(0.5, 0.5, s=f'Too much organisms in the ground truth ({len(ground_truth.index)})', weight='bold', va='center',
                    ha='center', color='r')
        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'classification_bar.svg')

        #############
        # ECDF PLOT #
        #############
        # BOTH READ AND ABUNDANCE CLASSIFIERS #

        fig = plt.figure(figsize=(9, 8))
        # Create gridspec
        G = GridSpec(len(classifiers_input), ROWS)
        # Add subplots from gridspec and make list
        ax = [fig.add_subplot(G[:, 0])]
        ax.extend([fig.add_subplot(G[i, 1], sharex=fig.axes[0]) for i in range(len(classifiers_input))])

        classifier_results = []
        for classifier_result in classifiers_input:
            classifier_results.append(classifier_result.only_classified(True))

        plot_ecdf(classifier_results, ax[0])

        for i, classifier_results in enumerate(classifiers_input):
            plot_ecdf_zoomed(classifier_results.only_classified(True),
                             plt.rcParams['axes.prop_cycle'].by_key()['color'][i], ax[i + 1])

        fig.supylabel('Count')
        fig.supxlabel('Relative abundance of reads')
        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'ecdf.svg')

        ###########
        # TPR FPS #
        ###########
        # BOTH READ AND ABUNDANCE CLASSIFIERS #

        # noinspection PyTypeChecker
        gs = GridSpec(ROWS, n_rows(len(classifiers_input), ROWS))
        fig = plt.figure(figsize=(16, 10))
        for i, classifier_result in enumerate(classifiers_input):
            ax = fig.add_subplot(gs[i])
            plot_tpr_fps(classifier_result.only_classified(True), ground_truth, classifier_result.classifier, ax)

        # set labels
        fig.supylabel(r'TPR $\left(\frac{TP}{TP + FN}\right)$')
        fig.supxlabel(r'FPs')
        # plt.show()
        # plt.tight_layout()
        plt.savefig(output_plots_organism / 'tpr_fps_curve.svg')

        ######
        # PR #
        ######
        # BOTH READ AND ABUNDANCE CLASSIFIERS #

        gs = GridSpec(ROWS, n_rows(len(classifiers_input), ROWS))
        fig = plt.figure(figsize=(15, 9))
        fig.suptitle(f'PR Curves')
        for i, classifier_result in enumerate(classifiers_input):
            ax = fig.add_subplot(gs[i])
            plot_precision_recall(classifier_result.only_classified(True), ground_truth, ax)
            ax.set_title(classifier_result.classifier)
        # set labels
        fig.supylabel(r'Precision $\left(\frac{TP}{TP + FP}\right)$')
        fig.supxlabel(r'Recall $\left(\frac{TP}{TP + FN}\right)$')
        # fig.subplots_adjust(wspace=2, hspace=2)
        # plt.show()
        plt.tight_layout()
        plt.savefig(output_plots_organism / 'pr_curve.svg')

        ######
        # L1 #
        ######
        # ONLY READ CLASSIFIERS #

        fig, ax = plt.subplots(figsize=(9, 6))
        if sum(ground_truth.values) != 0:
            classifier_results = [ground_truth, ground_truth_tax]
            for classifier_result in classifiers_input:
                classifier_results.append(classifier_result.only_classified(True))
            l1_matrix = generate_l1_matrix(classifier_results)
            plot_l1_heatmap(l1_matrix, ax)

        else:
            ax.text(0.5, 0.5, s='Unknown ground truth (abundance)', weight='bold', va='center', ha='center', color='r')

        plt.tight_layout()
        # plt.show()
        plt.savefig(output_plots_organism / 'l1_heatmap.svg')

    # status.stop()
