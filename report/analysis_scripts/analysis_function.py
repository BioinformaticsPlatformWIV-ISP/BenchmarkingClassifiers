import itertools
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from dataclasses import dataclass
from sklearn import metrics
from typing import List

COLOR_LIST = [
    '#1F77B3', '#ADBED9', '#FF7F0E', '#FFBA78', '#2CB05C',
    '#98DF8A', '#D62728', '#FF9960', '#9566BD', '#C69CDA',
    '#8C564B', '#C17C94', '#E377C2', '#FDBDF6', '#7F7F7F',
    '#CFCFCF', '#B9BD22', '#DBDB8E', '#17BECF', '#9ECBF3'
]

list_uc = ['no_hit',
           'ambigious',
           'unclassified',
           'missing_genus',
           'not_distributed',
           'rounding_error'
           ]


@dataclass
class ClassificationResults:
    classifier: str
    results: pd.Series
    names_unclassified: List[str]

    @staticmethod
    def _normalize_series(series):
        return series / series.sum()

    def excl_no_hit(self, normalized=False):
        results_excl_no_hit = self.results[self.results.index != 'no_hit']
        if normalized:
            return self._normalize_series(results_excl_no_hit)
        return results_excl_no_hit

    def only_classified(self, normalized=False):
        """
        Return the reads not in the list_uc.
        """
        only_classified = self.results[~self.results.index.isin(list_uc)]
        if normalized:
            return self._normalize_series(only_classified)
        return only_classified

    def bundled_misclassified_others(self, ground_truth, normalized=False):
        """
        'Others' --> list_uc
        'Misclassified' --> Not in the ground truth
        The rest are correctly classified and is as-is
        """
        groups = np.where(np.isin(self.results.index, list_uc), 'Others',
                          np.where(np.isin(self.results.index, ground_truth.index), self.results.index,
                                   'Misclassified'))
        groups_cat = pd.Categorical(groups,
                                    categories=list(ground_truth.index) + ['Misclassified', 'Others'])
        bundled_misclassified_others = self.results.groupby(groups_cat, observed=False).aggregate('sum')
        if normalized:
            return self._normalize_series(bundled_misclassified_others)
        return bundled_misclassified_others

    def bundled_misclassified_excl_others(self, ground_truth, normalized=False):
        """
        Return the reads not in the list_uc BUT bundle reads in the ground truth as 'Misclassified'
        """
        only_classified = self.only_classified(normalized)
        groups = np.where(np.isin(only_classified.index, ground_truth.index), only_classified.index, 'Misclassified')
        groups_cat = pd.Categorical(groups,
                                    categories=list(ground_truth.index) + ['Misclassified'])
        bundled_misclassified_excl_others = only_classified.groupby(groups_cat, observed=False).aggregate('sum')
        if normalized:
            return self._normalize_series(bundled_misclassified_excl_others)
        return bundled_misclassified_excl_others

    def bundled_classified(self, normalized=False):
        """
        'Classified' --> ground truth + classified but not in ground truth
        The rest are in list_uc and are as-is
        """
        groups = np.where(np.isin(self.results.index, list_uc), self.results.index, 'Classified')
        groups_cat = pd.Categorical(groups,
                                    categories=['Classified'] + list_uc)
        bundled_classified = self.results.groupby(groups_cat, observed=False).aggregate('sum')
        if normalized:
            return self._normalize_series(bundled_classified)
        return bundled_classified


def plot_classification_bar(multiple_results: List[pd.Series], ground_truth, ground_truth_tax, ax=None):
    ax = ax
    theoretical = ground_truth.copy(deep=True).rename('Theoretical')
    theoretical_tax = ground_truth_tax.copy(deep=True).rename('Theoretical taxonomy')
    theoretical['Misclassified'], theoretical_tax['Misclassified'] = 0, 0
    if 'Others' in multiple_results[0].index:
        colormap = mpl.colors.ListedColormap(COLOR_LIST[:len(multiple_results[0]) - 2] + ['k', 'gray'])
        theoretical['Others'], theoretical_tax['Others'] = 0, 0
    else:
        colormap = mpl.colors.ListedColormap(COLOR_LIST[:len(multiple_results[0]) - 1] + ['k'])
    pd.concat([theoretical, theoretical_tax, *multiple_results], axis=1).fillna(0).T.plot.bar(stacked=True,
                                                                                              rot=0,
                                                                                              ylabel='Relative Abundance',
                                                                                              ax=ax,
                                                                                              colormap=colormap
                                                                                              )
    ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right', va='top', rotation_mode='anchor')


def tp_fp_fn(results: pd.Series, ground_truth: pd.Series, threshold=None) -> pd.Series:
    """
    Calculates the true positives, false positives and false negatives.
    :return: dictionary of tps, fps and fns
    """
    if threshold is None:
        classifier_status = np.isin(results.index, ground_truth.index)
    else:
        # Relative abundance of each organism
        class_reads_above_thresh = results[results >= threshold].index.values
        classifier_status = np.isin(class_reads_above_thresh, ground_truth.index)

    tp, fp = sum(classifier_status), sum(~classifier_status)
    fn = len(ground_truth) - tp
    return pd.Series(dict(TP=tp, FP=fp, FN=fn))


def precision_recall_f1(results: pd.Series, ground_truth: pd.Series, threshold=None) -> pd.Series:
    """
    Calculates the precision, recall and f1
    :return: dictionary of precision, recall and f1
    """
    tp, fp, fn = tp_fp_fn(results, ground_truth, threshold=threshold)
    if tp + fp != 0:
        precision = tp / (tp + fp)
    else:
        precision = np.nan
    recall = tp / (tp + fn)
    if precision != np.nan and recall != 0:
        f1 = 2 / (precision ** -1 + recall ** -1)
    else:
        f1 = np.nan
    return pd.Series(dict(Precision=precision, Recall=recall, F1=f1))


def tps_fps(results: pd.Series, ground_truth: pd.Series) -> pd.Series:
    # Relative abundance of each organism
    y_score = results.to_numpy()
    # Make boolean array for y_score (1 if in ground truth, 0 if not)
    y_true = np.where(results.index.isin(ground_truth.index), 1, 0)
    # Indices of the sorted score (high to low)
    desc_score_indices = np.argsort(y_score, kind="mergesort")[::-1]
    # Sort the score and the truth
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]

    # Get indices of sorted score when there is a difference between the next neighbour
    distinct_value_indices = np.where(np.diff(y_score))[0]

    # We add 'y_true.size - 1' because we also want the values past the lowest threshold (the last index)
    threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]

    # true positive is y_score >= threshold, false positive if y_score < threshold
    tps = np.cumsum(y_true)[threshold_idxs]
    # With every change in the y_score, a new TP of FP is introduced. The number of changes is denoted by
    # 1+threshold_idxs, which is equal to tps + fps.
    fps = 1 + threshold_idxs - tps
    threshold = y_score[threshold_idxs]

    return pd.Series(dict(TPS=tps, FPS=fps, threshold=threshold))


def tpr_fps(results: pd.Series, ground_truth: pd.Series) -> pd.Series:
    tps, fps, threshold = tps_fps(results, ground_truth)
    fns = len(ground_truth) - tps
    tpr = tps / (tps + fns)

    return pd.Series(dict(TPR=tpr, FPS=fps, threshold=threshold))


def precision_recall_f1_threshold(results: pd.Series, ground_truth: pd.Series) -> pd.Series:
    tps, fps, threshold = tps_fps(results, ground_truth)
    fns = len(ground_truth) - tps
    precision = tps / (tps + fps)
    recall = tps / (tps + fns)
    with np.errstate(divide='ignore'):
        f1 = 2 / (precision ** -1 + recall ** -1)

    # reverse the outputs so recall is decreasing
    return pd.Series(dict(precision=np.r_[precision[::-1], 1], recall=np.r_[recall[::-1], 0],
                          threshold=threshold[::-1], f1=np.r_[f1[::-1], 0]))


def calculate_l1_distance(results1: pd.Series, results2: pd.Series) -> float:
    """
    Calculates the l1 (Manhattan) distance between organisms in N-dimension.
    :return: l1 distance
    """
    union_organisms = set().union(results1.index, results2.index)
    # Extend both lists with 0 for missing organisms
    reads_count_extended = {organism: results1.get(organism, 0) for organism in union_organisms}
    second_count_extended = {organism: results2.get(organism, 0) for organism in union_organisms}
    l1_distance = sum(abs(val1 - val2) for
                      val1, val2 in zip(reads_count_extended.values(), second_count_extended.values()))
    return l1_distance


def plot_top15(results: pd.Series, ground_truth: pd.Series, classifier_name: str, ax=None):
    ax = ax
    top15 = results.sort_values(ascending=False).iloc[:15]
    sns.barplot(x=top15.values,
                y=top15.index,
                ax=ax,
                color='grey'
                )
    ax.set_xlabel('Count')
    # ax.set_ylabel('')
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    ax.set_title(classifier_name)
    # Make organisms bold if they are in ground truth
    for i, label in enumerate(ax.yaxis.get_ticklabels()):
        if label.get_text() in ground_truth.index:
            label.set_fontweight('bold')


def plot_ecdf(multiple_results: List[pd.Series], ax=None):
    ax = ax
    for results in multiple_results:
        sns.ecdfplot(results,
                     stat='count',
                     label=results.name,
                     ax=ax
                     )

    ax.legend()
    ax.set_yscale('log')
    ax.set_ylabel(None)
    ax.set_xlabel(None)
    y_last = ax.get_ylim()[-1]
    ax.set_ylim(top=y_last + 10 ** int(math.log10(y_last)))


def plot_ecdf_zoomed(results: pd.Series, color, ax=None):
    ax = ax
    sns.ecdfplot(results,
                 stat='count',
                 label=None,
                 ax=ax,
                 color=color
                 )
    x_data = ax.get_lines()[0].get_data()[0]
    y_data = ax.get_lines()[0].get_data()[-1]
    # Determine at which y value the zoomed plot will start (when the difference of two x values is higher than 0.001)
    y_start = (np.diff(x_data)[1:] > 0.001).argmax()
    ax.set_ylim(bottom=y_data[y_start])

    # y_first will equal y_start (see above)
    y_first, y_last = ax.get_ylim()[0], ax.get_ylim()[-1]
    # set top limit so all zoomed plots have an equal amount of space at the top
    ax.set_ylim(top=y_last + 0.1 * (y_last - y_first))

    ax.set_ylabel(None)
    ax.set_xlabel(None)


def plot_tpr_fps(results: pd.Series, ground_truth: pd.Series, classifier_name: str, ax=None):
    ax = ax
    tpr, fps, threshold = tpr_fps(results, ground_truth)
    if np.sum(tpr) == 0:
        ax.text(0.5, 0.5, 'No TPs', ha='center', va='center', color='red')
    else:
        # get index where max TPR is reached
        # only plot to max TPR (no use to plot further)
        y_max_ind = tpr.argmax()
        ax.plot(fps[:y_max_ind + 1],
                tpr[:y_max_ind + 1],
                )
        # plot the point where max TPR is reached as red dot
        ax.scatter(fps[y_max_ind],
                   tpr[y_max_ind],
                   c='red',
                   label=f'{threshold[y_max_ind]:09.7f}',
                   s=4 ** 2,
                   zorder=2.5
                   )
        # add red horizontal line to red dot
        ax.axhline(y=tpr[y_max_ind],
                   color='red',
                   ls='dashed',
                   lw=0.75
                   )
        # If TPR makes a jump, annotate that point
        index_wo_annot = []
        for i in range(len(threshold[:y_max_ind])):
            if i == 0 or tpr[i] - tpr[i - 1]:
                ax.annotate(f'{threshold[i]:09.7f}',
                            (fps[i], tpr[i]),
                            xytext=(0, 2),
                            textcoords='offset points',
                            ha='center', va='bottom',
                            fontsize='x-small'
                            )
            else:
                index_wo_annot.append(i)
        # color points without annotations grey
        ax.scatter(fps[index_wo_annot],
                   tpr[index_wo_annot],
                   c='grey',
                   s=4 ** 2,
                   zorder=2.5
                   )
        # color points with annotation blue
        index_w_annot = [i for i in range(len(threshold[:y_max_ind])) if i not in index_wo_annot]
        ax.scatter(fps[index_w_annot],
                   tpr[index_w_annot],
                   c='blue',
                   s=4 ** 2,
                   zorder=2.5
                   )
        ax.legend(title='Cutoff threshold', loc='lower right', scatteryoffsets=[0.6],
                  fontsize='x-small', title_fontsize='small')
    ax.set_ylim(0, 1.05)
    ax.margins(x=0.3)
    ax.set_title(classifier_name)


def plot_precision_recall(results: pd.Series, ground_truth: pd.Series, ax=None):
    ax = ax
    precision, recall, threshold, f1 = precision_recall_f1_threshold(results, ground_truth)
    if np.sum(precision) == 0:
        ax.text(0.5, 0.5, 'No TPs', ha='center', va='center', color='red')
    else:
        ax.plot(recall,
                precision
                )
        # ax.plot(recall,
        #         f1)
        # Threshold has one value less than PR, because the endpoint for PR is (0, 1) with no threshold!
        # reminder: threshold is sorted from small to large, so P from s to l and vice versa for R
        index_w_annot = []
        j = 1
        for i in range(len(threshold)):
            # last point
            if i == len(threshold) - 1:
                ax.annotate(f'{threshold[i]:09.7f}',
                            (recall[i], precision[i]),
                            xytext=(0, 2),
                            textcoords='offset points',
                            ha='center', va='bottom',
                            fontsize='x-small',
                            color='red'
                            )
                index_w_annot.append(i)
            # first point
            elif i == 0:
                ax.annotate(f'{threshold[i]:09.7f}',
                            (recall[i], precision[i]),
                            xytext=(0, -4),
                            textcoords='offset points',
                            ha='center', va='top',
                            fontsize='x-small',
                            color='red'
                            )
                index_w_annot.append(i)

            # do not annotate if jump in recall is not significant (both next and previous)
            elif abs(recall[i + 1] - recall[i]) < 0.01 and abs(recall[i] - recall[i - 1]) < 0.01:
                continue
            # do not annotate if jump in recall is not significant (both next and previous)
            elif abs(precision[i + 1] - precision[i]) < 0.01 and abs(precision[i] - precision[i - 1]) < 0.01:
                continue
            # do not annotate if precision difference with last annotated point is less than 0.03
            elif abs(precision[i] - precision[index_w_annot[-1]]) < 0.03:
                continue
            else:
                # annotate with number
                ax.annotate(j,
                            (recall[i], precision[i]),
                            xytext=(0, 2),
                            textcoords='offset points',
                            fontsize='x-small',
                            ha='center', va='bottom')
                # increase annotation number
                j += 1
                # keep track of annotated index of PR
                index_w_annot.append(i)

        # Create string with track number (j) and threshold
        text_str = '\n'.join(f'{a:d}' + '. ' + f'{b:09.7f}' for a, b in
                             zip(range(1, len(index_w_annot) - 1), threshold[index_w_annot[1:-1]]))
        # print string next to plot
        ax.text(1.1, 0.5, r"$\bf{Threshold}$" + '\n' + text_str,
                va='center', ma='right',
                fontsize='x-small',
                bbox=dict(boxstyle='round', facecolor='none', edgecolor='black'))

        # mark points with annotation blue
        ax.scatter(recall[index_w_annot[1:-1]],
                   precision[index_w_annot[1:-1]],
                   c='blue',
                   s=4 ** 2,
                   zorder=2.6
                   )
        # mark first and last point red.
        # last point is -2, because -1 is the (0,1) point
        ax.scatter([recall[0], recall[-2]],
                   [precision[0], precision[-2]],
                   c='red',
                   s=4 ** 2,
                   zorder=2.6
                   )

        # create list of indices without annotation and mark them frey
        index_wo_annote = [i for i in range(len(threshold)) if i not in index_w_annot]
        ax.scatter(recall[index_wo_annote],
                   precision[index_wo_annote],
                   c='grey',
                   s=4 ** 2,
                   zorder=2.5
                   )

        ax.text(0.05, 0.05, f'AUC = {metrics.auc(recall, precision):.5f}',
                fontsize='x-small',
                ha="left", va="bottom",
                bbox=dict(fc='none', edgecolor='black', boxstyle='round,pad=0.3'))

    ax.set_xlim(0 - 0.05, 1 + 0.05)
    ax.set_ylim(0 - 0.05, 1 + 0.05)


def generate_l1_matrix(multiple_results: List[pd.Series]):
    all_l1_distances = []
    for results_pair in itertools.combinations_with_replacement(multiple_results, 2):
        if results_pair[0].equals(results_pair[1]):
            all_l1_distances.append(0)
        else:
            all_l1_distances.append(calculate_l1_distance(results_pair[0], results_pair[1]))
    l1_matrix = np.zeros((len(multiple_results), len(multiple_results)))
    l1_matrix[np.triu_indices(len(multiple_results))] = all_l1_distances
    labels = []
    for result in multiple_results:
        labels.append(result.name)
    return pd.DataFrame(l1_matrix, columns=labels, index=labels)


def plot_l1_heatmap(l1_matrix, ax=None):
    ax = ax
    mask = np.ones_like(l1_matrix)
    mask[np.triu_indices_from(mask)] = False
    sns.heatmap(l1_matrix,
                vmin=0,
                vmax=2,
                mask=mask,
                annot=True)

    # Rotate the tick labels and set their alignment.
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right",
                       rotation_mode="anchor")

    # Rotate the tick labels and set their alignment.
    ax.set_yticklabels(ax.get_yticklabels(), rotation=45, ha="right",
                       rotation_mode="anchor")
