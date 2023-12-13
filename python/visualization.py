import pandas as pd
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import *
import matplotlib as mpl
from pandas.api.types import CategoricalDtype
from sklearn.metrics import roc_curve, auc
from sklearn import metrics
from utils import get_prior, cal_pr_values, cal_metrics_visualization

# set plot style
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.style.use('ggplot')
sns.set_style('white')
tick_font = 20
label_font = 24
legend_font = 18

model_name_list = ['ClinPred', 'REVEL', 'VARITY', 'MetaSVM', 'MetaLR',
                   'VEST4', 'M-CAP', 'MutPred', 'PrimateAI', 'MutationAssessor',
                   'LIST-S2', 'SIFT4G', 'DANN', 'MutationTaster', 'MAGPIE']
model_list = ['ClinPred_score', 'REVEL_score', 'VARITY_R', 'MetaSVM_score', 'MetaLR_score',
              'VEST4_score', 'M-CAP_score', 'MutPred_score', 'PrimateAI_score', 'MutationAssessor_score',
              'LIST-S2_score', 'SIFT4G_score', 'DANN_rankscore', 'MutationTaster_score', 'MAGPIE_pred']
compare_list = ['func', 'CLASS'] + model_list
threshold_list = [0.5, 0.5, 0.5, 0.5, 0.5,
                  0.5, 0.025, 0.79, 0.8, 0.8,
                  0.85, 0.95, 0.5, 0.3, 0.5]
colors = ['rosybrown', 'Teal', 'steelblue', 'DarkSlateGrey', 'OliveDrab',
          'DarkSeaGreen', 'CadetBlue', 'CornflowerBlue', 'peru', 'DarkKhaki',
          'GoldenRod', 'DarkGrey', 'LightBlue', 'royalblue', 'crimson']
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
save_dir = os.path.join(root, 'data/output/visualization/')
cutoff = None
version = 'old'


# ROC plot
def multi_models_roc(names, colors, data, save=None, dpin=500):
    """
    将多个机器模型的roc图输出到一张图上

    Args:
        names: list, 多个模型的名称
        colors: list, 多个模型的颜色
        data: DataFrame, 用于画图的数据
        save: 选择是否将结果保存（默认为png格式）, 如果有值就作为文件名

    Returns:
        返回图片对象plt
    """
    #         plt.figure(figsize=(12, 10), dpi=dpin)
    plt.figure(figsize=(10, 10), dpi=dpin)
    plot_list = []
    label_list = []

    sns.set(context=None,
            style='white',
            palette='deep',
            font='Times New Roman',
            font_scale=1.3,
            color_codes=True,
            rc=None)

    order_dict = {}
    order_index = 0
    for (name, colorname) in zip(names, colors):
        df = data.loc[:, [name, 'CLASS']].dropna().astype('float64')
        fpr, tpr, thresholds = roc_curve(df.CLASS, df[name].to_list(), pos_label=1)
        label = '{} ({:.3f})'.format('_'.join(name.split('_')[:-1]), auc(fpr, tpr))
        order_dict[auc(fpr, tpr)] = order_index
        order_index += 1
        label_list.append(label)
        plt_temp, = plt.plot(fpr, tpr, lw=5, label=label, color=colorname)
        plot_list.append(plt_temp)
        #             plt.axis('square')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.xticks(fontsize=tick_font)
        plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0], fontsize=tick_font)
        plt.xlabel('False Positive Rate', fontsize=label_font)
        plt.ylabel('True Positive Rate', fontsize=label_font)
        #             plt.title('ROC Curve',fontsize=20)
    plt.plot([0, 1], [0, 1], '--', lw=5, color='black')
    # plt.legend(handles = list(reversed(plot_list)), labels = list(reversed(label_list)),
    #            loc='lower right',fontsize=18)
    order_list = pd.DataFrame.from_dict({'AUC': order_dict.keys(), 'order': order_dict.values()}).sort_values(by='AUC',
                                                                                                              ascending=False).order.tolist()
    plt.legend(handles=[plot_list[idx] for idx in order_list],
               labels=[label_list[idx] for idx in order_list],
               loc='lower right', fontsize=legend_font)

    if save is not None:
        plt.savefig(save, dpi=500, bbox_inches='tight', format='pdf')

    return plt


def multi_models_pr(names, colors, data, mode, save=None, dpin=500, loc='upper left'):
    sns.set(context=None,
            style='white',
            palette='deep',
            font='Times New Roman',
            font_scale=1.3,
            color_codes=True,
            rc=None)

    test_precision, test_recall = [0.9, 0.9]
    plt.figure(figsize=(10, 10), dpi=dpin)
    label_list = []
    plot_list = []
    order_dict = {}
    order_index = 0
    for (name, colorname) in zip(names, colors):

        df = data.loc[:, [name, 'CLASS']].dropna().astype('float64')
        if 'SIFT4G' in name:
            df.SIFT4G_score = [1 - i for i in df.SIFT4G_score.tolist()]
        #             fpr, tpr, thresholds = precision_recall_curve(df.CLASS, df[name + '_pred'].to_list(), pos_label=1)
        y = df.CLASS.tolist()
        y_predicted = df[name].to_list()
        prior = get_prior(y)
        precisions, recalls, prc_thresholds = metrics.precision_recall_curve(y, y_predicted)
        recalls = np.insert(recalls, 0, 1)
        precisions = np.insert(precisions, 0, prior)
        # balanced precision and recall (prior = 0.5)
        balanced_precisions = precisions * (1 - prior) / (precisions * (1 - prior) + (1 - precisions) * prior)
        balanced_recalls = recalls
        if mode == 'prc':
            [auprc, up_auprc, rfp, pfr] = cal_pr_values(precisions, recalls, test_precision, test_recall)
            label = '{} ({:.3f})'.format('_'.join(name.split('_')[:-1]), auprc)
            order_dict[auprc] = order_index
            order_index += 1
            label_list.append(label)
            plt_temp, = plt.plot(recalls, precisions, lw=5, label=label, color=colorname)
            plot_list.append(plt_temp)
        elif mode == 'bprc':
            [aubprc, up_aubprc, brfp, bpfr] = cal_pr_values(balanced_precisions, balanced_recalls, test_precision,
                                                            test_recall)
            label = '{} ({:.3f})'.format('_'.join(name.split('_')[:-1]), aubprc)
            order_dict[aubprc] = order_index
            order_index += 1
            label_list.append(label)
            plt_temp, = plt.plot(balanced_recalls, balanced_precisions, lw=5, label=label, color=colorname)
            plot_list.append(plt_temp)

        #             plt.plot(fpr, tpr, lw=5, label='{} (PR={:.3f})'.format(name, average_precision_score([0 if i == -1 else 1 for i in df.CLASS.to_list()], df[name + '_pred'].to_list(), pos_label=1)),color = colorname)
        #             plt.plot([0, 1], [0, 1], '--', lw=5, color = 'black')
        plt.axis('square')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.xticks(fontsize=tick_font)
        plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0], fontsize=tick_font)
        plt.xlabel('Precision', fontsize=label_font)
        plt.ylabel('Recall', fontsize=label_font)
    #             if mode == 'prc':
    #                 plt.title('PRC Curve',fontsize=20)
    #             elif mode == 'bprc':
    #                 plt.title('Balanced PRC Curve',fontsize=20)
    #             plt.legend(loc='upper right',fontsize=10)
    # plt.plot([0, 1], [0, 1], '--', lw=5, color = 'black')
    # plt.legend(handles = list(reversed(plot_list)), labels = list(reversed(label_list)),
    #            loc=loc,fontsize=18)
    order_list = pd.DataFrame.from_dict({'AUC': order_dict.keys(), 'order': order_dict.values()}).sort_values(
        by='AUC', ascending=False).order.tolist()
    plt.legend(handles=[plot_list[idx] for idx in order_list],
               labels=[label_list[idx] for idx in order_list],
               loc=loc, fontsize=legend_font)

    if save is not None:
        plt.savefig(save, dpi=500, bbox_inches='tight', format='pdf')

    return plt


def get_metrics_plot(df, type_name, name_color_dict, save=None, metric_list=None, height=8, n_cols=5):
    sns.set(context=None,
            style='white',
            palette='deep',
            font='Times New Roman',
            font_scale=3,
            color_codes=True,
            rc=None)
    if len(metric_list) == 10:
        fig, ax = plt.subplots(nrows=2, ncols=n_cols, sharex=True, sharey=False, figsize=(20, height), dpi=500)
    elif len(metric_list) == n_cols:
        fig, ax = plt.subplots(nrows=1, ncols=n_cols, sharex=True, sharey=False, figsize=(20, height),
                               dpi=500)  # (20, 4.2) for manuscript
    subplot_index = 1
    for metric in metric_list:
        if len(metric_list) == 10:
            plt.subplot(2, 5, subplot_index)
        elif len(metric_list) == n_cols:
            plt.subplot(1, n_cols, subplot_index)
        plt.title(metric, fontsize=24)
        plt.xlim(None)
        sns.scatterplot(data=df, x=metric, y='model_name',
                        c=[name_color_dict[i] for i in df.model_name.tolist()], s=100)
        if subplot_index % n_cols == 1:
            plt.ylabel("Model", fontsize=label_font)
            plt.yticks(fontsize=tick_font)
        else:
            plt.ylabel(None)
            plt.yticks([])
        if len(metric_list) == 10:
            if subplot_index >= 6:
                plt.xticks([0, 0.5, 1], fontsize=tick_font)
            else:
                plt.xticks(None)
        elif len(metric_list) == n_cols:
            plt.xticks([0, 0.5, 1], fontsize=tick_font)
        subplot_index += 1
        plt.xlabel(None)
    if save:
        plt.savefig(save, dpi=500, bbox_inches='tight', format='pdf')


def visualize(test, filename):
    print('---' + time.asctime(time.localtime(time.time())) + '--- visualizing result.\n')
    for col in model_list:
        test[col] = test[col].astype('float64').tolist()
    name_color_dict = {}
    for name, color in zip(model_name_list, colors):
        name_color_dict[name] = color

    model_order = CategoricalDtype(
        model_name_list,
        ordered=True
    )

    name_threshold_dict = {}
    for name, threshold in zip(model_name_list, threshold_list):
        name_threshold_dict[name] = threshold

    #  PRC
    multi_models_pr(model_list, [name_color_dict[i] for i in model_name_list], test,
                    save=f'{save_dir}{filename}_PRC.pdf', mode='prc', loc='lower left')

    #  AUC
    df_test_roc = test.copy()
    df_test_roc.SIFT4G_score = [1 - i if str(i) != 'nan' else np.nan for i in df_test_roc.SIFT4G_score]

    train_roc_graph = multi_models_roc(model_list, [name_color_dict[i] for i in model_name_list],
                                       df_test_roc, save=f'{save_dir}{filename}_AUC.pdf')

    #  other metrics
    if cutoff:
        df_metrics_test = test[(test.AF <= cutoff)][model_list + ['CLASS']]
    else:
        df_metrics_test = test[model_list + ['CLASS']]
    df_metrics_test.columns = model_name_list + ['CLASS']
    df_metrics_test.SIFT4G = [1 - i for i in df_metrics_test.SIFT4G.tolist()]
    df_performance_test = cal_metrics_visualization(df_metrics_test, version, threshold_list, model_name_list)

    get_metrics_plot(df_performance_test, filename, name_color_dict=name_color_dict,
                     metric_list=['Missing rate', 'MCC', 'Accuracy', 'Precision', 'Recall', 'F1-score', 'F_beta-score',
                                  'G-mean', 'AUPRC', 'AUBPRC'],
                     save=f'{save_dir}{filename}_performance.pdf')
    print('---' + time.asctime(time.localtime(time.time())) + f'--- MAGPIE prediction finished. Results are in {save_dir} folder.\n')
