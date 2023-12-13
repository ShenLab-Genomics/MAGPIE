import pandas as pd
import numpy as np
import time
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn import metrics

feature_list = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'CLASS', 'Gene.refGene', 'func', 'omim',
                'phastConsElements100way', 'phyloP100way_vertebrate', 'phyloP20way_mammalian',
                'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'SiPhy_29way_logOdds',
                'phyloP30way_mammalian', 'phastCons30way_mammalian', 'AF', 'AF_raw', 'AF_male',
                'AF_female', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj', 'AF_eas', 'AF_fin', 'AF_nfe',
                'AF_oth', 'gdi', 'gdi_phred', 'rvis1', 'rvis2', 'lof_score', 'molecular_weight',
                'equipotential_point', 'hydrophilic', 'hydrophobic', 'amphipathic ', 'cyclic',
                'essential', 'aromatic', 'aliphatic', 'nonpolar', 'polar_uncharged', 'acidic',
                'basic', 'sulfur', 'pka_cooh', 'pka_nh3', 'blosum100', 'DS_AG', 'DS_AL', 'DS_DG',
                'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'Gm12878', 'H1hesc', 'Hepg2',
                'Hmec', 'Hsmm', 'Huvec', 'K562', 'Nhek', 'Nhlf']


def binarySearch(list1, target):
    left = 0
    right = len(list1) - 1

    while left <= right:
        mid = (right + left) // 2
        if mid + 1 < len(list1) and list1[mid] <= target < list1[mid + 1]:
            return mid
        elif target > list1[mid]:
            left = mid + 1
        else:  # target < list[mid]
            right = mid - 1
    return -1


def MissingUniqueStatistics(df):
    total_entry_list = []
    total_missing_value_list = []
    missing_value_ratio_list = []
    data_type_list = []
    unique_values_list = []
    number_of_unique_values_list = []
    variable_name_list = []
    for col in df.columns:
        variable_name_list.append(col)
        missing_value_ratio = round((df[col].isna().sum() / len(df[col])), 4)
        total_entry_list.append(df[col].shape[0] - df[col].isna().sum())
        total_missing_value_list.append(df[col].isna().sum())
        missing_value_ratio_list.append(missing_value_ratio)
        data_type_list.append(df[col].dtype)
        unique_values_list.append(list(df[col].unique()))
        number_of_unique_values_list.append(len(df[col].unique()))
    data_info_df = pd.DataFrame({'Variable': variable_name_list,
                                 '#_Total_Entry': total_entry_list,
                                 '#_Missing_Value': total_missing_value_list,
                                 '%_Missing_Value': missing_value_ratio_list,
                                 'Data_Type': data_type_list, 'Unique_Values': unique_values_list,
                                 '#_Uniques_Values': number_of_unique_values_list})
    return data_info_df.sort_values(by='#_Missing_Value', ascending=False)


def preprocess(data):
    data_info = MissingUniqueStatistics(data).set_index('Variable')
    categorical_columns = data_info[data_info['Data_Type'] == 'object'].index.tolist()
    for col in categorical_columns:
        list1 = []
        for i in data[col].tolist():
            if i == '.' or i == '-':
                list1.append(np.nan)
            else:
                list1.append(i)
        data[col] = list1
    return data


def get_prior(y):
    y = np.array(y)
    return (y == 1).sum() / len(y)


def get_interpreted_x_from_y(new_y, x, y, type='first_intersection'):
    new_x = np.nan

    if type == 'last_intersection':
        x = x[::-1]
        y = y[::-1]

    for i in range(len(x)):
        if y[i] == new_y:
            new_x = x[i]
            break
        if i < len(x) - 1:
            if (new_y - y[i]) * (new_y - y[i + 1]) < 0:
                new_x = x[i] + (new_y - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
                break
    return new_x


def get_interpreted_y_from_x(new_x, x, y, type='first_intersection'):
    new_y = np.nan
    if type == 'last_intersection':
        x = x[::-1]
        y = y[::-1]

    for i in range(len(y)):
        if x[i] == new_x:
            new_y = y[i]
            break
        if i < len(y) - 1:
            if (new_x - x[i]) * (new_x - x[i + 1]) < 0:
                new_y = y[i] + (new_x - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                break
    return new_y


def cal_pr_values(precisions, recalls, test_precision=0.9, test_recall=0.9):
    precisions = np.array(precisions)
    recalls = np.array(recalls)
    auprc = -np.sum(np.diff(recalls) * np.array(precisions)[:-1])  # average precision
    up_precisions = precisions - test_precision
    new_up_precisions = up_precisions[up_precisions > 0]
    new_up_recalls = recalls[up_precisions > 0]
    up_auprc = -np.sum(np.diff(new_up_recalls) * np.array(new_up_precisions)[:-1])

    rfp = get_interpreted_x_from_y(test_precision, recalls, precisions, type='first_intersection')
    pfr = get_interpreted_y_from_x(test_recall, recalls, precisions, type='last_intersection')

    return [auprc, up_auprc, rfp, pfr]


def cal_performance_metrics(df_metrics, df_performance, model, label):
    test_precision, test_recall = [0.9, 0.9]

    df = df_metrics.loc[:, [model, label]].astype('float64')
    df_auc = df_metrics.loc[:, [model, label]].astype('float64')
    df[model] = [1 if i >= threshold else -1 for i in df[model].tolist()]
    df[model] = df[model].astype('int64')

#     print(df[label].value_counts())
#     print(df[model].value_counts())

    C = confusion_matrix(df[label].tolist(), df[model].tolist())
    prior = get_prior(df[label].tolist())
    TN, FP, FN, TP = C.ravel()
    accuracy = round((TP + TN) / (TP + FP + TN + FN), 3)
    precision = round(TP / (TP + FP), 3)
    recall = round(TP / (TP + FN), 3)
    specificity = round(TN / (TN + FP), 3)
    f1_score = round(2 * precision * recall / (precision + recall), 3)
    fbeta_score = round(((1 + beta_score ** 2) * precision * recall / (beta_score ** 2 * precision + recall)), 3)
    g_mean = round((recall * specificity) ** 0.5, 3)
    mcc = round((TP * TN - FP * FN) / (((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5), 3)
    AUC = round(roc_auc_score(df_auc[label].tolist(), df_auc[model].tolist(), average='macro'), 3)
    precisions, recalls, prc_thresholds = metrics.precision_recall_curve(df_auc[label].tolist(),
                                                                         df_auc[model].tolist())
    recalls = np.insert(recalls, 0, 1)
    precisions = np.insert(precisions, 0, prior)
    balanced_precisions = precisions * (1 - prior) / (precisions * (1 - prior) + (1 - precisions) * prior)
    balanced_recalls = recalls
    [AUPRC, UP_AURPC, rfp, pfr] = cal_pr_values(precisions, recalls, test_precision, test_recall)
    [AUBPRC, UP_AUBPRC, brfp, bpfr] = cal_pr_values(balanced_precisions, balanced_recalls, test_precision,
                                                    test_recall)
#     print([TP, FP, TN, FN, accuracy, precision, recall, specificity, f1_score, fbeta_score, g_mean, mcc, AUC, AUPRC,
#            AUBPRC, UP_AUBPRC, UP_AURPC])
    df_performance.loc[model] = {'AUPRC': AUPRC, 'AUBPRC': AUBPRC, 'UP_AUBPRC': UP_AUBPRC,
                                 'AUC': AUC, 'TP': TP, 'FP': FP, 'TN': TN, 'FN': FN,
                                 'Accuracy': accuracy, 'Precision': precision, 'Recall': recall,
                                 'Specificity': specificity, 'F1-score': f1_score,
                                 'F_beta-score': fbeta_score, 'G-mean': g_mean, 'MCC': mcc,
                                 'Missing rate': 1 - df_metrics[model].count() / df_metrics.shape[
                                     0] * 1.0, 'Sample': df_metrics[model].count()}
    return df_performance


def cal_metrics_visualization(df_metrics, version, inner_threshold_list, model_name_list, sort=None):
    beta_score = 0.5
    test_precision, test_recall = [0.9, 0.9]
    df_performance = pd.DataFrame(columns=['AUC', 'MCC', 'TP', 'FP', 'TN', 'FN', 'Accuracy',
                                           'Precision', 'Recall', 'Specificity', 'F1-score', 'F_beta-score',
                                           'G-mean', 'AUPRC', 'UP_AUBPRC', 'AUBPRC', 'Missing rate', 'Sample'])
    # feature_name = next(feature_iter)
    # data = data1[data1[feature_name] != '.']
    # feature_list.append(feature_name)
    missing_dict = {}
    for model, threshold in zip(model_name_list, inner_threshold_list):
        if version == 'old':
            #  The orginal version of performance, drop samples containing missing values
            df = df_metrics.loc[:, [model, 'CLASS']].dropna().astype('float64')
            df_auc = df.copy()
        elif version == 'new':
            #  New version of performance, use benign label to fill nan values.
            df = df_metrics.loc[:, [model, 'CLASS']].fillna(-1).astype('float64')
            df_auc = df.copy()
        else:
            assert 'Wrong parameter of version.'

        df[model] = [1 if i >= threshold else -1 for i in df[model].tolist()]
        df[model] = df[model].astype('int64')
        C = confusion_matrix(df['CLASS'].tolist(), df[model].tolist(), labels=[1, -1])
#         print(C)
#         print(C.ravel())
#         print(df[(df.CLASS == 1) & (df[model] == 1)].shape[0])

        prior = get_prior(df['CLASS'].tolist())
        TP, FN, FP, TN = C.ravel()
        accuracy = round((TP + TN) / (TP + FP + TN + FN), 3)
        precision = round(TP / (TP + FP), 3)
        recall = round(TP / (TP + FN), 3)
        specificity = round(TN / (TN + FP), 3)
        f1_score = round(2 * precision * recall / (precision + recall), 3)
        fbeta_score = round(((1 + beta_score ** 2) * precision * recall / (beta_score ** 2 * precision + recall)), 3)
        g_mean = round((recall * specificity) ** 0.5, 3)
        mcc = round((TP * TN - FP * FN) / (((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5), 3)
        if np.isnan(mcc):
            mcc = 0
        if np.isnan(g_mean):
            g_mean = 0
        try:
            AUC = round(roc_auc_score(df_auc['CLASS'].tolist(), df_auc[model].tolist(), average='macro'), 3)
            precisions, recalls, prc_thresholds = metrics.precision_recall_curve(df_auc['CLASS'].tolist(),
                                                                                 df_auc[model].tolist())
            recalls = np.insert(recalls, 0, 1)
            precisions = np.insert(precisions, 0, prior)
            balanced_precisions = precisions * (1 - prior) / (precisions * (1 - prior) + (1 - precisions) * prior)
            balanced_recalls = recalls
            [AUPRC, UP_AURPC, rfp, pfr] = cal_pr_values(precisions, recalls, test_precision, test_recall)
            [AUBPRC, UP_AUBPRC, brfp, bpfr] = cal_pr_values(balanced_precisions, balanced_recalls, test_precision,
                                                            test_recall)
        except:
            AUC, AUPRC, AUBPRC, UP_AUBPRC, UP_AURPC = [0] * 5
#         print([TP, FP, TN, FN, accuracy, precision, recall, specificity, f1_score, fbeta_score, g_mean, mcc, AUC, AUPRC,
#                AUBPRC, UP_AUBPRC, UP_AURPC])
        df_performance.loc[model.split('_p')[0]] = {'AUPRC': AUPRC, 'AUBPRC': AUBPRC, 'UP_AUBPRC': UP_AUBPRC,
                                                    'AUC': AUC, 'TP': TP, 'FP': FP, 'TN': TN, 'FN': FN,
                                                    'Accuracy': accuracy, 'Precision': precision, 'Recall': recall,
                                                    'Specificity': specificity, 'F1-score': f1_score,
                                                    'F_beta-score': fbeta_score, 'G-mean': g_mean, 'MCC': mcc,
                                                    'Missing rate': 1 - df_metrics[model].count() / df_metrics.shape[
                                                        0] * 1.0, 'Sample': df_metrics[model].count()}
    if sort:
        df_performance = df_performance.sort_values(by=sort, ascending=False)
    else:
        df_performance = df_performance.sort_values(by=['AUC', 'Accuracy'], ascending=False)
    df_performance.insert(0, 'model_name', df_performance.index.tolist())
    df_performance.insert(df_performance.shape[-1], 'variant_index',
                          list(reversed(range(1, df_performance.shape[0] + 1))))
    df_performance.insert(df_performance.shape[-1], 'Missing sample',
                          [df_metrics.shape[0] - i for i in df_performance.Sample.tolist()])
    # df_performance['model_name'] = df_performance['model_name'].astype(model_order)
    return df_performance
