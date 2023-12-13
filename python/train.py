import os.path
import joblib
import numpy as np
import pandas as pd
import time
import re
import itertools as it
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from openfe import openfe, tree_to_formula, transform
import lightgbm as lgb
from utils import get_prior, get_interpreted_x_from_y, get_interpreted_y_from_x, cal_pr_values, cal_performance_metrics
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
feature_list_new = ['phastConsElements100way', 'phyloP100way_vertebrate',
                    'phyloP20way_mammalian', 'phastCons100way_vertebrate',
                    'phastCons20way_mammalian', 'SiPhy_29way_logOdds',
                    'phyloP30way_mammalian', 'phastCons30way_mammalian', 'AF', 'AF_raw',
                    'AF_male', 'AF_female', 'AF_afr', 'AF_ami', 'AF_amr', 'AF_asj',
                    'AF_eas', 'AF_fin', 'AF_nfe', 'AF_oth', 'gdi', 'gdi_phred', 'rvis1',
                    'rvis2', 'lof_score', 'molecular_weight', 'equipotential_point',
                    'hydrophilic', 'hydrophobic', 'amphipathic ', 'cyclic', 'essential',
                    'aromatic', 'aliphatic', 'nonpolar', 'polar_uncharged', 'acidic',
                    'basic', 'sulfur', 'pka_cooh', 'pka_nh3', 'blosum100', 'DS_AG', 'DS_AL',
                    'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'Gm12878',
                    'H1hesc', 'Hepg2', 'Hmec', 'Hsmm', 'Huvec', 'K562', 'Nhek', 'Nhlf',
                    'func_frameshift', 'func_nonframeshift', 'func_nonsynonymous SNV',
                    'func_startloss', 'func_stopgain', 'func_stoploss',
                    'omim_Autosomal_dominant', 'omim_Autosomal_recessive',
                    'omim_X_linked_dominant', 'omim_X_linked_recessive', 'omim_other']
selection_params = {
    'boosting_type': 'gbdt',
    'objective': 'binary',
    'metric': {'binary_logloss', 'auc'},  # 二进制对数损失
    'num_leaves': 40,
    'max_depth': 6,
    'max_bin': 255,
    'min_data_in_leaf': 101,
    'learning_rate': 0.01,
    'feature_fraction': 1.0,
    'bagging_fraction': 1.0,
    'bagging_freq': 45,
    'lambda_l1': 0.001,
    'lambda_l2': 0.4,  # 越小l2正则程度越高
    'min_split_gain': 0.0,
    'verbose': 5,
    'is_unbalance': True
}
original_addon_list = ['Huvec', 'K562', 'Gm12878', 'H1hesc', 'Hsmm', 'Nhlf', 'Hmec', 'Hepg2', 'DS_AG', 'DS_AL', 'DS_DG',
                       'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']
threshold = 0.5
beta_score = 0.5


def load_data(df_train, y_train):
    train_x, train_y = df_train, y_train
    # 用sklearn.cross_validation进行训练数据集划分，这里训练集和交叉验证集比例为7：3，可以自己根据需要设置
    X, val_X, y, val_y = train_test_split(
        train_x,
        train_y,
        test_size=0.1,
        random_state=1,
        stratify=train_y  # 这里保证分割后y的比例分布与原数据一致
    )
    lgb_train = lgb.Dataset(X, y)
    lgb_eval = lgb.Dataset(val_X, val_y, reference=lgb_train)
    return lgb_train, lgb_eval


def train_model(df_train, y_train):
    lgb_train, lgb_eval = load_data(df_train, y_train)
    gbm = lgb.train(selection_params,
                    lgb_train,
                    num_boost_round=40,
                    valid_sets=[lgb_eval],
                    early_stopping_rounds=5)
    return gbm


def autofe(data, filename, n_blocks = 2, min_candidate_features = 2000):
    print('---' + time.asctime(time.localtime(time.time())) + '--- autoFE\n')
    train_info = data.iloc[:, :7]
    X_train = data[feature_list_new].astype('float64')
    _, X_test = train_test_split(X_train, test_size=0.1)
    Y_train = np.array(data['CLASS']).reshape(len(data['CLASS']), )

    ofe = openfe()
    # generate new features
    features = ofe.fit(data=X_train, label=Y_train, n_jobs=30,
                       # categorical_features = cell_types + ['func', 'omim'],
                       n_data_blocks=n_blocks, min_candidate_features=min_candidate_features, feature_boosting=True)

    X_train_tr, X_test_tr = transform(X_train, X_test, features, n_jobs=30)
    joblib.dump(features, os.path.join(root, f'data/result/openFE_{filename}.features'))

    X_train_tr.index = list(range(X_train_tr.shape[0]))

    with open(os.path.join(root, f'data/result/feature_name_{filename}.txt', 'w')) as f_write:
        f_write.write('name\n')
        for feature in features:
            f_write.write(tree_to_formula(feature) + '\n')
    feature_name_list = pd.read_csv(os.path.join(root, f'data/result/feature_name_{filename}.txt'), sep='\t').name.tolist()

    feature_num = X_train_tr.shape[0]
    train_round = 1

    while feature_num > 200:
        gbm = train_model(X_train_tr, Y_train)
        feature_imp = pd.DataFrame({'Value': gbm.feature_importance(), 'Feature': X_train.columns})
        if train_round >= 50:
            break
        train_round += 1
        feature_sum = feature_imp.Value.sum()
        drop_list = feature_imp[feature_imp.Value / feature_sum < 1e-3].Feature.tolist()
        X_train_tr.drop(drop_list, axis=1, inplace=True)
        feature_num = X_train_tr.shape[0]

    core_feature_list = []
    addon_feature_list = []
    for feature in X_train_tr.columns:
        if 'autoFE_f_' in feature:
            new_feature = re.split('[+*/,-]',
                                   feature_name_list[int(feature.replace('autoFE_f_', ''))].split('(')[1][:-1])
            if len(new_feature) == 1:
                if new_feature in original_addon_list:
                    addon_feature_list.append(feature)
                else:
                    core_feature_list.append(feature)
            elif len(new_feature) == 2:
                if new_feature[0] not in original_addon_list and new_feature[1] not in original_addon_list:
                    core_feature_list.append(feature)
                else:
                    addon_feature_list.append(feature)
            else:
                assert 'wrong feature format'
        else:
            if feature in original_addon_list:
                addon_feature_list.append(feature)
            else:
                core_feature_list.append(feature)

    train_X, val_X, train_y, val_y = train_test_split(
        X_train_tr,
        Y_train,
        test_size=0.1,
        random_state=214,
        stratify=Y_train  # 这里保证分割后y的比例分布与原数据一致
    )
    best_AUC_val = 0
    best_combination_val = []
    for combination_num in range(1, len(addon_feature_list)):
        for e in it.combinations(addon_feature_list, combination_num):
            X_train_it = train_X[core_feature_list + list(e)]
            X_val_it = val_X[core_feature_list + list(e)]
            gbm = train_model(X_train_it, train_y)
            val_pred = gbm.predict(X_val_it)
            AUC_val = round(roc_auc_score(list(val_y), val_pred, average='macro'), 3)

            if AUC_val > best_AUC_val:
                best_AUC_val = AUC_val
                best_combination_val = core_feature_list + list(e)

    pd.DataFrame(best_combination_val, columns=['feature']).to_csv(os.path.join(root, f'data/result/selection_{filename}.csv'),
                                                                   index=False)
    X_train_tr = X_train_tr[best_combination_val]
    return pd.concat([train_info, X_train_tr, pd.DataFrame(Y_train)], axis=1)


def step_training(data, filename):
    X_train = data.iloc[:, 7:]
    Y_train = data['CLASS']
    lgb_train = lgb.Dataset(X_train, Y_train, free_raw_data=False)

    # set init params excluding CV params
    params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': 'auc',
        'nthread': 4,
        'learning_rate': 0.1
    }
    max_auc = float('0')
    best_params = {}

    # Imporve accuracy
    for num_leaves in range(5, 100, 5):
        for max_depth in range(3, 8, 1):
            params['num_leaves'] = num_leaves
            params['max_depth'] = max_depth

            cv_results = lgb.cv(
                params,
                lgb_train,
                seed=1,
                nfold=5,
                metrics=['auc'],
                early_stopping_rounds=10,
                verbose_eval=True
            )

            mean_auc = pd.Series(cv_results['auc-mean']).max()
            boost_rounds = pd.Series(cv_results['auc-mean']).idxmax()

            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['num_leaves'] = num_leaves
                best_params['max_depth'] = max_depth
    if 'num_leaves' and 'max_depth' in best_params.keys():
        params['num_leaves'] = best_params['num_leaves']
        params['max_depth'] = best_params['max_depth']

    # Avoid over-fitting
    for max_bin in range(5, 256, 10):
        for min_data_in_leaf in range(1, 102, 10):
            params['max_bin'] = max_bin
            params['min_data_in_leaf'] = min_data_in_leaf

            cv_results = lgb.cv(
                params,
                lgb_train,
                seed=1,
                nfold=5,
                metrics=['auc'],
                early_stopping_rounds=10,
                verbose_eval=True
            )

            mean_auc = pd.Series(cv_results['auc-mean']).max()
            boost_rounds = pd.Series(cv_results['auc-mean']).idxmax()

            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['max_bin'] = max_bin
                best_params['min_data_in_leaf'] = min_data_in_leaf
    if 'max_bin' and 'min_data_in_leaf' in best_params.keys():
        params['min_data_in_leaf'] = best_params['min_data_in_leaf']
        params['max_bin'] = best_params['max_bin']

    for feature_fraction in [0.6, 0.7, 0.8, 0.9, 1.0]:
        for bagging_fraction in [0.6, 0.7, 0.8, 0.9, 1.0]:
            for bagging_freq in range(0, 50, 5):
                params['feature_fraction'] = feature_fraction
                params['bagging_fraction'] = bagging_fraction
                params['bagging_freq'] = bagging_freq

                cv_results = lgb.cv(
                    params,
                    lgb_train,
                    seed=1,
                    nfold=5,
                    metrics=['auc'],
                    early_stopping_rounds=10,
                    verbose_eval=True
                )

                mean_auc = pd.Series(cv_results['auc-mean']).max()
                boost_rounds = pd.Series(cv_results['auc-mean']).idxmax()

                if mean_auc >= max_auc:
                    max_auc = mean_auc
                    best_params['feature_fraction'] = feature_fraction
                    best_params['bagging_fraction'] = bagging_fraction
                    best_params['bagging_freq'] = bagging_freq

    if 'feature_fraction' and 'bagging_fraction' and 'bagging_freq' in best_params.keys():
        params['feature_fraction'] = best_params['feature_fraction']
        params['bagging_fraction'] = best_params['bagging_fraction']
        params['bagging_freq'] = best_params['bagging_freq']

    for lambda_l1 in [1e-5, 1e-3, 1e-1, 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
        for lambda_l2 in [1e-5, 1e-3, 1e-1, 0.0, 0.1, 0.4, 0.6, 0.7, 0.9, 1.0]:
            params['lambda_l1'] = lambda_l1
            params['lambda_l2'] = lambda_l2
            cv_results = lgb.cv(
                params,
                lgb_train,
                seed=1,
                nfold=5,
                metrics=['auc'],
                early_stopping_rounds=10,
                verbose_eval=True
            )

            mean_auc = pd.Series(cv_results['auc-mean']).max()
            boost_rounds = pd.Series(cv_results['auc-mean']).idxmax()

            if mean_auc >= max_auc:
                max_auc = mean_auc
                best_params['lambda_l1'] = lambda_l1
                best_params['lambda_l2'] = lambda_l2
    if 'lambda_l1' and 'lambda_l2' in best_params.keys():
        params['lambda_l1'] = best_params['lambda_l1']
        params['lambda_l2'] = best_params['lambda_l2']

    for min_split_gain in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        params['min_split_gain'] = min_split_gain

        cv_results = lgb.cv(
            params,
            lgb_train,
            seed=1,
            nfold=5,
            metrics=['auc'],
            early_stopping_rounds=10,
            verbose_eval=True
        )

        mean_auc = pd.Series(cv_results['auc-mean']).max()
        boost_rounds = pd.Series(cv_results['auc-mean']).idxmax()

        if mean_auc >= max_auc:
            max_auc = mean_auc

            best_params['min_split_gain'] = min_split_gain
    if 'min_split_gain' in best_params.keys():
        params['min_split_gain'] = best_params['min_split_gain']

    lgb_train, lgb_eval = load_data(X_train, Y_train)
    params = {
        'boosting_type': 'gbdt',
        'objective': 'binary',
        'metric': {'binary_logloss', 'auc'},  # 二进制对数损失
        'num_leaves': 40,
        'max_depth': 6,
        'max_bin': 255,
        'min_data_in_leaf': 101,
        'learning_rate': 0.01,
        'feature_fraction': 1.0,
        'bagging_fraction': 1.0,
        'bagging_freq': 45,
        'lambda_l1': 0.001,
        'lambda_l2': 0.4,  # 越小l2正则程度越高
        'min_split_gain': 0.0,
        'verbose': 5,
        'is_unbalance': True
    }
    for key in best_params.keys():
        if key == 'max_depth':
            params[key] = best_params[key]
        elif key == 'max_leaves':
            params[key] = best_params[key]
        else:
            params[key] = best_params[key]

    gbm = lgb.train(params,
                    lgb_train,
                    num_boost_round=1000,
                    valid_sets=lgb_eval,
                    early_stopping_rounds=50)
    joblib.dump(gbm, os.path.join(root, f'data/result/MAGPIE_{filename}.model'))


def train(train_file):
    data = pd.read_csv(train_file)
    filename = os.path.splitext(os.path.basename(train_file))[0]
    data_bpca = pd.read_csv(os.path.join(os.path.abspath(train_file), f'{os.path.splitext(os.path.basename(train_file))[0]}_bpca.csv'), header = None)
    data_bpca.columns = data.columns[7:]
    data = pd.concat([data.iloc[:, :7], data_bpca], axis = 1)
    data = autofe(data, filename)
    step_training(data, filename)


