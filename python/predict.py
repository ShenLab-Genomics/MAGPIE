import joblib
import pandas as pd
import numpy as np
import os
import time
from openfe import openfe, transform

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def predict(test, autoFE_features, selection, model_file, filename):
    feature_list = ['phastConsElements100way', 'phyloP100way_vertebrate',
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

    for col in feature_list:
        test[col] = test[col].replace('-', np.nan).replace('.', np.nan)
    features = joblib.load(autoFE_features)
    X_test = test[feature_list].astype('float64')

    print('---' + time.asctime(time.localtime(time.time())) + '--- transforming dataset\n')
    _, X_test_tr = transform(X_test, X_test, features, n_jobs=30)
    feature_list = pd.read_csv(selection).feature.tolist()
    X_test_filtered = X_test_tr[feature_list]

    print('---' + time.asctime(time.localtime(time.time())) + '--- predicting\n')
    model = joblib.load(model_file)
    test_pred = model.predict(X_test_filtered)
    pd.concat([test, pd.DataFrame(test_pred, columns=['MAGPIE_pred'])], axis=1).to_csv(os.path.join(root, f'data/result/{filename}.csv'), index=False)
    return pd.concat([test, pd.DataFrame(test_pred, columns=['MAGPIE_pred'])], axis=1)
