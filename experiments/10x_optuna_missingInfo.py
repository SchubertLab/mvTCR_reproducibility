"""
python -u 10x_optuna.py --model poe --donor 1 --split 0
"""
# comet-ml must be imported before torch and sklearn
import comet_ml

import sys
sys.path.append('../mvTCR/')

from tcr_embedding.models.model_selection import run_model_selection
import tcr_embedding.utils_training as utils
from tcr_embedding.utils_preprocessing import group_shuffle_split

import scanpy as sc

import os
import argparse
import config.constants_10x as const
from sklearn.preprocessing import OneHotEncoder


utils.fix_seeds(42)

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='poe')
parser.add_argument('--donor', type=str, default=None)
parser.add_argument('--filter_non_binder', type=bool, default=True)
parser.add_argument('--split', type=int, default=0)
parser.add_argument('--gpus', type=int, default=1)
args = parser.parse_args()


adata = utils.load_data('10x')
if str(args.donor) != 'None':
    adata = adata[adata.obs['donor'] == f'donor_{args.donor}']
else:
    enc = OneHotEncoder(sparse=False)
    enc.fit(adata.obs['donor'].to_numpy().reshape(-1, 1))
    adata.obsm['donor'] = enc.transform(adata.obs['donor'].to_numpy().reshape(-1, 1))

if args.filter_non_binder:
    adata = adata[adata.obs['binding_name'].isin(const.HIGH_COUNT_ANTIGENS)]

# subsample to get statistics
random_seed = args.split
train_val, test = group_shuffle_split(adata, group_col='clonotype', val_split=0.20, random_seed=random_seed)
train, val = group_shuffle_split(train_val, group_col='clonotype', val_split=0.25, random_seed=random_seed)

pad_len = adata.obsm['beta_seq'].shape[1]
for chain in ['alpha', 'beta']:
    mask = np.random.choice(len(adata), replace=False, size=int(len(adata) * 0.15))
    adata.obsm[f'{chain}_seq'][mask] = np.zeros(pad_len)
    adata.obs[f'{chain}_len'][mask] = 0

mask = np.random.choice(len(adata), replace=False, size=int(len(adata) * 0.15))
adata.X[mask] = np.zeros(adata.X.shape[1])

adata.obs['set'] = None
adata.obs.loc[train.obs.index, 'set'] = 'train'
adata.obs.loc[val.obs.index, 'set'] = 'val'
adata.obs.loc[test.obs.index, 'set'] = 'test'
adata = adata[adata.obs['set'].isin(['train', 'val'])].copy()

n_samples_prefix = '' if args.n_samples is None else f'_{args.n_samples}'

params_experiment = {
    'study_name': f'10x_{args.donor}_{args.model}_filtered_{args.filter_non_binder}_split_{args.split}{n_samples_prefix}',
    'comet_workspace': None,
    'model_name': args.model,
    'balanced_sampling': 'clonotype',
    'metadata': ['binding_name', 'clonotype', 'donor'],
    'save_path': os.path.join(os.path.dirname(__file__), '..', 'optuna',
                              f'10x_missing_{args.donor}_{args.model}_split_{args.split}'),
    'conditional': 'donor' if str(args.donor) == 'None' else None
}

if args.model == 'rna':
    params_experiment['balanced_sampling'] = None

params_optimization = {
    'name': 'knn_prediction',
    'prediction_column': 'binding_name',
}

timeout = (2 * 24 * 60 * 60) - 300
run_model_selection(adata, params_experiment, params_optimization, None, timeout, args.gpus)
