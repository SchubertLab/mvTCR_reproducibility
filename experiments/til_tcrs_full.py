"""
python -u haniffa_tcr_optuna.py --model moe
"""
# comet-ml must be imported before torch and sklearn
import comet_ml

import sys
sys.path.append('../mvTCR/')

from tcr_embedding.models.model_selection import run_model_selection
import tcr_embedding.utils_training as utils
from tcr_embedding.utils_preprocessing import group_shuffle_split, encode_tcr

import os
import argparse


utils.fix_seeds(42)

parser = argparse.ArgumentParser()
parser.add_argument('--rna_weight', type=int, default=1)
parser.add_argument('--model', type=str, default='moe')
parser.add_argument('--gpus', type=int, default=1)
parser.add_argument('--conditional', type=str, default='patient')

args = parser.parse_args()


adata = utils.load_data('TIL/processed.h5ad')

len_beta = adata.obs['IR_VDJ_1_junction_aa'].str.len().max()
len_alpha= adata.obs['IR_VJ_1_junction_aa'].str.len().max()
pad = max(len_beta, len_alpha)

encode_tcr(adata, 'IR_VJ_1_junction_aa', 'IR_VDJ_1_junction_aa', pad)

if args.conditional != 'None':
    from sklearn.preprocessing import OneHotEncoder
    enc = OneHotEncoder(sparse=False)
    enc.fit(adata.obs[args.conditional].to_numpy().reshape(-1, 1))
    adata.obsm[args.conditional] = enc.transform(adata.obs[args.conditional].to_numpy().reshape(-1, 1))

random_seed = 42
train, val = group_shuffle_split(adata, group_col='clonotype', val_split=0.20, random_seed=random_seed)
adata.obs['set'] = 'train'
adata.obs.loc[val.obs.index, 'set'] = 'val'
adata = adata[adata.obs['set'].isin(['train', 'val'])]


params_experiment = {
    'study_name': f'til_full_{args.model}_{args.rna_weight}_{args.conditional}',
    'comet_workspace': None,
    'model_name': args.model,
    'balanced_sampling': 'clonotype',
    'metadata': ['clonotype', 'meta.cluster'],
    'save_path': os.path.join(os.path.dirname(__file__), '..', 'optuna',
                              f'til_full_{args.model}_{args.rna_weight}_{args.conditional}'),
    'conditional': args.conditional,
}
if args.model == 'rna':
    params_experiment['balanced_sampling'] = None

params_optimization = {
    'name': 'pseudo_metric',
    'prediction_labels':
        {'clonotype': 1,
         'meta.cluster': args.rna_weight}
}

timeout = (2 * 24 * 60 * 60) - 300
run_model_selection(adata, params_experiment, params_optimization, None, timeout, args.gpus)
