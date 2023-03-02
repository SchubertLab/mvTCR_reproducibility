"""
python -u haniffa_tcr_optuna.py --model moe
"""
# comet-ml must be imported before torch and sklearn
import comet_ml

import scanpy as sc
import sys
sys.path.append('../mvTCR/')

from tcr_embedding.models.model_selection import run_model_selection
import tcr_embedding.utils_training as utils
from tcr_embedding.utils_preprocessing import group_shuffle_split

import os
import argparse


utils.fix_seeds(42)

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='moe')
parser.add_argument('--gpus', type=int, default=1)
parser.add_argument('--seed', type=int, default=0)

args = parser.parse_args()


adata = utils.load_data('haniffa')

# subsample to get statistics
random_seed = args.seed
sc.pp.subsample(adata, n_obs=20000, random_state=random_seed)
train, val = group_shuffle_split(adata, group_col='cdr3_beta', val_split=0.25, random_seed=random_seed)

adata.obs['set'] = 'train'
adata.obs.loc[val.obs.index, 'set'] = 'val'
adata = adata[adata.obs['set'].isin(['train', 'val'])]

params_experiment = {
    'study_name': f'haniffa_beta_{args.model}_split_{args.seed}',
    'comet_workspace': None,  # 'Covid',
    'model_name': args.model,
    'balanced_sampling': 'cdr3_beta',
    'metadata': ['clonotype', 'full_clustering'],
    'save_path': os.path.join(os.path.dirname(__file__), '..', 'optuna',
                              f'haniffa_beta_{args.model}_split_{args.seed}'),
    'conditional': 'patient_id',
    'beta_only': True
}
if args.model == 'rna':
    params_experiment['balanced_sampling'] = None

params_optimization = {
    'name': 'pseudo_metric',
    'prediction_labels':
        {'clonotype': 1,
         'full_clustering': 1}
}

timeout = (2 * 24 * 60 * 60) - 300
run_model_selection(adata, params_experiment, params_optimization, None, timeout, args.gpus)
