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

import os
import argparse


utils.fix_seeds(42)

parser = argparse.ArgumentParser()
parser.add_argument('--model', type=str, default='poe')
parser.add_argument('--split', type=int, default=0)
parser.add_argument('--gpus', type=int, default=1)
parser.add_argument('--n_samples', type=int, default=None)
args = parser.parse_args()

adata = utils.load_data('minervina/01_annotated_data.h5ad')

# subsample to get statistics
random_seed = args.split
adata.obs['group_col'] = [seq[1:-1] for seq in adata.obs['IR_VDJ_1_junction_aa']]
train_val, test = group_shuffle_split(adata, group_col='group_col', val_split=0.20, random_seed=random_seed)
train, val = group_shuffle_split(train_val, group_col='group_col', val_split=0.25, random_seed=random_seed)

adata.obs['set'] = None
adata.obs.loc[train.obs.index, 'set'] = 'train'
adata.obs.loc[val.obs.index, 'set'] = 'val'
adata.obs.loc[test.obs.index, 'set'] = 'test'
adata = adata[adata.obs['set'].isin(['train', 'val'])]


params_experiment = {
    'study_name': f'minervina_{args.model}_split_{args.split}',
    'comet_workspace': None,
    'model_name': args.model,
    'balanced_sampling': 'clonotype',
    'metadata': [],
    'save_path': os.path.join(os.path.dirname(__file__), '..', 'optuna',
                              f'minervina_{args.model}_split_{args.split}'),
    'beta_only': True,
    'conditional': 'donor'
}

if args.model == 'rna':
    params_experiment['balanced_sampling'] = None

params_optimization = {
    'name': 'knn_prediction',
    'prediction_column': 'epitope',
}

timeout = (2 * 24 * 60 * 60) - 300
run_model_selection(adata, params_experiment, params_optimization, None, timeout, args.gpus)
