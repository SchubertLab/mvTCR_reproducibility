# Reproducibility
This repo contains the code to reproduce the results from XXX

# Setup
To reproduce the results, clone this repo. Additionally, you will need to clone the mvTCR code from the release XYZ.

```
clone command todo
```

## TODO: MAC => ask yang
## Linux
Please run the following commands for a Linux-based OS:
```
conda create --name mvTCR python=3.8.8 -y 
conda activate mvTCR 
pip install -r requirements.txt 
conda install nb_conda_kernels -y 
conda install bottleneck -y
```

## Windows
Please comment torch from the requirements.txt, i.e. write a # before torch. Then execute the line to install all the requirements except PyTorch:

``` 
conda create --name mvTCR python=3.8.8 -y
conda activate mvTCR
pip install -r requirements.txt 
conda install nb_conda_kernels -y
```

Then install PyTorch 1.8.0 with the correct CUDA Version following the command here: https://pytorch.org/get-started/previous-versions/

# Preprocessing
TODO: copy over once done

The folder preprocessing contains jupyter notebooks for preprocessing the single cell data. The links to obtain the data are referenced in the corresponding publiction and in our paper.

# Experiments
TODO: copy over once done 

The folder experiments contains the settings to train the mvTCR model on all datasets mentioned in the paper. Typically, the hyperparameter optimization was conducted for 48 GPU-hours. After HPO, you will need to copy the best resulting model indicated by the console output, or SQL database, to the corresponding folder, as indicated in the evaluation files. Based on different runtimes on the compute clusters, different hardware components, and CUDA-versions, the resulting models might diver from the models used in the publication. We therefore also provide all trained models under TODO: zenodo link.

# Evaluation
TODO: copy over once done

The notebooks in the folder evaluation are used to derive the shared embedding and therefore requiere the trained models and the processed data. Since these notebooks save intermediate results for creating the figures, you will need to run these notebooks before you can recreate the figures.

# Figures
TODO: copy over once done

All figures (except concept figure 1a) of the manuscript can be recreated by running the notebooks in the figure-Folder. This requieres to run the evaluation first on the preprocessed data and the models downloaded from Zenodo. Note, that due to reproducability issues in the underlying library, the UMAP visualisation might not be fully recreated. However, the should qualitatively match the results shown in the paper.
