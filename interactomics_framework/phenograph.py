import numpy as np
import scanpy as sc
import pandas as pd
import phenograph
import random
import csv

random.seed(42)

# from guided pbmc3k tutorial of PhenoGraph 
# https://github.com/dpeerlab/PhenoGraph/blob/master/examples/tutorial_pbmc3k.ipynb


# load data for clustering 


k = 5 # same as k for the SNN graph in bluster from Domi

communities, graph, Q = phenograph.cluster(pd.DataFrame(adata.obsm['X_pca']),k=k) # run PhenoGraph

# store the results in adata:
adata.obs['PhenoGraph_clusters'] = pd.Categorical(communities)
adata.uns['PhenoGraph_Q'] = Q
adata.uns['PhenoGraph_k'] = k
