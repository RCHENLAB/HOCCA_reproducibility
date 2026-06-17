import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import sys
#scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


dir1=sys.argv[1]
bname=sys.argv[2]

import os
if not os.path.exists(dir1):
    os.makedirs(dir1)


adata=sc.read(f'{dir1}/{bname}_velocity.h5ad')

scv.pl.velocity_graph(adata, threshold=.1, save=f'{dir1}/{bname}_velocity_graph.png', color="celltype")

df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
df.to_csv(f'{dir1}/{bname}_velocity_gene.csv')

scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], save=f'{dir1}/{bname}_cell_cycle.png')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],save=f'{dir1}/{bname}_velocity_length_confidence.png')

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save=f'{dir1}/{bname}_cell_transitions.png')

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,save=f'{dir1}/{bname}_paga_transition.png')

