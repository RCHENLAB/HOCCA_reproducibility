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

sp_loom_h5ad=sys.argv[1]
h5ad=sys.argv[2]
bname=sys.argv[3]
hvg=int(sys.argv[4])
dir1=f"/dfs3b/ruic20_lab/tingty7/projects/ocular_surface/scRNA_wo_fetal/RNA_velocity/3scvelo_general/{bname}"

import os
if not os.path.exists(dir1):
    os.makedirs(dir1)

adata=sc.read(f"{h5ad}")

ldata = sc.read(f"{sp_loom_h5ad}")

adata = scv.utils.merge(adata, ldata)

scv.pl.proportions(adata, save=f'{dir1}/{bname}_splice_unsplice.png')

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=hvg)

scv.pp.moments(adata, n_pcs=30, n_neighbors=30, use_rep='X_scVI')

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

adata.write(f'{dir1}/{bname}_velocity_pre.h5ad')

scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'{dir1}/{bname}_velocity.png', color="celltype")

df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

plt.tight_layout()
plt.style.use('default')  # Reset to default style
plt.gca().set_facecolor('white')  # Set subplot background to white
plt.savefig(f'{dir1}/{bname}_Kinetic_rate_paramters_histograms.png', dpi=300, bbox_inches='tight', transparent=True)

df1=scv.get_df(adata, 'fit*', dropna=True).head()
df1.to_csv(f'{dir1}/{bname}_Kinetic_rate_paramters.csv')


scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save=f'{dir1}/{bname}_pseudo_time.png')

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save=f'{dir1}/{bname}_Latent_time.png')

adata.write(f'{dir1}/{bname}_velocity.h5ad')

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100, save=f'{dir1}/{bname}_top_gene_heatmap.png')


top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False,save=f'{dir1}/{bname}_top_gene_curve.png',color="celltype")

scv.tl.rank_dynamical_genes(adata, groupby='celltype')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.to_csv(f'{dir1}/{bname}_rank_dynamical_genes.png')
