#sbatch -p standard --mem=100GB --time=0-12 --account=ruic20_lab --error=${err}/${h5ad}.err --output=${err}/${h5ad}.out $comd $h5ad $indir

import scanpy as sc
import scvelo as scv
import os
import pandas as pd
import anndata as ad
import sys

h5ad=sys.argv[1]
indir=sys.argv[2]
outdir=sys.argv[3]

adata=sc.read(f'{indir}/{h5ad}.h5ad')
samplelist="/dfs3b/ruic20_lab/tingty7/projects/ocular_surface/scRNA/RNA_velocity/loom_path.csv"

sl=pd.read_csv(samplelist, index_col=0)

loom_list=[]
n=0
#loom_full=ad.AnnData()
for i in adata.obs["sampleid"].value_counts().index:
    n=n+1
    file1=sl.loc[i, sl.columns[0]]
    adata_sp=sc.read(file1)
    adata_sp.obs.index = [i + x.split(":", 1)[1][:-1] + "-1" for x in adata_sp.obs.index]
    idx=adata_sp.obs.index.intersection(adata.obs.index)
    adata_sp=adata_sp[idx].copy()

    adata_sp.var_names_make_unique()
    loom_list.append(adata_sp)

loom_full=ad.concat(loom_list)
loom_full.write(f"{outdir}/{h5ad}_sp_loom.h5ad")