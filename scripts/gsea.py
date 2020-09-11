#!/nfs/production3/ma/home/manikg/miniconda3/envs/scEnv/bin python3

import sys
import os
import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import re
from gprofiler import GProfiler
import gseapy

########################################################

pathToData = sys.argv[1]
pathToDEResults = sys.argv[2]
study = sys.argv[3]
cell = sys.argv[4]
conditionComparison = sys.argv[5]
outEnrich = sys.argv[6]
#outHeat = sys.argv[7]
#outDot = sys.argv[8]
outHallmark = sys.argv[7]
outReactome = sys.argv[8]
outKegg = sys.argv[9]
pathToHallmarkGmt = sys.argv[10]
pathToReatomeGmt = sys.argv[11]
pathToKeggGmt = sys.argv[12]

#Split condition comparison to the two pairs
print(conditionComparison)
conditionComparison = conditionComparison.split("_")

########################################################
#Load the data
ad = sc.read(pathToData)

## For plotting later: Remove punctuations marks from the 'cell' and 'condition' obs
#Note that - has been ommitted as it was used to denote CD4- T cells
typeList = []
for element in ad.obs['cell']:
    typeList.append(re.sub('[^A-Za-z0-9.-]+', '_', element))
ad.obs['cellRe'] = typeList
ad.obs['cellRe'] = ad.obs['cellRe'].astype('category')

typeList = []
for element in ad.obs['condition']:
    typeList.append(re.sub('[^A-Za-z0-9.-]+', '_', element))
ad.obs['conditionRe'] = typeList  
ad.obs['conditionRe'] = ad.obs['conditionRe'].astype('category')

#Read the DE results
deAll = pd.read_csv(pathToDEResults, sep = "\t")

de = deAll[deAll.FDR<0.05].sort_values('FDR')
rankedGeneList = pd.DataFrame({'0': deAll.primerid,'1':deAll.FDR})
rankedGeneList = rankedGeneList.sort_values('1').reset_index(drop=True)

#Perform GSEA using gene sets downloaded from MSigDB https://www.gsea-msigdb.org/gsea/downloads.jsp

if len(rankedGeneList)>0:
    pre_res = gseapy.prerank(rnk=rankedGeneList, 
                             gene_sets=pathToHallmarkGmt,
                             processes=4,
                             outdir=os.path.dirname(outHallmark), format='png', seed=6)
    pre_res.res2d[pre_res.res2d.fdr<0.05].to_csv(outHallmark, sep = "\t")

    pre_res = gseapy.prerank(rnk=rankedGeneList, 
                             gene_sets=pathToReatomeGmt,
                             processes=4,
                             outdir=os.path.dirname(outReactome), format='png', seed=6)
    pre_res.res2d[pre_res.res2d.fdr<0.05].to_csv(outReactome, sep = "\t")

    pre_res = gseapy.prerank(rnk=rankedGeneList, 
                             gene_sets=pathToKeggGmt,
                             processes=4,
                             outdir=os.path.dirname(outKegg), format='png', seed=6)
    pre_res.res2d[pre_res.res2d.fdr<0.05].to_csv(outKegg, sep = "\t")

    #Perform function enrichment using Gene Ontology Biological Process
    #Taken from tutorial: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
    #Interpretation of differentially expressed genes in paneth cells - g:profiler
    if len(de)>0:
        gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

        enrichmentResults = gp.profile(organism='hsapiens', sources=['GO:BP'], user_threshold=0.05,
                                       significance_threshold_method='fdr', 
                                       background=ad.var_names.tolist(), 
                                       ordered=True,
                                       query=de['primerid'].tolist())
        enrichmentResults.to_csv(outEnrich, sep = "\t")

########################################################
#Plot and save the DE results
#ad1 = ad[(ad.obs['study'] == study) & 
#         (ad.obs['cellRe'] == cell) & 
#         (ad.obs['condition'].isin(conditionComparison))]

#sc.pl.heatmap(ad1, list(de['primerid']), groupby='condition', use_raw=False, save=outHeat.split("heatmap")[1],
#             cmap = "bwr", vmin=-3, vmax=3, show=False, show_gene_labels=True, swap_axes=True)

#sc.pl.dotplot(ad1, list(de['primerid']), groupby='condition', standard_scale='var',
#             save=outDot.split("dotplot")[1], show=False)
