#!/usr/bin/env python
# coding: utf-8

# In[1]:

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc  # v1.4.3
import pyvdj as pv
import random
from scipy import stats
import upsetplot
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import seaborn as sns
import helper
from statannot import add_stat_annotation
import itertools
import os

study = str(sys.argv[1])
pathToData = str(sys.argv[2])
receptor = str(sys.argv[3])
currentConditionColumn=str(sys.argv[4])
outFile = str(sys.argv[5])

# In[2]:

if ((study=='Liao')and(receptor=='BCR')):
    with open(outFile, 'w') as fp: 
        print('No data!')
else:
    
    sc.settings.set_figure_params(dpi=300)


    # In[3]:
    fileName = pathToData+study+'.h5'
    adata_all = sc.read(fileName)
    #adata_all.shape


    # In[4]:

    pyvdj = receptor+'vdj'
    has_vdjdata_col = receptor+'_has_vdjdata'

    # In[7]:


    adata_all.obs['sample'].isin(adata_all.uns[pyvdj]['df']['sample'].unique()).value_counts()


    # In[8]:


    #For TCR analysis
    if receptor=='TCR':
        adata = adata_all[#(adata_all.obs['cell'] != 'doublets') & 
                          (adata_all.obs.cell.str.contains(' T')) & 
                          (adata_all.obs['sample'].isin(adata_all.uns[pyvdj]['df']['sample'].unique()))]
    elif receptor=='BCR':
        adata = adata_all[#(adata_all.obs['cell'] != 'doublets') & 
                          (adata_all.obs.cell.str.contains('^B', regex=True) | 
                          adata_all.obs.cell.str.contains('Plasma')) & 
                          (adata_all.obs['sample'].isin(adata_all.uns[pyvdj]['df']['sample'].unique()))]
    else:
        print("Undefined receptor")

    adata.obs['cell'] = adata.obs['cell'].cat.remove_unused_categories()
    adata.obs[currentConditionColumn] = adata.obs[currentConditionColumn].cat.remove_unused_categories()
    adata.obs['sample'] = adata.obs['sample'].cat.remove_unused_categories()

    # In[9]:


    adata.obs.cell.unique()


    # In[10]:


    #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)
    ##adata = adata[:, adata.var.highly_variable]
    ##sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
    ##sc.pp.scale(adata, max_value=10)
    #sc.tl.pca(adata, svd_solver='arpack')
    #sc.pl.pca_variance_ratio(adata, log=True)
    #sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    #sc.tl.umap(adata)


    # In[11]:


    has_vdjdata_col


    # In[12]:


    sc.settings.set_figure_params(dpi=300)
    sc.pl.umap(adata, color=['cell'], title=['Cell-type: '+study+' et al.'], save='_'+study+'_'+receptor+'_cellType.png')


    # In[13]:


    adata = helper.add_obs(adata, obs=['clonotype', 'is_clone', 'clone_count'], pyvdj = pyvdj, has_vdjdata_col = has_vdjdata_col)


    # In[14]:


    #adata.obs[has_vdjdata_col] = adata.obs[has_vdjdata_col].astype('category')
    #adata.obs['vdj_is_clone'] = adata.obs['vdj_is_clone'].astype('category')


    # In[15]:


    sc.settings.set_figure_params(dpi=300)
    sc.pl.umap(adata, color=[has_vdjdata_col, 'vdj_is_clone'], 
               title=[receptor+' has V(D)J data', 'Is clone'], 
               palette=['grey', 'blue'],
               save='_'+study+'_'+receptor+'_has_vdj_'+currentConditionColumn+'.png')


    # In[16]:


    allcell = len(adata.obs[has_vdjdata_col])
    vdjcell = sum(adata.obs[has_vdjdata_col])
    print('We have %d cells with VDJ data, out of %d cells.' % (vdjcell, allcell))


    # In[17]:


    allcell = len(adata.obs['vdj_is_clone'])
    vdjcell = sum(adata.obs['vdj_is_clone'])
    print('We have %d cells with clone, out of %d cells.' % (vdjcell, allcell))


    # In[18]:


    #adata.obs[has_vdjdata_col] = adata.obs[has_vdjdata_col].astype(int)
    adata = helper.add_obs(adata, obs=['all_productive', 'any_productive'], pyvdj = pyvdj, has_vdjdata_col = has_vdjdata_col)


    # In[19]:


    #sc.pl.umap(adata, color=['vdj_all_productive', 'vdj_any_productive'])


    # Note that upon plotting, the pandas Series in the adata.obs dataframe are turned into a categorical type by scanpy. This can important, as you see below we compare with 'True' instead of True:

    # In[20]:


    n_prod = sum(adata.obs['vdj_all_productive'] == 'True')
    print('There are %d cells which have only productive chains.' % (n_prod))


    # In[21]:


    n_prod = sum(adata.obs['vdj_any_productive'] == 'True')
    print('There are %d cells with at least one productive chain.' % (n_prod))


    # In[22]:


    repl_values = {
        'True': 'False',
        'False': 'True',
        }
    adata.obs['vdj_none_productive'] = adata.obs['vdj_any_productive'].replace(to_replace=repl_values, inplace=False)
    n_prod = sum(adata.obs['vdj_none_productive'] == 'True')
    print('There are %d cells with no productive chains.' % (n_prod))


    # In[23]:


    sc.pl.umap(adata, color=['vdj_all_productive', 'vdj_any_productive','vdj_none_productive'], title=['All productive chains', 'Any productive chain', 'No productive chain'],
               frameon=True, save='_'+study+'_'+receptor+'_productive_chain.png')


    # The following command adds one metadata column for each type of chain found in the V(D)J data:

    # ## Clone count

    # The clone count (of a clonotype) is defined as the number of clones in the clonotype (as determined by Cell Ranger). We can flag which cells are members of a clonotype with not fewer than n clones:

    # In[27]:


    adata = helper.add_gt_n(adata, n=2, has_vdjdata_col=has_vdjdata_col)


    # In[28]:


    for condition in adata.obs[currentConditionColumn].unique().to_list():
        adata_sub = adata[adata.obs[currentConditionColumn]==condition]
        adata_sub = helper.add_obs(adata_sub, obs=['clone_count'], pyvdj = pyvdj, has_vdjdata_col = has_vdjdata_col)
        sc.pl.umap(adata_sub, color=['vdj_clone_count'], ncols = 1, use_raw=False, save = '_'+study+'_'+condition+'_'+receptor+'_clone_count_'+currentConditionColumn+'.png', 
               title = [condition])


    # #### Subset the anndata to cells having atleast one productive TRA chain and atleast one productive TRB chain and calculate the clone_counts again on this subset.

    # In[29]:


    adata_sub = adata[adata.obs[has_vdjdata_col] == True]
    adata_sub.obs[currentConditionColumn+'_sample'] = adata_sub.obs[currentConditionColumn+'_sample'].cat.remove_unused_categories()
    adata_sub.obs.loc[:, 'barcode_meta'] = [x+'_'+y for x,y in zip(adata_sub.obs.loc[:, 'barcode'], 
                                                                   adata_sub.obs.loc[:, 'sample'])]
    adata_sub.uns[pyvdj]['df'] = adata_sub.uns[pyvdj]['df'][(adata_sub.uns[pyvdj]['df']['barcode_meta']).isin(adata_sub.obs['barcode_meta'])]
    adata_sub.uns[pyvdj]['df'].shape

    #Subset cells
    include=[]
    for barcode_meta in adata_sub.uns[pyvdj]['df'].barcode_meta.unique().tolist():
        sub = adata_sub.uns[pyvdj]['df'].loc[adata_sub.uns[pyvdj]['df'].barcode_meta==barcode_meta, ['barcode_meta', 'chain']]
        if receptor=='TCR':
            isCellOfInterest=set(['TRA', 'TRB']).issubset(set(sub.chain.unique().tolist()))
        else:
            isCellOfInterest=(set(['IGH', 'IGK']).issubset(set(sub.chain.unique().tolist())) or set(['IGH', 'IGL']).issubset(set(sub.chain.unique().tolist())))
        if isCellOfInterest:
            include.append(barcode_meta)
    adata_sub.uns[pyvdj]['df'] = adata_sub.uns[pyvdj]['df'].loc[(adata_sub.uns[pyvdj]['df'].barcode_meta.isin(include)) 
                                                        & (adata_sub.uns[pyvdj]['df'].productive_any==True), :]
    adata_sub = adata_sub[(adata_sub.obs.barcode_meta.isin(include)) & (adata_sub.obs.vdj_any_productive=='True')]
    print(adata_sub.uns[pyvdj]['df'].shape)
    print(adata_sub.obs.shape)

    adata = helper.add_obs(adata_sub, obs=['clonotype', 'is_clone', 'clone_count'], pyvdj = pyvdj, has_vdjdata_col = has_vdjdata_col)


    # In[30]:


    #sc.settings.set_figure_params(dpi=300)
#     sc.pl.umap(adata, color=[has_vdjdata_col, 'vdj_is_clone'], 
#                title=[receptor+' has V(D)J data', 'Is clone'], 
#                palette=['grey', 'black'])


    # ### Clonotype distribution

    # Now we obtain statistics for each specified AnnData metadata category:

    # In[31]:


    clonotypes = adata.obs['vdj_clonotype_n_2'].unique()[2:].tolist()
    cdr3_dict = helper.get_spec(adata, clonotypes, pyvdj = pyvdj)


    # In[32]:


    meta = currentConditionColumn+'_sample'
    #~ adata.obs['metadata']  # one category for each 10x channel
    adata = helper.up_stats(adata, meta, pyvdj = pyvdj, has_vdjdata_col=has_vdjdata_col)


    # In[33]:


    adata.uns[pyvdj]['stats'][meta].keys()


    # In[34]:


    # Now we can plot cell numbers:
#     helper.plot_cells(adata.uns[pyvdj]['stats'][meta])
#     plt.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/has_vdj_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # Plot clonotype distribution for each sample:

    # In[35]:


    sc.settings.set_figure_params(dpi=300)
    clonotype_dist = adata.uns[pyvdj]['stats'][meta]['clonotype_dist']
    for meta_item in adata.obs[meta].unique().tolist():
        print(meta_item)
        clonotype_dist[meta_item].sort_index().plot(kind = 'bar', figsize=(7, 5))
        plt.subplots_adjust(bottom=0.3, left=0.2)
        plt.grid(False)
        plt.title(meta_item)
        plt.ylabel('Number of cells')
        plt.show()
        plt.close()


    # In[ ]:





    # ## Diversity indices

    # Which sample is the most diverse? We can make a simplified attempt to answer this by calculating the Shannon and the Simpson diversity indices.

    # ### Corrected for sample size

    # In[36]:


    adata_obs = adata.obs[adata.obs[has_vdjdata_col] == True]


    # In[37]:


    meta = currentConditionColumn+'_sample'
    category = 'vdj_clonotype'


    # In[38]:


    shannon_index_dict = helper.downsampling(adata_obs, meta, category)


    # In[39]:


    for s, i in shannon_index_dict.items():
        print(s)
        print(np.median(i))

    shi_dict = {}
    for s, i in shannon_index_dict.items():
        shi_dict[s] = np.median(i)

    shi_df = pd.DataFrame({
        'name': list(shi_dict),
        'shi': list(shi_dict.values())#,
        #'celltype': d #No idea what 'd' is here
        })
    #shi_df_all = pd.concat([shi_df_all, shi_df], ignore_index=True)


    # In[40]:


    adata.obs[currentConditionColumn+'_sample'].cat.categories


    # In[41]:


    shi_df.index = shi_df.name
    shi_df=shi_df.reindex(adata.obs[currentConditionColumn+'_sample'].cat.categories)
    shi_df.name = shi_df.index


    # Each sample was randomly downsampled to the size of the smallest sample 100 times, and the median of the 100 calculated indices is used.

    # In[42]:


    x_names = shi_df.name
    x_axis = range(1, x_names.shape[0]+1)  # bar positions on x: CHANGE THIS DEPENDING ON THE NUMBER OF SAMPLES IN EACH STUDY
    plt.bar(x=x_axis, height=shi_df.shi.to_list(), width=0.7)
    plt.xticks(x_axis, x_names, rotation='vertical')
    plt.grid(False)
    plt.ylabel('Shannon index ('+receptor+' diversity)')
    figure_div = plt.gcf()
    figure_div.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/ShI_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # In[43]:


    df = shi_df
    df[currentConditionColumn] = pd.Categorical(list(x.split('_')[0] for x in shi_df['name'].unique()),
                                    categories=adata.obs[currentConditionColumn].cat.categories,
                                    ordered=True)
    df.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/results/ShI_'+receptor+'_'+study+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)


    # In[44]:


    df


    # ## Expanded clonotypes

    # We can plot the proportion of expanded clonotypes. Expanded clonotypes are defined as having at least two clones. Alternatively, we could use the thresholds calculated by powerBCR, as mentioned above.

    # #### By sample

    # In[45]:


    n = 2  # clonotypes with at least 2 clones are considered expanded
    adata.obs['vdj_expanded_n_2'] = adata.obs['vdj_clone_count']
    adata.obs['vdj_expanded_n_2'][adata.obs['vdj_expanded_n_2'] < n] = 0
    adata.obs['vdj_expanded_n_2'][adata.obs['vdj_expanded_n_2'] >= n] = 1


    # In[46]:


    obs_vdj_df = adata.obs[adata.obs[has_vdjdata_col] == True]


    # In[47]:


    #obs_vdj_df[currentConditionColumn+'_sample'] = obs_vdj_df[currentConditionColumn+'_sample'].cat.remove_unused_categories()


    # In[48]:


    #obs_vdj_df.loc[:, 'barcode_meta'] = [x+'_'+y for x,y in zip(obs_vdj_df.loc[:, 'barcode'], obs_vdj_df.loc[:, 'sample'])]


    # In[49]:


    uns_vdj_df = adata.uns[pyvdj]['df'][(adata.uns[pyvdj]['df']['barcode_meta']).isin(obs_vdj_df['barcode_meta'])]
    uns_vdj_df.shape


    # In[51]:


    uns_vdj_df_complete = pd.merge(uns_vdj_df, obs_vdj_df.loc[:, ['barcode_meta', 'cell', currentConditionColumn+'_sample', currentConditionColumn]], how='outer', on='barcode_meta')


    # In[52]:


    uns_vdj_df_complete.productive_any.value_counts()


    # In[53]:


    #uns_vdj_df_complete


    # In[54]:


    obs_vdj_df.vdj_any_productive.value_counts()


    # In[55]:


    uns_vdj_df_complete.shape


    # In[56]:


    #uns_vdj_df_complete = uns_vdj_df_complete.loc[~uns_vdj_df_complete.chain.isin(['TRA', 'TRB']), :]


    # In[57]:


    exp_counts = obs_vdj_df.groupby([currentConditionColumn+'_sample']).vdj_expanded_n_2.value_counts()
    exp_counts.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/ExpandedProportions/expanded_proportions_'+receptor+'_'+study+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)


    # In[58]:


    exp_counts


    # In[59]:


    c_df = pd.read_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/ExpandedProportions/expanded_proportions_'+receptor+'_'+study+'_'+currentConditionColumn+'.tsv', sep='\t', header=0)


    # In[60]:


    if ((study=='Li')&(receptor=='BCR')):
        if (currentConditionColumn=='severity'):
            c_df = c_df.append({currentConditionColumn+'_sample':'Moderate_P5', 'vdj_expanded_n_2':1, 'vdj_expanded_n_2.1':0},
                               ignore_index=True).sort_values(by=[currentConditionColumn+'_sample'])
        else:
            c_df = c_df.append({currentConditionColumn+'_sample':'Convalescent_P5', 'vdj_expanded_n_2':1, 'vdj_expanded_n_2.1':0},
                               ignore_index=True).sort_values(by=[currentConditionColumn+'_sample'])
        c_df[currentConditionColumn+'_sample'] = c_df[currentConditionColumn+'_sample'].astype('category')
        c_df[currentConditionColumn+'_sample'] = c_df[currentConditionColumn+'_sample'].cat.reorder_categories(obs_vdj_df[currentConditionColumn+'_sample'].cat.categories)

    # In[61]:


    width = 0.5

    expanded_raw = c_df[c_df['vdj_expanded_n_2'] == 1]['vdj_expanded_n_2.1'].to_list()
    notexpan_raw = c_df[c_df['vdj_expanded_n_2'] == 0]['vdj_expanded_n_2.1'].to_list()
    total_raw = [x+y for x, y in zip(expanded_raw, notexpan_raw)]

    expanded = [x/y*100 for x, y in zip(expanded_raw, total_raw)]
    notexpan = [x/y*100 for x, y in zip(notexpan_raw, total_raw)]

    N = len(c_df[currentConditionColumn+'_sample'].unique())
    ind = np.arange(N)  # the x locations for the groups

#     p1 = plt.bar(ind, expanded, width)#, color=['#3B76AF', '#EF8536', '#EF8536'])
#     p2 = plt.bar(ind, notexpan, width, bottom=expanded, color='grey')
#     plt.grid(False)
#     plt.ylabel('Cells [%]')
#     # plt.title('')
#     plt.xticks(ind, c_df[currentConditionColumn+'_sample'].unique(), rotation='vertical')
#     # plt.legend((p1[0], p2[0]), ('Expanded', 'Not expanded'))
#     # plt.subplots_adjust(bottom=0.7, left=0.25, top=0.9)
#     figure_exp = plt.gcf()
#     figure_exp.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/expanded_proportions_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # This plot shows the proportion of cells that are clones of expanded clonotypes (grey: unique cells of not expanded clonotypes).

    # In[63]:


    df = pd.DataFrame({'names':c_df[currentConditionColumn+'_sample'].unique(),
                  'expanded':expanded,
                 'condition':(x.split('_')[0] for x in c_df[currentConditionColumn+'_sample'].unique()),
                      'type':'Expanded'})
    df.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/results/expanded_proportions_'+receptor+'_'+study+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)


    # In[64]:


    obs_vdj_df['vdj_clone_count_group'] = obs_vdj_df['vdj_clone_count']
    obs_vdj_df['vdj_clone_count_group'][obs_vdj_df['vdj_clone_count'] == 0] = "0"
    obs_vdj_df['vdj_clone_count_group'][obs_vdj_df['vdj_clone_count'] == 1] = "1"
    obs_vdj_df['vdj_clone_count_group'][obs_vdj_df['vdj_clone_count'] > 5] = ">5"
    obs_vdj_df['vdj_clone_count_group'][(obs_vdj_df['vdj_clone_count']>=2) & 
                                       (obs_vdj_df['vdj_clone_count']<=5)] = "2-5"


    # In[65]:


    meta = 'cell'
    df1 = pd.DataFrame(obs_vdj_df.groupby([meta, currentConditionColumn+'_sample', 'vdj_clone_count_group']).size(),
                 columns=['Cell count']).reset_index() 
    #print(df1.tail(10))
    tf = pd.DataFrame(df1.groupby([currentConditionColumn+'_sample'])['Cell count'].sum()).reset_index()
    #print(tf.tail(10))

    df = pd.merge(df1, tf,  how='left', on=[currentConditionColumn+'_sample'])
    #print(df)

    df['Cell count'] = (df['Cell count_x'].div(df['Cell count_y']))*100
    #print(df)

    df = df.loc[df.vdj_clone_count_group!='0']
    df['vdj_clone_count_group'] = df['vdj_clone_count_group'].astype('category')
    #print(df.head())

    df[currentConditionColumn] = pd.Categorical(list(x.split('_')[0] for x in df[currentConditionColumn+'_sample']),
                                    categories=adata.obs[currentConditionColumn].cat.categories,
                                    ordered=True)

    df.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/results/Cell_count_'+receptor+'_'+study+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)


    # In[66]:


    meta = 'cell'
    #df1 = pd.DataFrame(obs_vdj_df.groupby([meta, currentConditionColumn+'_sample', 'vdj_clone_count_group', 'vdj_clonotype']).size(),
    #             columns=['Cell count']).reset_index() 
    #print(df1.tail(10))
    #tf = pd.DataFrame(df1.groupby([meta, currentConditionColumn+'_sample'])['Cell count'].sum()).reset_index()
    #print(tf.tail(10))

    #meta = 'cell'
    df1 = pd.DataFrame(obs_vdj_df.groupby([meta, currentConditionColumn+'_sample', 'vdj_clone_count_group'])['vdj_clonotype'].nunique()).reset_index() 
    #print(df1.tail(10))
    tf = pd.DataFrame(df1.groupby([currentConditionColumn+'_sample'])['vdj_clonotype'].sum()).reset_index()
    #print(tf.tail(10))

    df = pd.merge(df1, tf,  how='left', on=[currentConditionColumn+'_sample'])
    #print(df)

    df['Cell count'] = (df['vdj_clonotype_x'].div(df['vdj_clonotype_y']))*100
    #print(df)

    df = df.loc[df.vdj_clone_count_group!='0']
    df['vdj_clone_count_group'] = df['vdj_clone_count_group'].astype('category')
    #print(df.head())

    df[currentConditionColumn] = list(x.split('_')[0] for x in df[currentConditionColumn+'_sample'])

    df.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/results/Clonotype_count_'+receptor+'_'+study+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)


    # In[67]:


    #Top 20 expanded clonotypes 
    #Top 20 expanded clonotypes 
    df1 = obs_vdj_df.loc[:, [currentConditionColumn+'_sample', 'cell', 'vdj_clone_count', 'vdj_clonotype']]
    tf = pd.DataFrame(df1.groupby([currentConditionColumn+'_sample'])['vdj_clone_count'].sum()).reset_index()

    df = pd.merge(df1, tf,  how='left', on=[currentConditionColumn+'_sample'])
    df['vdj_clone_count'] = (df['vdj_clone_count_x'].div(df['vdj_clone_count_y']))
    df[currentConditionColumn] = pd.Categorical(list(x.split('_')[0] for x in df[currentConditionColumn+'_sample']),
                                    categories=adata.obs[currentConditionColumn].cat.categories,
                                    ordered=True)

    #df = pd.DataFrame(df.groupby(['vdj_clonotype', 'condition', 'vdj_clone_count', 'cell']).size()).reset_index()
    for condition in df[currentConditionColumn].unique().tolist():
        sub = df.loc[df[currentConditionColumn]==condition, :]
        sub = sub.sort_values(by='vdj_clone_count', ascending=False)
        top_20_clonotypes = sub.vdj_clonotype.unique().tolist()[0:10]
        sub_sub = sub[sub.vdj_clonotype.isin(top_20_clonotypes)]
        sub_sub.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/Top20clonotypes/'+study+'_'+receptor+'_'+condition+'_'+currentConditionColumn+'.csv')
        #sns.barplot(x='vdj_clonotype', y='vdj_clone_count', data=sub_sub, hue='cell')
        #plt.show()


    # ## Public, private and condition-specific clonotypes

    # See definitions of public and private clonotypes.
    # 
    # Let's find public CDR3s! (i.e. CDR3s shared between donors or conditions)
    # There are two ways to approach this:
    # 
    # - We can obtain all CDR3 sequences for each condition, then find shared sequences. We provide an example in a section further below.
    # - A much stricter approach is to find CDR3-combinations (clonotypes) shared between conditions. This is more accurate, because the combination of two CDR3 sequences define the specificity of the BCR. One major disadvantage is that it will return very few positives. The implementation of this approach is imperfect here as it considers the set of CDR3s for each clonotype, instead of the pairs of CDR3s that make up a BCReceptor. For this method, we calculate statistics grouped by donor (as opposed to sample, i.e. 10x channel):

    # ### Finding public clonotypes

    # #### By sample

    # In[70]:


    meta = currentConditionColumn+'_sample'
    adata = helper.up_stats(adata, meta, pyvdj = pyvdj, has_vdjdata_col=has_vdjdata_col)
    cdr3 = adata.uns[pyvdj]['stats'][meta]['cdr3']


    # In[71]:


    #set(cdr3['C141']).intersection(*cdr3)


    # In[72]:


    cdr3_codes_rev = adata.uns[pyvdj]['cdr3']['cdr3_codes_rev']

    cdr3_coded = {}
    for key, value in cdr3.items():
        print(key)
        cdr3_coded[key] = [cdr3_codes_rev[cdr3] for cdr3 in value]

    cdr3_coded.keys()


    # We create a multi-indexed table that can be used for the plotting function:

    # In[73]:


    contents = upsetplot.from_contents(cdr3_coded)
    contents


    # In[74]:


    upsetplot.plot(contents, subset_size='count', show_counts=True)
    plt.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/upsetplot_public_clonotypes_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # #### By condition

    # In[75]:


    meta = currentConditionColumn
    adata = helper.up_stats(adata, meta, pyvdj = pyvdj, has_vdjdata_col=has_vdjdata_col)
    cdr3 = adata.uns[pyvdj]['stats'][meta]['cdr3']


    # In[76]:


    #set(cdr3['C141']).intersection(*cdr3)


    # In[77]:


    cdr3_codes_rev = adata.uns[pyvdj]['cdr3']['cdr3_codes_rev']

    cdr3_coded = {}
    for key, value in cdr3.items():
        print(key)
        cdr3_coded[key] = [cdr3_codes_rev[cdr3] for cdr3 in value]

    cdr3_coded.keys()


    # We create a multi-indexed table that can be used for the plotting function:

    # In[78]:


    contents = upsetplot.from_contents(cdr3_coded)
    contents


    # In[79]:


    upsetplot.plot(contents, subset_size='count', show_counts=True)
    plt.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/upsetplot_public_clonotypes_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # ### Finding public CDR3 sequences

    # We can use get_spec() with a list of clonotypes for each category, or obtain sequences directly from the data, filtered for cells as shown here:

    # #### By sample

    # In[80]:


    meta = currentConditionColumn+'_sample'
    vdjdf = adata.uns[pyvdj]['df']
    cdr3_simple_dict = {}


    # In[81]:


    obs_col = adata.uns[pyvdj]['obs_col']
    for m in adata.obs[meta].unique():
        print(m)
        cells = adata.obs.loc[adata.obs[meta] == m][obs_col].tolist()
        cdr3_m = vdjdf.loc[vdjdf['barcode_meta'].isin(cells)]['cdr3']
        cdr3_simple_dict[m] = set(cdr3_m)


    # In[82]:


    cdr3_simple_dict.keys()


    # In[83]:


    #!pip install venn
    #from venn import venn
    #%matplotlib inline

    #venn(cdr3_simple_dict)


    # In[84]:


    contents = upsetplot.from_contents(cdr3_simple_dict)
    contents
    #contents.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/upsetplot_public_clonotypes_2_'+receptor+'_'+study+'_'+currentConditionColumn+'.csv')


    # In[85]:


    #upsetplot.plot(contents, subset_size='count', show_counts=True)
    #plt.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/upsetplot_public_cdr3_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # In[86]:


    #We can list the public CDR3s (i.e. the ones found in more than one donor) with the following code:
    allcdr3 = set()
    for k, v in cdr3_simple_dict.items():
        allcdr3 = allcdr3 | v


    # In[87]:


    len(allcdr3)


    # We have made a set containing all CDR3s and now we list which donor each CDR3 can be found in.

    # In[88]:


    public_cdr3_table = pd.DataFrame(index=allcdr3)
    for k, v in cdr3_simple_dict.items():
        bool_cdr3 = [True if cdr3_val in v else False for cdr3_val in allcdr3]
        public_cdr3_table[k] = bool_cdr3


    # In[89]:


    public_cdr3_table['sum'] = public_cdr3_table.sum(axis=1) > 1


    # In[90]:


    print('There are %d public CDR3s in the dataset.' % sum(public_cdr3_table['sum']))


    # To be thorough in finding public CDR3s, we also compare with public databases: download and unzip this file: https://github.com/antigenomics/vdjdb-db/releases/tag/2019-08-08, then

    # In[91]:


    #!unzip -q /nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj//nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj//nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/vdjdb-2019-08-08.zip -d /nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/vdjdb-2019-08-08


    # In[92]:


    vdjdb_file = '/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/vdjdb-2019-08-08/vdjdb.txt'
    vdjdb = pd.read_csv(vdjdb_file, sep='\t')


    # In[93]:


    public_vdjdb = allcdr3 & set(vdjdb['cdr3'])


    # In[94]:


    print('There are %i CDR3s in our data that were registered in VDJdb.' % len(public_vdjdb))


    # In[95]:


    len(public_vdjdb & set(public_cdr3_table[public_cdr3_table['sum']].index))


    # In[96]:


    public_cdr3_table['in_vdjdb'] = 'False'
    public_cdr3_table.loc[list(public_vdjdb & set(public_cdr3_table[public_cdr3_table['sum']].index)), 'in_vdjdb'] = 'True'  


    # In[97]:


    public_cdr3_table[public_cdr3_table['sum']].to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/PublicCDR3Sequences/public_cdr3_table_sample_'+receptor+'_'+study+'_'+currentConditionColumn+'.csv')


    # ## CDR3 sequences to align

    # Pick the public CDR3 sequences shared between different conditions (except healthy). Instead of finding the exact same CDR-H3 sequences as above, the criteria below ais to find the CDR-H3 sequences with Hamming distance of less than, equal to 15% of the length of the CDR-H3 sequence
    # Follow the paper: https://www.biorxiv.org/content/10.1101/2020.07.08.194456v1.full.pdf
    # 
    # 460 IGH sequences from single cells with paired productive heavy and light chains were searched\
    # 461 against COVID-19 patient bulk IGH repertoires to identify convergent sequences according to\
    # 462 the following criteria:
    # 
    # - utilization of the same IGHV and IGHJ genes; 
    # - same CDR-H3 lengths;
    # - CDR-H3 amino acid sequences that were within a Hamming distance cutoff of 15% of the length of the CDR-H3. 
    # 
    # 464 Two native heavy and light chain pairs, designated mAb2A and mAb4A, which were found in convergent clusters characterized by 
    # - low-to-mid SHM frequencies
    # - included at least one class-switched member
    # 
    # were selected for cloning and expression.

    # In[98]:


    #keep all the cdr3 sequences, not in the healthy samples
    #healthy_cols = uns_vdj_df_complete[uns_vdj_df_complete[currentConditionColumn]=="Healthy"][currentConditionColumn+'_sample'].unique().tolist()
    #df = public_cdr3_table[public_cdr3_table[healthy_cols].sum(axis=1)==0]
    #cdr3_seqs_to_include = df[df['sum']].index.tolist()

    if receptor=='TCR':
        chains = ['TRA', 'TRB']
    else:
        chains = ['IGH', 'IGK', 'IGL']

    for chain in chains:
        df1 = uns_vdj_df_complete.loc[((uns_vdj_df_complete.cdr3!="None")&
                                      (uns_vdj_df_complete.v_gene!="None")&
                                      (uns_vdj_df_complete.j_gene!="None"))&
                                      (uns_vdj_df_complete.chain==chain),
                                      ['v_gene', 'j_gene', 'd_gene', 'cdr3', 'clonotype_meta', currentConditionColumn,
                                         'cdr3_nt', currentConditionColumn+'_sample']]
        df1['cdr3_len'] = df1.cdr3.apply(func=len)

        df1.drop_duplicates(inplace=True)

        df=pd.DataFrame(columns = ['v_gene', 'j_gene', 'd_gene',
                                   'cdr3_len', 'hamming_distance', 'cutOff', 
                                   'cdr3', 'cdr3_nt', 'condition','clonotype_meta',
                                  'sample', 'comparison'])

        #df1 = df1.iloc[0:5000, :]

        hdThresh = 0.15

        for index, row in df1.loc[:,['v_gene','j_gene','cdr3_len']].drop_duplicates().iterrows():
            #print(row)
            v_gene = row['v_gene'] 
            j_gene = row['j_gene'] 
            cdr3_len = row['cdr3_len']
            sub = df1[(df1.v_gene==v_gene)&(df1.j_gene==j_gene)&(df1.cdr3_len==cdr3_len)]
            if len(sub.cdr3.unique().tolist())>1:
                #if((len(sub[currentConditionColumn].unique())>1) & (sub[currentConditionColumn].unique().tolist() not in ['Healthy'])):
                cdr3_seq_pairs = list(itertools.combinations(sub.cdr3.unique().tolist(), 2)) 
                for cdr3_seq_pair in cdr3_seq_pairs:
                    h_dist = helper.hamming_distance(cdr3_seq_pair[0], cdr3_seq_pair[1])
                    cutOff = hdThresh*cdr3_len
                    if (h_dist <= cutOff):
                        #print(cdr3_seq_pair)
                        temp1 = {'v_gene':v_gene,
                               'j_gene':j_gene,
                                 'd_gene': tuple(sub[(sub.v_gene==v_gene) & (sub.j_gene==j_gene) & (sub.cdr3==cdr3_seq_pair[0])].d_gene.unique().tolist()),
                                'cdr3_len':cdr3_len,
                                'hamming_distance': h_dist,
                                'cutOff': cutOff,
                                 'cdr3':cdr3_seq_pair[0],
                                'cdr3_nt':tuple(sub[sub.cdr3==cdr3_seq_pair[0]].cdr3_nt.unique().tolist()),
                                'condition':tuple(sub[sub.cdr3==cdr3_seq_pair[0]][currentConditionColumn].unique().tolist()),
                                'clonotype_meta':tuple(sub[sub.cdr3==cdr3_seq_pair[0]].clonotype_meta.unique().tolist()),
                                'sample':tuple(sub[sub.cdr3==cdr3_seq_pair[0]][currentConditionColumn+'_sample'].unique().tolist()),
                                'comparison':tuple(sub[sub.cdr3.isin(cdr3_seq_pair)][currentConditionColumn+'_sample'].unique().tolist())}

                        temp2 = {'v_gene':v_gene,
                                 'j_gene':j_gene,
                                 'd_gene': tuple(sub[(sub.v_gene==v_gene) & (sub.j_gene==j_gene) & (sub.cdr3==cdr3_seq_pair[1])].d_gene.unique().tolist()),
                                'cdr3_len':cdr3_len,
                                'hamming_distance': h_dist,
                                'cutOff': cutOff,
                                'cdr3':cdr3_seq_pair[1], 
                                'cdr3_nt':tuple(sub[sub.cdr3==cdr3_seq_pair[1]].cdr3_nt.unique().tolist()),
                                'condition':tuple(sub[sub.cdr3==cdr3_seq_pair[1]][currentConditionColumn].unique().tolist()),
                                'clonotype_meta':tuple(sub[sub.cdr3==cdr3_seq_pair[1]].clonotype_meta.unique().tolist()),
                                 'sample':tuple(sub[sub.cdr3==cdr3_seq_pair[1]][currentConditionColumn+'_sample'].unique().tolist()),
                                 'comparison':tuple(sub[sub.cdr3.isin(cdr3_seq_pair)][currentConditionColumn+'_sample'].unique().tolist())}


                        #Remove this entry if this sequence is present in healthy samples
                        if (set(['Healthy']).issubset(set(list(temp1['condition']))))|(set(['Healthy']).issubset(set(list(temp2['condition'])))):
                            pass #shared with healthy samples
                        elif (temp1['sample']==temp2['sample']):
                            pass #not a public cdr3 sequence
                        else:
                            df = df.append(temp1, ignore_index=True)
                            #print('temp1')
                            df = df.append(temp2, ignore_index=True)
                            #print('temp2')
                            #write to a text file for alignment

        #print(df)
        #Remove duplicates from the dataframe
        fileName = '/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/PublicCDR3Sequences/CDR3_seq_'+study+'_'+receptor+'_'+chain+'_'+str(hdThresh)
        df.drop_duplicates(inplace=True)
        df.to_csv(fileName+'_'+currentConditionColumn+'.csv')

        #Write CDR3 amino acid sequence to a file
        if os.path.exists(fileName+'_aa.txt'):
            os.remove(fileName+'_aa.txt')

        keys=['v_gene', 'j_gene', 'd_gene', 'clonotype_meta']
        for current_len in df.cdr3_len.unique():
            f=open(fileName+'_'+str(current_len)+'_aa.txt', 'a')
            for index, row in df.loc[df.cdr3_len==current_len, :].iterrows():
                f.write('>'+'_'.join([str(elem) for elem in list(itertools.chain(*[list(row[k] for k in keys), [study], [index]]))])+"\n"+row['cdr3']+"\n")
            f.close()
            #Writte CDR3 nucleotide sequence to a file
            if os.path.exists(fileName+'_nt.txt'):
                os.remove(fileName+'_nt.txt')
            f=open(fileName+'_'+str(current_len)+'_nt.txt', 'a')
            for index, row in df.loc[df.cdr3_len==current_len, :].iterrows():
                iteration=0
                for i in range(0,len(row['cdr3_nt'])):
                    f.write('>'+'_'.join([str(elem) for elem in list(itertools.chain(*[list(row[k] for k in keys), [study], [index], [iteration]]))])+"\n")
                    f.write(row['cdr3_nt'][i]+"\n")
                    iteration=iteration+1
            f.close()


    # In[99]:


    #sub


    # #### By condition (to include in paper)

    # In[100]:


    meta = currentConditionColumn
    vdjdf = adata.uns[pyvdj]['df']
    cdr3_simple_dict = {}


    # In[101]:


    obs_col = adata.uns[pyvdj]['obs_col']
    for m in adata.obs[meta].unique():
        print(m)
        cells = adata.obs.loc[adata.obs[meta] == m][obs_col].tolist()
        cdr3_m = vdjdf.loc[vdjdf['barcode_meta'].isin(cells)]['cdr3']
        cdr3_simple_dict[m] = set(cdr3_m)


    # In[102]:


    cdr3_simple_dict.keys()


    # In[103]:


    #!pip install venn
    #from venn import venn
    #%matplotlib inline

    #venn(cdr3_simple_dict)


    # In[104]:


    contents = upsetplot.from_contents(cdr3_simple_dict)
    contents
    #contents.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/upsetplot_public_clonotypes_2_'+receptor+'_'+study+'_'+currentConditionColumn+'.csv')


    # In[105]:


    upsetplot.plot(contents, subset_size='count', show_counts=True)
    plt.savefig('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/scripts/figures/upsetplot_public_cdr3_'+receptor+'_'+study+'_'+currentConditionColumn+'.png')


    # In[106]:


    #We can list the public CDR3s (i.e. the ones found in more than one donor) with the following code:
    allcdr3 = set()
    for k, v in cdr3_simple_dict.items():
        allcdr3 = allcdr3 | v


    # In[107]:


    len(allcdr3)


    # We have made a set containing all CDR3s and now we list which donor each CDR3 can be found in.

    # In[108]:


    public_cdr3_table = pd.DataFrame(index=allcdr3)
    for k, v in cdr3_simple_dict.items():
        bool_cdr3 = [True if cdr3_val in v else False for cdr3_val in allcdr3]
        public_cdr3_table[k] = bool_cdr3


    # In[109]:


    public_cdr3_table['sum'] = public_cdr3_table.sum(axis=1) > 1


    # In[110]:


    print('There are %d public CDR3s in the dataset.' % sum(public_cdr3_table['sum']))


    # To be thorough in finding public CDR3s, we also compare with public databases: download and unzip this file: https://github.com/antigenomics/vdjdb-db/releases/tag/2019-08-08, then

    # In[111]:


    #!unzip -q /nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj//nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj//nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/vdjdb-2019-08-08.zip -d /nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/vdjdb-2019-08-08


    # In[112]:


    vdjdb_file = '/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/vdjdb-2019-08-08/vdjdb.txt'
    vdjdb = pd.read_csv(vdjdb_file, sep='\t')


    # In[113]:


    public_vdjdb = allcdr3 & set(vdjdb['cdr3'])


    # In[114]:


    print('There are %i CDR3s in our data that were registered in VDJdb.' % len(public_vdjdb))


    # In[115]:


    len(public_vdjdb & set(public_cdr3_table[public_cdr3_table['sum']].index))


    # In[116]:


    public_cdr3_table['in_vdjdb'] = 'False'
    public_cdr3_table.loc[list(public_vdjdb & set(public_cdr3_table[public_cdr3_table['sum']].index)), 'in_vdjdb'] = 'True'  


    # In[117]:


    public_cdr3_table[public_cdr3_table['sum']].to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/PublicCDR3Sequences/public_cdr3_table_sample_'+receptor+'_'+study+'_'+currentConditionColumn+'.csv')


    # ## VDJ chain usage

    # In[118]:

    if receptor=='TCR':
        chains = ['TRA', 'TRB']
    else:
        chains = ['IGH', 'IGK', 'IGL']
        
    for x_gene in ['v_gene', 'j_gene']:
        for chain in chains:
    #x_gene = 'v_gene'
    #chain = 'TRA'
            sub = uns_vdj_df_complete[uns_vdj_df_complete.chain==chain]
            print(sub.shape)
            df1 = pd.DataFrame(sub.groupby([currentConditionColumn+'_sample', x_gene]).size(),
                         columns=['Gene count']).reset_index() 
            #print(df1.head())
            #print(df1.tail())

            tf = pd.DataFrame(df1.groupby([currentConditionColumn+'_sample'])['Gene count'].sum()).reset_index()
            #print(tf.head())
            #print(tf.tail())

            df = pd.merge(df1, tf,  how='left', on=[currentConditionColumn+'_sample'])
            #print(df)

            df['Gene count'] = (df['Gene count_x'].div(df['Gene count_y']))*100
            #print(df)

            df = df.loc[df[x_gene]!='None']
            df[x_gene] = df[x_gene].astype('category')
            #print(df.head())

            df[currentConditionColumn] = pd.Categorical(list(x.split('_')[0] for x in df[currentConditionColumn+'_sample']),
                                    categories=adata.obs[currentConditionColumn].cat.categories,
                                    ordered=True)
            df.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/results/boxplot_'+receptor+'_'+study+'_'+chain+'_'+x_gene+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)


    # #### Export data frames for plotting chain usage heatmap in R

    # In[119]:


    df1 = pd.DataFrame(uns_vdj_df_complete.groupby([currentConditionColumn+'_sample', 'v_gene', 'j_gene']).size(), columns=['Freq']).reset_index() 
    #print(df1.head())
    #print(df1.tail())

    tf = pd.DataFrame(df1.groupby([currentConditionColumn+'_sample'])['Freq'].sum()).reset_index()
    #print(tf.head())
    #print(tf.tail())

    df = pd.merge(df1, tf,  how='left', on=[currentConditionColumn+'_sample'])
    #print(df)

    df['Freq'] = (df['Freq_x'].div(df['Freq_y']))*100
    #print(df)
    df[currentConditionColumn] = pd.Categorical(list(x.split('_')[0] for x in df[currentConditionColumn+'_sample']),
                                    categories=adata.obs[currentConditionColumn].cat.categories,
                                    ordered=True)

    df.to_csv('/nfs/leia/research/ma/manikg/COVID-19sc/meta-analysis/vdj/data/ChainUsage/'+receptor+'/ChainUsage_'+study+'_'+receptor+'_'+currentConditionColumn+'.tsv', sep='\t', header=True)

    with open(outFile, 'w') as fp: 
        print('Done!')

