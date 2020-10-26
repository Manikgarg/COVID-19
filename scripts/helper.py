#
# This file is part ofpv.
#
# pyVDJ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyVDJ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along withpv.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
import random
import pyvdj as pv
from scipy import stats
import igraph
import re
import itertools

def find_matches(item, d):
    for k in d:
        if re.search(k, item):
            return re.sub(k, d[k], item)


def hamming_distance(chain1, chain2):
    #Return the Hamming distance between equal-length sequences
    if len(chain1) != len(chain2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(c1 != c2 for c1, c2 in zip(chain1, chain2))


def add_has_vdjdata(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    df = adata.uns[pyvdj]['df']

    adata.obs[has_vdjdata_col] = adata.obs[obs_col].isin(df['barcode_meta'])
      # 'barcode_meta' is in the df and corresponds to 'obs_col' in anndata.
    return adata


def add_clonotype(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    df = adata.uns[pyvdj]['df']

    adata.obs['vdj_clonotype'] = adata.obs[obs_col]
    # Remove missing ones:
    adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_clonotype'] = None

    BCR_dict = dict(zip(df['barcode_meta'], df['clonotype_meta']))
    adata.obs['vdj_clonotype'].replace(to_replace=BCR_dict, inplace=True)
    adata.obs['vdj_clonotype'] = adata.obs['vdj_clonotype'].astype('category')

    return adata


def add_is_clone(adata, pyvdj, has_vdjdata_col):
    if 'vdj_clone_count' not in adata.obs.columns:
        adata = add_clone_count(adata, pyvdj, has_vdjdata_col)
    adata.obs['vdj_is_clone'] = adata.obs['vdj_clone_count'] > 1

    return adata


def add_all_productive(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    df = adata.uns[pyvdj]['df']

    adata.obs['vdj_all_productive'] = adata.obs[obs_col]
    adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_all_productive'] = None
    prod_dict = dict(zip(df['barcode_meta'], df['productive_all']))
    adata.obs['vdj_all_productive'].replace(to_replace=prod_dict, inplace=True)

    return adata


def add_any_productive(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    df = adata.uns[pyvdj]['df']

    adata.obs['vdj_any_productive'] = adata.obs[obs_col]
    adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_any_productive'] = None
    prod_dict = dict(zip(df['barcode_meta'], df['productive_any']))
    adata.obs['vdj_any_productive'].replace(to_replace=prod_dict, inplace=True)

    return adata


def add_chains(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    chains = adata.uns[pyvdj]['df']['chain'].unique()
    chains = [x for x in chains if str(x) != 'nan']

    chain_nested_dict = dict() # for making one .obs column for each chain type
    grouped_cells = adata.uns[pyvdj]['df'].groupby('barcode_meta')
    for c in chains:
        print(c)
        adata.uns[pyvdj]['df']['chain_' + c] = adata.uns[pyvdj]['df']['chain'] == c
        has_chain = grouped_cells['chain_' + c].apply(any)
        chain_nested_dict[c] = dict(zip(has_chain.index, has_chain))

    for c in chain_nested_dict.keys():
        adata.obs['vdj_chain_' + c] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_chain_' + c] = 'No_data'

        adata.obs['vdj_chain_' + c].replace(to_replace=chain_nested_dict[c], inplace=True)

    return adata


def add_genes(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    constant_genes = adata.uns[pyvdj]['df']['c_gene'].unique()
    constant_genes = [x for x in constant_genes if str(x) != 'nan']

    constant_genes_nested_dict = dict()
    grouped_cells = adata.uns[pyvdj]['df'].groupby('barcode_meta')
    for c in constant_genes:
        print(c)
        adata.uns[pyvdj]['df']['vdj_constant_' + c] = adata.uns[pyvdj]['df']['c_gene'] == c
        has_c_gene = grouped_cells['vdj_constant_' + c].apply(any)
        constant_genes_nested_dict[c] = dict(zip(has_c_gene.index, has_c_gene))

    for c in constant_genes_nested_dict.keys():
        print('Preparing annotation for %s' % c)
        adata.obs['vdj_constant_' + c] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_constant_' + c] = 'No_data'
        adata.obs['vdj_constant_' + c].replace(to_replace=constant_genes_nested_dict[c], inplace=True)

    return adata


def add_v_genes(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    v_genes = adata.uns[pyvdj]['df']['v_gene'].unique()
    v_genes = [x for x in v_genes if str(x) != 'nan']

    v_genes_nested_dict = dict()
    grouped_cells = adata.uns[pyvdj]['df'].groupby('barcode_meta')
    for v in v_genes:
        print(v)
        adata.uns[pyvdj]['df']['vdj_v_' + v] = adata.uns[pyvdj]['df']['v_gene'] == v
        has_v_gene = grouped_cells['vdj_v_' + v].apply(any)
        v_genes_nested_dict[v] = dict(zip(has_v_gene.index, has_v_gene))

    for v in v_genes_nested_dict.keys():
        adata.obs['vdj_v_' + v] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_v_' + v] = 'No_data'
        adata.obs['vdj_v_' + v].replace(to_replace=v_genes_nested_dict[v], inplace=True)

    return adata


def add_j_genes(adata, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    j_genes = adata.uns[pyvdj]['df']['j_gene'].unique()
    j_genes = [x for x in j_genes if str(x) != 'nan']

    j_genes_nested_dict = dict()
    grouped_cells = adata.uns[pyvdj]['df'].groupby('barcode_meta')
    for j in j_genes:
        print(j)
        adata.uns[pyvdj]['df']['vdj_j_' + j] = adata.uns[pyvdj]['df']['j_gene'] == j
        has_j_gene = grouped_cells['vdj_j_' + j].apply(any)
        j_genes_nested_dict[j] = dict(zip(has_j_gene.index, has_j_gene))

    for j in j_genes_nested_dict.keys():
        adata.obs['vdj_j_' + j] = adata.obs[obs_col]
        adata.obs.loc[(-adata.obs[has_vdjdata_col]), 'vdj_j_' + j] = 'No_data'
        adata.obs['vdj_j_' + j].replace(to_replace=j_genes_nested_dict[j], inplace=True)

    return adata


def add_clone_count(adata, pyvdj, has_vdjdata_col):
    # Number of clones in clonotype
    if 'vdj_clonotype' not in adata.obs.columns:
        adata = add_clonotype(adata, pyvdj, has_vdjdata_col)

    clone_count = adata.obs['vdj_clonotype'].value_counts()
    clone_count_dict = dict(zip(clone_count.index, clone_count))

    adata.obs['vdj_clone_count'] = adata.obs['vdj_clonotype']
    adata.obs['vdj_clone_count'].replace(to_replace=clone_count_dict, inplace=True)

    adata.obs.loc[(adata.obs[has_vdjdata_col] == False), 'vdj_clone_count'] = 0
      # no data
    adata.obs['vdj_clone_count'] = adata.obs['vdj_clone_count'].astype(int)

    return adata


def add_obs(adata, obs, pyvdj, has_vdjdata_col):
    # obs: list. which of the below metadata to add?
    adder_functions = {
        'has_vdjdata': add_has_vdjdata,  # load_vdj() adds this by default
        'clonotype': add_clonotype,
        'is_clone': add_is_clone,
        'all_productive': add_all_productive,
        'any_productive': add_any_productive,
        'chains': add_chains,
        'genes': add_genes,
        'v_genes': add_v_genes,
        'j_genes': add_j_genes,
        'clone_count': add_clone_count,
        }

    for e in obs:
        func = adder_functions[e]
        adata = func(adata, pyvdj, has_vdjdata_col)

    return adata

def add_gt_n(adata, n, has_vdjdata_col):
    # Clonotypes with not less than n clones
    if (not isinstance(n, int) or n < 0):
        raise ValueError('n must be a non-negative integer')

    # These two adata.obs columns will be used:
    if 'vdj_clonotype' not in adata.obs.columns:
        raise Exception('Please add \'clonotype\' with add_obs()')
    if 'vdj_clone_count' not in adata.obs.columns:
        raise Exception('Please add \'clone_count\' with add_obs()')

    clon_bool = adata.obs['vdj_clonotype'].value_counts()
    txt_less = 'Less_than_n'
    clon_bool[clon_bool < n] = txt_less
    clon_bool[clon_bool != txt_less] = clon_bool[clon_bool != txt_less].index.tolist()
    clon_dict = dict(zip(clon_bool.index, clon_bool))

    newcol = 'vdj_clonotype_n_' + str(n)
    adata.obs[newcol] = adata.obs['vdj_clonotype']
    adata.obs[newcol] = adata.obs[newcol].cat.add_categories(['No_data', txt_less])
    adata.obs.loc[(-adata.obs[has_vdjdata_col]), newcol] = 'No_data'
    adata.obs[newcol].replace(to_replace=clon_dict, inplace=True)

    return adata

def get_spec(adata, clonotypes, pyvdj):
    df = adata.uns[pyvdj]['df']
    clonotypes = list(set(clonotypes))

    cdr3_dict = {}
    for c in clonotypes:
        print(c)
        cdr3 = df.loc[df['clonotype_meta'] == c, 'cdr3']
        cdr3 = cdr3.unique()
        cdr3 = cdr3.tolist()
        try:
            cdr3.remove('None')
        except ValueError:
            pass
        cdr3_dict[c] = cdr3

    return cdr3_dict

def make_cdr3set(df):  # for finding public CDR3s
    # df is adata.uns[pyvdj]['df']
    # Creates a set of productive CDR3s for each cell in df
    df = df[df['productive']]
    cdr3sets = df.groupby('barcode_meta')['cdr3'].apply(frozenset)

    cdr3_values = cdr3sets.unique()
    cdr3_codes = dict(zip(range(0, len(cdr3_values) ), cdr3_values))
    cdr3_codes_rev = dict(zip(cdr3_values, range(0, len(cdr3_values) )))

    cdr3set_dict = {
        'cdr3sets': cdr3sets,
        'cdr3_codes': cdr3_codes,
        'cdr3_codes_rev': cdr3_codes_rev}

    return cdr3set_dict


def up_stats(adata, meta, pyvdj, has_vdjdata_col):
    obs_col = adata.uns[pyvdj]['obs_col']
    # This function returns a dictionary of statistics on the VDJ data,
    # grouped by categories in the adata.obs[meta] column. Keys:
    # 'meta' stores the adata.obs columname
    # 'cells' count of cells, and cells with VDJ data per category
    # 'clonotype_counts' number of different clonotypes per category
    # 'clonotype_dist' clone count distribution
    # 'shared_cdr3' dictionary of cdr3 - cell
    stats_dict = {}
    stats_dict['meta'] = meta

    n_cells = adata.obs.groupby(meta).size()
    n_hasvdj = adata.obs.groupby(meta)[has_vdjdata_col].apply(sum)
    stats_dict['cells'] = [n_cells, n_hasvdj]

    n_cl = adata.obs.groupby(meta)['vdj_clonotype'].unique()
    n_cl = n_cl.apply(lambda x: ([z for z in x if str(z) != 'nan']))
    dict_cl = dict(zip(adata.obs[meta].unique(), n_cl))
    stats_dict['clonotype_unique'] = dict_cl
    len_cl = [len(y) for y in n_cl]
    dict_cl = dict(zip(adata.obs[meta].unique(), len_cl))
    stats_dict['clonotype_counts'] = dict_cl

    dist_cl = adata.obs.groupby(meta)['vdj_clone_count'].value_counts()
    stats_dict['clonotype_dist'] = dist_cl

    if 'cdr3' not in adata.uns[pyvdj]:
        adata.uns[pyvdj]['cdr3'] = make_cdr3set(adata.uns[pyvdj]['df'])
        print('Adding CDR3 to adata.uns')

    # 'shared_cdr3'
    cdr3set_dict = adata.uns[pyvdj]['cdr3']
    cdr3sets = cdr3set_dict['cdr3sets']
    # Group anndata cells (obs_col) by meta cats
    grouped = adata.obs.groupby(meta)
    # For each meta cat cellist: subset cdr3df to these cells, then
    cdr3_dict = {}
    for name, group in grouped:
        group[obs_col].values
        cdr3_list = cdr3sets[cdr3sets.index.isin(group[obs_col].values)]

        cdr3_dict[name] = cdr3_list.unique()  # keys of cdr3_dict
# and their copynumber (cdr3 clones) within groups -- postponed for later
        # Return simplified sets:
        cdr3_dict[name] = cdr3_list.unique()

    stats_dict['cdr3'] = cdr3_dict

    if 'stats' not in adata.uns[pyvdj]:
        adata.uns[pyvdj]['stats'] = {}
    adata.uns[pyvdj]['stats'][meta] = stats_dict

    return adata

# We define a plotting function. For easy customization, it is not included in the package.
def plot_cells(stats_dict):
    hasvdj = stats_dict['cells'][1]
    novdj = stats_dict['cells'][0] - hasvdj

    bars = np.add(hasvdj, novdj).tolist()
     
    x_axis = list(range(0, len(bars)))  # bar positions on x
    x_names = list(hasvdj.index)
    barWidth = 1

    p1 = plt.bar(x_axis, hasvdj, color='#553864', edgecolor='white', width=barWidth)
    p2 = plt.bar(x_axis, novdj, bottom=hasvdj, color='#81a75d', edgecolor='white', width=barWidth)
    plt.legend((p1[0], p2[0]), ('Has VDJ', 'No VDJ'), ncol=2, bbox_to_anchor=(0, 1), loc='lower left', fontsize='small')
    plt.grid(False)
    plt.xticks(x_axis, x_names, rotation='vertical')
    plt.subplots_adjust(bottom=0.4, left=0.2)
    plt.xlabel(stats_dict['meta'])
    plt.ylabel('Number of cells')
    
    
def downsampling(collection_df, meta, category, cutoff=0):
    # collection_df: df of annotation (celltype, species) and sample (donor, region) for each measurement (cell, specimen)
    # meta: group by this column of collection_df
    # category: count by this column of collection_df
    collection_dict = dict()
    subset_lengths = []
    for s in collection_df[meta].unique():
        subset_s = collection_df.loc[collection_df[meta] == s][category]

        s_length = len(subset_s)
        if s_length > cutoff:
            subset_lengths.append(s_length)
            collection_dict[s] = subset_s.astype(str)
        else:
            print('Sample "%s" skipped (count too low: %s)' % (s, s_length))

    shannon_index_dict = dict()
    x = min(subset_lengths)

    for s, subset_s in collection_dict.items():
        shannon_index_dict[s] = []

        l = len(subset_s)
        print()
        print('Subsampling sample %s' % s)
        print()
        for r in range(0, 100):
            print('Iteration %d ' % r, end='\r')
            x_rand_index = random.sample(range(l), x)
            subset_s_x = subset_s.iloc[x_rand_index]
            sd = dict(subset_s_x.value_counts())
            shannon_index =pv.shannon(sd)
            shannon_index_dict[s].append(shannon_index)

        plt.hist(shannon_index_dict[s])
        plt.title(s)
        plt.xlabel('Shannon index')
        plt.show()
        plt.close()

    return shannon_index_dict

try:
    import numpy as np
    from Levenshtein import distance
    from scipy.spatial.distance import pdist, squareform
except ImportError:
    has_pkgs = False
else:
    has_pkgs = True

try:
    import igraph
except ImportError:
    has_igraph = False
else:
    has_igraph = True


def get_cdr3(adata, pyvdj):
    df = adata.uns[pyvdj]['df']
    samples = adata.obs['sample'].unique()

    cdr3_dict = {}
    for s in samples:
        cdr3 = df.loc[(df['sample'] == s), 'cdr3']
        cdr3 = cdr3.unique()
        cdr3 = cdr3.tolist()
        try:
            cdr3.remove('None')
        except ValueError:
            pass
        cdr3_dict[s] = cdr3

    return cdr3_dict


def get_dist(cdr3_dict, sample):
    """Requires numpy, scipy and Levenshtein."""
    if not has_pkgs:
        raise ImportError("Function requires numpy, scipy and python-Levenshtein packages.")

    cdr3 = cdr3_dict[sample]
    cdr3 = np.array(cdr3).reshape(-1, 1)
    dist = pdist(cdr3, lambda x,y: distance(x[0], y[0]))

    dist = (squareform(dist))

    return dist


def graph_cdr3(dist):
    """Requires igraph."""
    if not has_igraph:
        raise ImportError("Function requires the igraph-python package.")

    dist[dist > 1] = 0 # use Levenshtein distance of 1
    g = igraph.Graph.Adjacency((dist > 0).tolist())
    
    return g