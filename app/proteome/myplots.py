import base64
from io import BytesIO
from itertools import  groupby
import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from venn import venn

from  matplotlib.colors import LinearSegmentedColormap

from matplotlib.lines import Line2D
import matplotlib

# import ptitprince as pt

from upsetplot import from_contents
from upsetplot import plot as UPSET_PLOT

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from scipy.cluster.hierarchy import dendrogram 
from scipy.special import expit

from . utils import cobmine_samp_control, columnsforbox


my_colors  = ['cyan','gold','hotpink','peru','red','navy','purple','grey','tan','steelblue',
    'peachpuff','palegreen','orchid','darkred','black','orange','teal','firebrick','indigo','orchid','darkred','black',
    'orange','purple','grey','tan','brown','forestgreen',
    'cyan','gold','hotpink','peru','red','navy','purple','grey','tan','steelblue',
    'cyan','gold','hotpink','peru','red','navy','purple','grey','tan','steelblue',
    'peachpuff','palegreen','orchid','darkred','black','orange','teal','firebrick','indigo','orchid','darkred','black',
    'orange','purple','grey','tan','brown','forestgreen',
    'cyan','gold','hotpink','peru','red','navy','purple','grey','tan','steelblue','hotpink','peru','red','navy']


def get_graph():
    buffer = BytesIO()
    plt.savefig(buffer, format = 'svg')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()
    return graph


def get_auto_figsize(df_shape):
    if df_shape[0] < 50:
        return (df_shape[1], int(df_shape[0]/5)+2)
    else:          
        return (df_shape[1], int(df_shape[0]/10)+1)


def plot_pca(df_for_pca,sample_columns,control_columns,title, before,sname,cname):

    df = df_for_pca.transpose()

    # scaler = StandardScaler()
    # scaler.fit(df)
    # df_np_array = scaler.transform(df)

    df_np_array = df
    
    sample_list = list()
    i = 0

    for samples in sample_columns:
        for each_sample in samples:
            sample_list.append(sname[i])
        i += 1

    i = 0
    for control in control_columns:
        for each_control in control:
            sample_list.append(cname[i])
        i += 1



    pca = PCA(n_components=2)
    components = pca.fit_transform(df_np_array)

    percentagee=pca.explained_variance_ratio_
    per = [i * 100 for i in percentagee]
    per = ["{:.2f}".format(i) for i in per]

    pca_df = pd.DataFrame(components, columns = ['x','y'])
    pca_df['samples'] = pd.Series(sample_list)

    plt.switch_backend('AGG')


    if before:
        sns.scatterplot(data=pca_df,x=pca_df['x'],y=pca_df['y'],hue=pca_df['samples'], s = 80)
        # sns.scatterplot(data=pca_df,x=pca_df['x'],y=pca_df['y'],hue=pca_df['samples'], legend = False, s = 80)

    else:
        sns.scatterplot(data=pca_df,x=pca_df['x'],y=pca_df['y'],hue=pca_df['samples'],s=80)

    plt.legend(fontsize=6)

    plt.xlabel('PC1 ('+str(per[0])+'%)')

    plt.ylabel('PC2 ('+str(per[1])+'%)')

    plt.title(title)

    plt.axvline(x=0, linestyle='--', color='#7d7d7d', linewidth=1)
    plt.axhline(y=0, linestyle='--', color='#7d7d7d', linewidth=1)

    plt.tight_layout()

    pcaplot = get_graph()
    return pcaplot


def scale_number(unscaled, to_min, to_max, from_min, from_max):
    return (to_max-to_min)*(unscaled-from_min)/(from_max-from_min)+to_min

def scale_list(l, to_min, to_max):
    return [scale_number(i, to_min, to_max, min(l), max(l)) for i in l]


def plot_heatmap(df,both,fc_left,fc_right,lg2cut , z_score):

    auto_fig_size = get_auto_figsize(df.shape)

    cols = list()
    for x in df.columns:
        x = x.replace('LOG2 foldchange of','')
        x = x.replace('FOLDCHANGE_','')
        x = x.strip()
        cols.append(x)

    max_val = df.melt().value.max()
    min_val = df.melt().value.min()

    if both:
        df.fillna(0, inplace = True)
        cbar_kws={"ticks":[min_val,-lg2cut,lg2cut,max_val] }
        v = scale_list([ min_val, (min_val+(-lg2cut))/2,-lg2cut,0, lg2cut ,(lg2cut+max_val)/2, max_val] , 0,1)
        

    else:
        df.fillna(1, inplace = True)
        cbar_kws={"ticks":[min_val, fc_left,fc_right, max_val] }
        v = scale_list([ min_val, (min_val+(fc_left))/2,fc_left,1, fc_right ,(fc_right+max_val)/2, max_val] , 0,1)
    
    v = sorted(v)

    c = ["darkgreen","green","palegreen","black","lightcoral","red","darkred"]
    l = list(zip(v,c))

    cmap = LinearSegmentedColormap.from_list('rg',l, N=256)

    df.columns = cols
    plt.switch_backend('AGG')
    sns.set(font_scale=0.5)
    try:
        ax = sns.clustermap(df,yticklabels=True,cmap=cmap,cbar_kws=cbar_kws, ) 
    except:
        ax = sns.clustermap(df,yticklabels=True,) 

    heatmap = get_graph()
               
    plt.tight_layout()

    matplotlib.rc_file_defaults()
    return heatmap

def plot_heatmap_new(df,index,fc_left,fc_right,cutoff, plot_width , plot_height):

    df.set_index(index, inplace = True)
    df.fillna(0, inplace = True)

    if plot_width != None:
        auto_fig_size = (int(plot_width),int(plot_height))
    else:
        auto_fig_size = get_auto_figsize(df.shape)
    
    max_val = df.melt().value.max()
    min_val = df.melt().value.min()

    if cutoff == 'true':
        cbar_kws={"ticks":[min_val, fc_left,fc_right, max_val], 'label':'protein expression'}


        v = scale_list([ min_val, (min_val+(fc_left))/2,fc_left,fc_right ,(fc_right+max_val)/2, max_val] , 0,1)
        v = sorted(v)
        c = ["darkgreen","green","palegreen","black","lightcoral","red","darkred"]
        l = list(zip(v,c))
        cmap = LinearSegmentedColormap.from_list('rg',l, N=256)

        plt.switch_backend('AGG')
        sns.set(font_scale=0.5)
        cm = sns.clustermap(df,yticklabels=True,cmap=cmap, cbar_kws=cbar_kws, figsize=auto_fig_size)

    else:
        cm = sns.clustermap(df,yticklabels=True, figsize=auto_fig_size)

    cm.cax.set_aspect(0.05)


    plt.yticks(fontsize=5)
    heatmap = get_graph()
    
    plt.tight_layout()

    matplotlib.rc_file_defaults()
    return heatmap,auto_fig_size

def plot_volcano(df,lfc, pv, lg2cut, pvalue_cut_off, genes, title , both):
    
    difex_ptn = dict()

    plt.switch_backend('AGG')

    pv_thr = -(np.log10(pvalue_cut_off))
    lfc_thr = ( lg2cut,-lg2cut)

    color=("red", "grey", "green","black")

    df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] > pv_thr), 'color_add_axy'] = color[0]  # upregulated
    df.loc[(df[lfc] <= lfc_thr[1]) & (df[pv] > pv_thr), 'color_add_axy'] = color[2]  # downregulated
    df.loc[(df[lfc] > lfc_thr[1]) & (df[pv] > pv_thr) & (df[lfc] < lfc_thr[0]), 'color_add_axy'] = color[3]
    df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate

    filt = df['color_add_axy'] == color[0]
    to_lable_up =  df.loc[filt]
    
    to_lable_up.sort_values([lfc, pv], ascending=[False, False], inplace=True)
    to_lable_up = to_lable_up.head(5)


    filt1 = df['color_add_axy'] == color[2]
    to_lable_down =  df.loc[filt1]
    

    to_lable_down.sort_values([lfc, pv], ascending=[True, False], inplace=True)
    to_lable_down = to_lable_down.head(5)



    difex_ptn['Upregulated'] = df.loc[df['color_add_axy'] == color[0], genes].tolist()
    difex_ptn['Downregulated'] = df.loc[df['color_add_axy'] == color[2], genes].tolist()

    df['color_add_axy'].replace('nan', np.nan, inplace=True)                       #edit
    df['color_add_axy'].fillna('gray', inplace=True)                               #edit

    plt.scatter(x = df[lfc], y = df[pv], s = 26 , c = df['color_add_axy'], alpha=0.8)

    pv_max = df[pv].max()+1
    pv_min = df[pv].min()
    lfc_max = df[lfc].max()
    
    if lfc_max > 4:
        lfc_max = lfc_max + 1
    else:
        lfc_max = 4
    
    try:
        plt.xlim(-lfc_max,lfc_max)
        plt.ylim(pv_min,pv_max)
    
    except:
        plt.xlim(-4,4)
                          
    plt.axvline(x=lfc_thr[0], linestyle='--', color='#7d7d7d', linewidth=1,)
    plt.axvline(x=lfc_thr[1], linestyle='--', color='#7d7d7d', linewidth=1)
    plt.axhline(y=pv_thr, linestyle='--', color='#7d7d7d', linewidth=1)
    plt.ylabel('-log10 (p-value)', fontsize = 8)
    
    if both:
        plt.xlabel('log2 fold-change', fontsize = 8)
    else:
        plt.xlabel('fold-change', fontsize = 8)

    plt.title(title, fontsize = 10)


    for index, row in to_lable_up.iterrows():
        x = row[lfc]
        y = row[pv]
        gene_label = row[genes]

        plt.annotate(gene_label , xy = (x,y), xycoords='data', size= 8,
                                # xytext = (random.randint(0, 20),random.randint(0, 20)),
                                xytext = (30,22),

                                textcoords='offset pixels', ha='center', va='bottom',
                                # bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.3),
                                arrowprops=dict(arrowstyle='->',
                                                color='grey'))


    for index, row in to_lable_down.iterrows():
        x = row[lfc]
        y = row[pv]
        gene_label = row[genes]

        plt.annotate(gene_label , xy = (x,y), xycoords='data', size= 8,
                                # xytext = (random.randint(0, 20),random.randint(0, 20)),
                                xytext = (-30, 22),

                                textcoords='offset pixels', ha='center', va='bottom',
                                # bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.3),
                                   arrowprops=dict(arrowstyle='->',
                                                color='grey'))
    volcano_plt = get_graph()
    return volcano_plt, difex_ptn




def plot_volcano_new(df,index, pvalCutoff,fcCutoff, pv,lfc, lableProtiens, topInput , new_colours):
    genes = index
    plt.switch_backend('AGG')

    pv_thr = -(np.log10(pvalCutoff))
    #pv_thr = pvalCutoff
    lfc_thr = (fcCutoff,-fcCutoff)

    color=("red", "grey", "green","black")

    df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] > pv_thr), 'color_add_axy'] = color[0]  # upregulated
    df.loc[(df[lfc] <= lfc_thr[1]) & (df[pv] > pv_thr), 'color_add_axy'] = color[2]  # downregulated
    df.loc[(df[lfc] > lfc_thr[1]) & (df[pv] > pv_thr) & (df[lfc] < lfc_thr[0]), 'color_add_axy'] = color[3]
    df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate

    df['color_add_axy'].replace('nan', np.nan, inplace=True)                       #edit
    df['color_add_axy'].fillna('gray', inplace=True)                               #edit


    

    if new_colours:
        new_colours = new_colours.split(',')
        color=(new_colours[0], "grey", new_colours[1],"black")
    else:   
        color=("red", "grey", "green","black")

    df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] > pv_thr), 'color_add_axy'] = color[0]  # upregulated
    
    df.loc[(df[lfc] <= lfc_thr[1]) & (df[pv] > pv_thr), 'color_add_axy'] = color[2]  # downregulated

    df.loc[(df[lfc] > lfc_thr[1]) & (df[pv] > pv_thr) & (df[lfc] < lfc_thr[0]), 'color_add_axy'] = color[3]

    df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate

    filt = df['color_add_axy'] == color[0]
    to_lable_up =  df.loc[filt]

    filt1 = df['color_add_axy'] == color[2]
    to_lable_down =  df.loc[filt1]

    if lableProtiens == 'top':

        to_lable_up.sort_values([lfc, pv], ascending=[False, False], inplace=True)
        to_lable_up = to_lable_up.head(topInput)

        to_lable_down.sort_values([lfc, pv], ascending=[True, False], inplace=True)
        to_lable_down = to_lable_down.head(topInput)
    plt.scatter(x = df[lfc], y = df[pv], s = 32 , c = df['color_add_axy'], alpha=0.7)

    # plt.xlim(-4,4)

    plt.axvline(x=lfc_thr[0], linestyle='--', color='#7d7d7d', linewidth=1,)
    plt.axvline(x=lfc_thr[1], linestyle='--', color='#7d7d7d', linewidth=1)
    plt.axhline(y=pv_thr, linestyle='--', color='#7d7d7d', linewidth=1)
    plt.ylabel('-log10 (p-value)', fontsize = 8)
    plt.xlabel('log2 fold-change', fontsize = 8)

    for index, row in to_lable_up.iterrows():
        x = row[lfc]
        y = row[pv]
        gene_label = row[genes]

        plt.annotate(gene_label , xy = (x,y), xycoords='data', size= 8,
                                # xytext = (random.randint(0, 20),random.randint(0, 20)),
                                xytext = (40,20),

                                textcoords='offset pixels', ha='center', va='bottom',
                                # bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.3),
                                arrowprops=dict(arrowstyle='->',
                                                color='grey'))


    for index, row in to_lable_down.iterrows():
        x = row[lfc]
        y = row[pv]
        gene_label = row[genes]

        plt.annotate(gene_label , xy = (x,y), xycoords='data', size= 8,
                                # xytext = (random.randint(0, 20),random.randint(0, 20)),
                                xytext = (-40,20),

                                textcoords='offset pixels', ha='center', va='bottom',
                                # bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.3),
                                arrowprops=dict(arrowstyle='->',
                                                color='grey'))


    volcano_plt = get_graph()

    return volcano_plt



# https://stackoverflow.com/questions/9074996/matplotlib-how-to-annotate-point-on-a-scatter-automatically-placed-arrow

def box_plot_sns(df, title, samp_col, control_col, isbio):
    if isbio:
        new_cols, fordf =  cobmine_samp_control(samp_col, control_col)

        df = df[fordf]

    df.columns = columnsforbox(df.columns)

    df = np.log2(df)
    
    plt.switch_backend('AGG')

    flierprops = dict(marker='o', markerfacecolor='white', markersize=3,
                  linestyle='none', markeredgecolor='black')

    ax = sns.boxplot(data = df, notch=True, flierprops = flierprops, linewidth = 0.5, width = 0.5)

    colourdict = dict()

    if isbio:

        i = 0
        for sample in new_cols:
            for s in sample:
               colourdict[s] = i
            i += 1

    else:
        i = 0
        for sample in samp_col:
            for s in sample:
               colourdict[s] = i
            i += 1

        # i = 0
        for cntrl in control_col:
            for c in cntrl:
               colourdict[c] = i
            i += 1


    for patch, c_dict in zip(ax.patches,colourdict):
        patch.set_facecolor(my_colors[colourdict[c_dict]])

    plt.xticks(fontsize=6, rotation=90)
    plt.ylabel('log2 of Abundances')
    plt.title(title)
    plt.tight_layout()
    box_plot = get_graph()
    return box_plot



def plot_bar_grap(df):

    # uplist_bar = []
    # dowlist_bar = []
    # samples_bar = list(df.columns)
    # samples_bar = [i.replace('LOG2FC-Expression_','') for i in samples_bar]
    # samples_bar = [i.replace('FC-Expression_','') for i in samples_bar]

    # for column in df:
        
    #     uplist_bar.append(len(df.loc[df[column] == 'Upregulated']))

    #     dowlist_bar.append(len(df.loc[df[column] == 'Downregulated']))


    # uplist_bar = json.dumps(uplist_bar)
    # dowlist_bar = json.dumps(dowlist_bar)
    # return uplist_bar , dowlist_bar , samples_bar


    plt.switch_backend('AGG')

    cols = list()
    for x in df.columns:
        if 'LOG2FC-Expression_' in x:
            x = x.replace('LOG2FC-Expression_','')
        if 'FC-Expression_' in x:
            x = x.replace('FC-Expression_','')
        x = x.strip()
        cols.append(x)
    df.columns = cols

    updict = dict()
    downdict = dict()


    for column in df:
        updict[column] = len(df.loc[df[column] == 'Upregulated'])
        downdict[column] = len(df.loc[df[column] == 'Downregulated'])

    upmax = max(updict.values())
    downmax = max(downdict.values())

    fig, ax = plt.subplots()
    plt.xticks(rotation=90, fontsize = 6)

    ax.bar(updict.keys(),updict.values(), color = 'red')
    ax.bar(downdict.keys(),[x*-1 for x in downdict.values()], color = 'green')

    ax.set_ylim(-downmax-5, upmax+5)

    upareg = list(updict.values())
    downreg = [x*-1 for x in downdict.values()]


    for i in range(len(upareg)):
        plt.text(i,upareg[i], upareg[i], ha = 'center')


    for j in range(len(downreg)):
        plt.text(j,downreg[j]-2, abs(downreg[j]), ha = 'center', va = 'bottom')

    ticks =  ax.get_yticks()
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])
    plt.axhline(y= 0 , linestyle='-', color='#7d7d7d', linewidth=1)

    custom_lines = [Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='green', lw=4)]
    plt.ylabel("Number of proteins")
    ax.legend(custom_lines, ['Upregulated', 'Downregulated'])

    plt.tight_layout()
    bar_graph = get_graph()
    return bar_graph



#all code below are for circular bargraph
#---------------------------------------------------------------------
def getgroupsize(golist):
    sep_go = [list(grp) for k, grp in groupby(golist)]
    counts = [len(x) for x in sep_go]
    return counts

def getfunctions(golist):
    sep_go = [list(grp) for k, grp in groupby(golist)]
    functions = [x[0] for x in sep_go]
    return functions
# The following is a helper function that given the angle at which the bar is positioned and the
# offset used in the barchart, determines the rotation and alignment of the labels.

def get_label_rotation(angle, offset):
    # Rotation must be specified in degrees :(
    rotation = np.rad2deg(angle + offset)
    if angle <= np.pi:
        alignment = "right"
        rotation = rotation + 180
    else:
        alignment = "left"
    return rotation, alignment


def add_labels(angles, values, labels, offset, ax):
    # This is the space between the end of the bar and the label
    padding = 4

    # Iterate over angles, values, and labels, to add all of them.
    for angle, value, label, in zip(angles, values, labels):
        angle = angle

        # Obtain text rotation and alignment
        rotation, alignment = get_label_rotation(angle, offset)

        # And finally add the text
        ax.text(
            x=angle,
            y=value + padding,
            s=label,
            ha=alignment,
            va="center",
            rotation=rotation,
            rotation_mode="anchor",
            size = 8
        )

def getcircbar(df):

    plt.switch_backend('AGG')
    df["strvalue"] = df["value"].astype(str)
    df["strvalue"] = df["strvalue"].apply(lambda x: '('+x+')')
    df["name"]  = df[['strvalue', 'name']].agg('-'.join, axis=1)


    OFFSET = np.pi / 2

    VALUES = np.log2(df["value"].values)
    LABELS = df["name"].values
    GROUP = df["group"].values

    PAD = 3
    ANGLES_N = len(VALUES) + PAD * len(np.unique(GROUP))
    ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
    WIDTH = (2 * np.pi) / len(ANGLES)

    offset = 0
    IDXS = []

    GROUPS_SIZE = getgroupsize(df['group'].tolist())

    for size in GROUPS_SIZE:
        IDXS += list(range(offset + PAD, offset + size + PAD))
        offset += size + PAD



    fig, ax = plt.subplots(figsize=(20, 10), subplot_kw={"projection": "polar"})
    # fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

    ax.set_ylim(-100, max(VALUES))
    ax.set_theta_offset(OFFSET)
    ax.set_frame_on(False)
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    GROUPS_SIZE = getgroupsize(df['group'].tolist())

    COLORS = [f"C{i}" for i, size in enumerate(GROUPS_SIZE) for _ in range(size)]

    ax.bar(
        ANGLES[IDXS], VALUES, width=WIDTH, color=COLORS,
        edgecolor="white", linewidth=2
    )

    add_labels(ANGLES[IDXS], VALUES, LABELS, OFFSET, ax)

    # Extra customization below here --------------------

    # This iterates over the sizes of the groups adding reference
    # lines and annotations.
# Set the coordinates limits


    offset = 0
    funcs = getfunctions(df['group'].tolist())

    for group, size in zip(funcs, GROUPS_SIZE):
        # Add line below bars
        x1 = np.linspace(ANGLES[offset + PAD], ANGLES[offset + size + PAD - 1], num=50)
        ax.plot(x1, [-5] * 50, color="#333333")

        # Add text to indicate group
        ax.text(
            np.mean(x1), -20, group, color="#333333", fontsize=8,
            fontweight="bold", ha="center", va="center"
        )

        # Add reference lines at 20, 40, 60, and 80
        x2 = np.linspace(ANGLES[offset], ANGLES[offset + PAD - 1], num=50)
        ax.plot(x2, [20] * 50, color="#bebebe", lw=0.8)
        ax.plot(x2, [40] * 50, color="#bebebe", lw=0.8)
        ax.plot(x2, [60] * 50, color="#bebebe", lw=0.8)
        ax.plot(x2, [80] * 50, color="#bebebe", lw=0.8)

        offset += size + PAD

    plt.tight_layout()
    get_circbar = get_graph()
    return get_circbar

#------------------------------------------------ cirbar done

def getElbowPlot(df,key):
    
    n_c = 15

    df.set_index(key, inplace = True)

    df = df.select_dtypes(exclude=['object'])

    cols = list()
    for x in df.columns:
        x = x.replace('LOG2 foldchange of','')
        x = x.replace('FOLDCHANGE_','')
        x = x.strip()
        cols.append(x)

    df.fillna(0, inplace = True)
    df.columns = cols
    plt.switch_backend('AGG')

    if df.shape[0] <= 15:
         n_c = df.shape[0]
    wcss=[]
    for i in range(1,n_c):
        kmeans=KMeans(n_clusters=i,init="k-means++",random_state=42)
        kmeans.fit(df)
        wcss.append(kmeans.inertia_)
    clusters=np.arange(1, n_c) 

    plt.plot(clusters,wcss,marker ='x', markeredgecolor = 'blue')

    plt.tight_layout()
    elbowmap = get_graph()
    return elbowmap


def getkClusterMap(df,key,n_clusters,fc_left,fc_right,lg2cut,both):

    auto_fig_size = get_auto_figsize(df.shape)
    plt.switch_backend('AGG')

    df.set_index(key, inplace = True)
    
    df = df.select_dtypes(exclude=['object'])

    cols = list()

    for x in df.columns:
        x = x.replace('LOG2 foldchange of','')
        x = x.replace('FOLDCHANGE_','')
        x = x.strip()
        cols.append(x)
    df.columns = cols

    max_val = df.melt().value.max()
    min_val = df.melt().value.min()

    # max_val = float("{:.2f}".format(max_val))
    if both:
        df.fillna(0, inplace = True)
        cbar_kws={"ticks":[min_val,-lg2cut,lg2cut,max_val] , "orientation": "horizontal"}
        v = scale_list([ min_val, (min_val+(-lg2cut))/2,-lg2cut,0, lg2cut ,(lg2cut+max_val)/2, max_val] , 0,1)

    else:
        df.fillna(1, inplace = True)
        cbar_kws={"ticks":[min_val, fc_left,fc_right, max_val], "orientation": "horizontal" }
        v = scale_list([ min_val, (min_val+(fc_left))/2,fc_left,1, fc_right ,(fc_right+max_val)/2, max_val] , 0,1)

    c = ["darkgreen","green","palegreen","black","lightcoral","red","darkred"]
    l = list(zip(v,c))

    cmap = LinearSegmentedColormap.from_list('rg',l, N=256)

    kmeans=KMeans(n_clusters=n_clusters,max_iter=1000,random_state=42)
    kmeans.fit(df)

    df['clusters']=kmeans.labels_

    df=df.sort_values('clusters')

    data=df.loc[:, df.columns != 'clusters']
    
    sns.set(font_scale=0.5)

    fig, axes = plt.subplots(figsize=(auto_fig_size))

    ax = sns.heatmap(data,yticklabels=True,cmap=cmap,cbar_kws=cbar_kws)

    clusters=list(df['clusters'])

    hline_list=[0]
    points=0
    for i in range(0,n_clusters):
        counts=clusters.count(i)
        points=points+counts
        hline_list.append(points)

    # ax.hlines(hline_list,*ax.get_xlim(),colors='white', linewidth = 3)

    cmap= matplotlib.cm.viridis

    bounds = hline_list

    bounds=list(reversed(bounds))
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    ax1=fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        ax=ax,  # Specify the axes to use for the colorbar
        ticks=bounds,
        spacing='proportional',
        orientation='vertical',
        #label='Discrete intervals, some other units',
        )


    plt.tight_layout()
    kmap = get_graph()
    
    matplotlib.rc_file_defaults()

    return kmap , df


def get_upset_plot(df):

    plt.switch_backend('AGG')
    columns={col: pd.unique(df[col].dropna().to_list()) for col in df}  #remove nan from list that comes since the columns have diff number of values
    data=from_contents(columns)
    UPSET_PLOT(data)
    
    plt.tight_layout()
    plott = get_graph()
    return plott

def get_bubble_plot(df,xAxis, yAxis,sizeCol, category):

    plt.switch_backend('AGG')
    sns.scatterplot(data=df, x=df[xAxis], y=df[yAxis], size=df[sizeCol], legend=False, hue=df[category], sizes=(50, 1000),linewidth = 0)
    plt.title("Bubble plot")
    plt.tight_layout()
    plott = get_graph()
    return plott

def get_box_plot(df,fields, colors):
    plt.switch_backend('AGG')

    df1=df.filter(fields)
    
    if colors:
        colors = colors.split(',')
        sns.boxplot(x="variable", y="value", data=pd.melt(df1), palette = colors)
    else:
        sns.boxplot(x="variable", y="value", data=pd.melt(df1))

    plt.xticks(rotation=90, fontsize=6)
    plt.tight_layout()

    boxplot = get_graph()
    return boxplot


def get_density_plot(df,fields):

    plt.switch_backend('AGG')
    newdf=df[fields]
    sns.set(style="darkgrid")
    #for f in fields:
    sns.kdeplot(data=newdf, fill=True, common_norm=False,alpha=0.2,legend=True)    
    plt.tight_layout()
    density_plott = get_graph()
    return density_plott

def get_histogram(df,fields):
    plt.switch_backend('AGG')
    sns.set(style="darkgrid")

    df1=df.filter(fields)
    ax=sns.histplot(x="value", hue="variable", data=pd.melt(df1), kde=True)
    # sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    
    plt.tight_layout()
    plott = get_graph()
    return plott

    
def get_pca_plot(df1):

    plt.switch_backend('AGG')

    df = df1.iloc[: , 1:]

    df = df.transpose()

    pca = PCA(n_components=2)
    components = pca.fit_transform(df)
    
    pca_df = pd.DataFrame(components, columns = ['x','y'])
    value=df.index.values


    sns.scatterplot(data=pca_df,x=pca_df['x'],y=pca_df['y'],hue= value ,s=150)
        
    plt.legend(fontsize=6)
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    plt.tight_layout()
    pcaplot = get_graph()
    return pcaplot


def plot_kmeans_new(df,index, fc_left, fc_right,n_clusters):
    auto_fig_size = get_auto_figsize(df.shape)

    plt.switch_backend('AGG')

    df.set_index(index, inplace = True)

    max_val = df.melt().value.max()
    min_val = df.melt().value.min()

    # max_val = float("{:.2f}".format(max_val))
    df.fillna(0, inplace = True)
   

    
    cbar_kws={"ticks":[min_val, fc_left,fc_right, max_val] ,"orientation": "horizontal"}
    v = scale_list([ min_val, (min_val+(fc_left))/2,fc_left,fc_right ,(fc_right+max_val)/2, max_val] , 0,1)
    v = sorted(v)
    c = ["darkgreen","green","palegreen","black","lightcoral","red","darkred"]
    l = list(zip(v,c))

    cmap = LinearSegmentedColormap.from_list('rg',l, N=256)

    kmeans=KMeans(n_clusters=n_clusters,max_iter=1000,random_state=42)
    kmeans.fit(df)

    df['clusters']=kmeans.labels_

    df=df.sort_values('clusters')
    
    data=df.loc[:, df.columns != 'clusters']

    sns.set(font_scale=0.5)

    fig, axes = plt.subplots(figsize = auto_fig_size)

    ax = sns.heatmap(data,yticklabels=True,cmap=cmap,cbar_kws=cbar_kws)

    clusters=list(df['clusters'])

    hline_list=[]
    hline_list.append(0)
    points=0
    for i in range(0,n_clusters):
        counts=clusters.count(i)
        points=points+counts
        hline_list.append(points)

    cmap= matplotlib.cm.viridis
    bounds = hline_list
    bounds=list(reversed(bounds))
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    ax1=fig.colorbar(
        matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm),
        ax=ax,  # Specify the axes to use for the colorbar
        ticks=bounds,
        spacing='proportional',
        orientation='vertical',
        )
    
    plt.tight_layout()
    kmap = get_graph()
    return kmap


# def get_raincloud_plot(input_df,fieldsRain):
#     df =  pd.DataFrame()
   
#     for col in fieldsRain:
#         values = input_df[col].tolist()
#         group = [col for i in range(0,len(values))]
#         dff = pd.DataFrame(values,group)
#         df = pd.concat([df,dff])

#     df.reset_index(inplace = True)
#     df.columns = ['group','score']
#     plt.switch_backend('AGG')
#     f, ax = plt.subplots(figsize=(7, 5))
#     dy = 'group'
#     dx= 'score'
#     ort="h"
#     pal = sns.color_palette(n_colors=1)

#     pal = "Set2"
#     f, ax = plt.subplots(figsize=(7, 5))
#     ax=pt.half_violinplot( x = dx, y = dy, data = df, palette = pal, bw = .2, cut = 0.,
#     scale = "area", width = .6, inner = None, orient = ort)
#     ax=sns.stripplot( x = dx, y = dy, data = df, palette = pal, edgecolor = "white",
#     size = 3, jitter = 1, zorder = 0, orient = ort)
#     ax=sns.boxplot( x = dx, y = dy, data = df, color = "black", width = .15, zorder = 10,\
#     showcaps = True, boxprops = {'facecolor':'none', "zorder":10},\
#     showfliers=True, whiskerprops = {'linewidth':2, "zorder":10},\
#     saturation = 1, orient = ort)

#     rainplot = get_graph()
#     return rainplot
def get_raincloud_plot(input_df, fieldsRain):

    df = pd.DataFrame()

    for col in fieldsRain:
        values = input_df[col].tolist()
        group = [col for i in range(0, len(values))]
        dff = pd.DataFrame({'group': group, 'score': values})
        df = pd.concat([df, dff])

    plt.switch_backend('AGG')
    f, ax = plt.subplots(figsize=(7, 5))

    dy = 'group'
    dx = 'score'
    ort = "h"
    pal = sns.color_palette("Set2", n_colors=len(fieldsRain))

    sns.violinplot(x=dx, y=dy, data=df, palette=pal, bw=.2, cut=0.,
                   scale="area", width=.6, inner=None, orient=ort, ax=ax)
    
    sns.stripplot(x=dx, y=dy, data=df, palette=pal, edgecolor="white",
                  size=3, jitter=1, zorder=0, orient=ort, ax=ax)
    
    sns.boxplot(x=dx, y=dy, data=df, color="black", width=.15, zorder=10,
                showcaps=True, boxprops={'facecolor': 'none', "zorder": 10},
                showfliers=True, whiskerprops={'linewidth': 2, "zorder": 10},
                saturation=1, orient=ort, ax=ax)

    rainplot = get_graph()
    return rainplot


def get_violin_plot(df,fields, colors):
    plt.switch_backend('AGG')

    df1=df.filter(fields)

    if colors:
        colors = colors.split(',')
        sns.violinplot(x="variable", y="value", data=pd.melt(df1), palette = colors )
    else:
        sns.violinplot(x="variable", y="value", data=pd.melt(df1))

    plt.xticks(rotation=90, fontsize=6)
    plt.tight_layout()

    violinplot = get_graph()
    return violinplot


def get_scurve_plot(df, axis):
    value=axis
    plt.switch_backend('AGG')

    x_data = df[axis]
    sigmoid_x=expit(x_data)
    sns.set(style="darkgrid")
    ax=sns.scatterplot(x=x_data,y=sigmoid_x,linewidth = 0.1)
    ax.set(xlabel=value,ylabel="Sigmoid (" + value + ")")

    plt.tight_layout()
    scurve = get_graph()
    return scurve

def  get_venn_plot(df,sets):
    
    plt.switch_backend('AGG')
    allsets=dict()
    for s in sets:
        #bla={s:set(df[s].dropna())}
        allsets.update({s:set(df[s].dropna())})
    venn(allsets)

    ven = get_graph()
    return ven

def get_ma_plot(df,fc_cutoff,sample,control,log2fc):
    plt.switch_backend('AGG')
    cutoff=np.log2(fc_cutoff)
    # log value calcultions
    log2_sample=np.log2(df[sample])
    log2_control=np.log2(df[control])  
    
    x_data=(log2_sample+log2_control)/2 #average of log of normalized counts/abundance of control and sample i.e A
    x_data=x_data.dropna()
    
    y_data=df[log2fc]  #log2fc i.e M i.e (log(sample count/control count))
    y_data=y_data.dropna()
    
    # colur
    max_val = y_data.max()
    
    min_val = y_data.min()
    if (-cutoff) <= min_val:
        min_val=(-cutoff)-1
    
    colors = pd.cut(
    y_data.dropna(),
    bins=[min_val, -cutoff,+cutoff,(max_val+1)],
    labels=["red", "grey", "green"],
    right=False
    )
    #plot

    #fig, ax = plt.subplots()
    #ax.scatter(x=x_data, y=y_data, c=colors)
    sns.scatterplot(x=x_data,y=y_data,linewidth = 0.05,c=colors)
    plt.axhline(y=0, linestyle='--', color='#7d7d7d', linewidth=1.5)
    plt.xlabel("Log2 mean expression")
    plt.ylabel("Log2 fold-change")
    plt.tight_layout()
    maplot = get_graph()
    return maplot
