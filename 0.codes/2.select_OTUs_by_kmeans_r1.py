from operator import index
from matplotlib import projections
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from matplotlib.ticker import MultipleLocator
from collections import OrderedDict
__author__ = 'Qijie Guan'

def get_plot_index():
    # 取CCl4喂食樟芝菌粉importance前15与Control前15的并集
    df_JFD2CCl4 = pd.read_excel('../2.output_files/RandomForestTree_featureImportance_JFD2CCl4.xlsx', index_col=0)
    df_JFH2CCl4 = pd.read_excel('../2.output_files/RandomForestTree_featureImportance_JFH2CCl4.xlsx', index_col=0)
    df_CTL2CCl4 = pd.read_excel('../2.output_files/RandomForestTree_featureImportance_CTL2CCl4.xlsx', index_col=0)
    index_using = df_JFD2CCl4.index.values
    # k-means cluster对并集进行聚类，取8类，输出为图像
    X = []
    for index_name in index_using:
        value_JFD2CCl4 = df_JFD2CCl4.loc[index_name, 'Importance']
        value_JFH2CCl4 = df_JFH2CCl4.loc[index_name, 'Importance']
        value_CTL2CCl4 = df_CTL2CCl4.loc[index_name, 'Importance']
        X.append([value_JFD2CCl4, value_JFH2CCl4, value_CTL2CCl4])
    # print(X)
    # print(np.array(X))
    X = np.array(X)
    s1 = StandardScaler()
    X = s1.fit_transform(X)
    kmeans = KMeans(n_clusters=8, random_state=0).fit(X)
    km_labels = kmeans.labels_
    fig1 = plt.figure()
    # ax = fig1.gca(projection='3d')
    ax = Axes3D(fig1)
    plot_x = X[:,0]
    plot_y = X[:,1]
    plot_z = X[:,2]
    # Data output
    df_feature_clusters = pd.DataFrame({'Features':index_using, 'Importance_JFD2CCl4':plot_x,\
        'Importance_JFDH2CCl4':plot_y, 'Importance_CTL2CCl4':plot_z, 'Cluster':km_labels})
    df_feature_clusters.to_excel('../3.figures_files/Feature_clusters.xlsx', index=False)
    colors_using = ['#000000', '#be0027', '#cf8d2e', '#e4e932',\
    '#2c9f45', '#371777', '#037ef3', '#00a4e4']
    for i in range(len(km_labels)):
        ax.scatter(plot_x[i], plot_y[i], plot_z[i], c=colors_using[km_labels[i]], label='Cluster %i'%km_labels[i])
    ax.set_zlabel('JFD vs CCl4')
    ax.set_ylabel('JFH vs CCl4')
    ax.set_xlabel('CTL vs CCl4')
    # ax.legend()
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    # ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
    plt.savefig('../3.figures_files/Feature_clusters_20210608.pdf')
    # plt.show()
    plt.close()
    return 0


def Core_cluster_selection(cluster_select=[3,4,7]):
    df_feature_clusters = pd.read_excel('../3.figures_files/Feature_clusters.xlsx', index_col=0)
    index_important = []
    for index_names in df_feature_clusters.index.values:
        for cluster_id in cluster_select:
            if df_feature_clusters.loc[index_names,'Cluster'] == cluster_id:
                index_important.append(index_names)
                break
    df_important = df_feature_clusters.loc[index_important]
    df_important.to_excel('../3.figures_files/Feature_clusters_important_20210612.xlsx', index=True)


def plot_core_clusters():
    df_important = pd.read_excel('../3.figures_files/Feature_clusters_important_20210612.xlsx', index_col=0)
    index_important = df_important.index.values
    # get dictionary, keys = OTU number, value = taxonomy
    df_all_data = pd.read_csv('../1.rawdata/OTU_table.from_biom_w_taxonomy_filtered.csv', index_col=0)
    d_OTU2tax = {}
    for index_name in df_all_data.index.values:
        d_OTU2tax[index_name] = df_all_data.loc[index_name, 'taxonomy']
    # print(d_OTU2tax)
    del(df_all_data)
    taxonomy_important = []
    for index_values in index_important:
        taxonomy_important.append(d_OTU2tax[str(index_values)])
    x_labels = []
    for treatment in ['CTL', 'CCl4', 'JFH-CCl4', 'JFD-CCl4']:
        for days in ['1', '14', '28']:
            x_labels.append('%s_day%s'%(treatment, days))
    df_OTU_table = pd.read_csv('../2.output_files/summarized_OTUs_20210612.csv', index_col=0)
    Z = []
    for index_name in index_important:
        for treatment in ['CTL', 'CCl4', 'JFH-CCl4', 'JFD-CCl4']:
            for days in ['1', '14', '28']:
                Z.append(df_OTU_table.loc[index_name, 'D%s-%s_values'%(days, treatment)])
                X = np.arange(12)
    X = np.arange(12)
    Y = np.arange(len(index_important))
    Z = np.array(Z)
    # 网格化坐标
    xx, yy = np.meshgrid(X, Y)
    # 矩阵扁平化
    X, Y = xx.ravel(), yy.ravel()
    # 
    height = np.zeros_like(Z)
    width = depth = 0.6
    fig = plt.figure(figsize=(30,30))
    ax = fig.gca(projection='3d')
    ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([0.7, 2, 1, 1]))
    df_FMT_top50 = pd.read_excel('../1.rawdata/FMT_OTUs_20210612.xlsx', index_col=0)
    indexes_FMT_top50 = list(df_FMT_top50.index.values)[:20]
    taxonomy_FMT_top50 = []
    for index_name in indexes_FMT_top50:
        taxonomy_FMT_top50.append(d_OTU2tax[str(index_name)])
    colors_for_plots = []
    print(taxonomy_important)
    print(taxonomy_FMT_top50)
    for tax_names in taxonomy_important:
        if tax_names in taxonomy_FMT_top50:
            colors_for_plots = colors_for_plots + ['#F5474890']*12
            print('match')
        else:
            colors_for_plots = colors_for_plots + ['#f5e6ca90']*12
            print('NaN')
    # ax.bar3d(X, Y-0.6, height, width, depth, Z,  color='#F5474890', shade=True)
    ax.bar3d(X, Y-0.6, height, width, depth, Z,  color=colors_for_plots, shade=True)
    y_axis = []
    for treatment in ['CTL', 'CCl4', 'JFH-CCl4', 'JFD-CCl4']:
        for days in ['1', '14', '28']:
            y_axis.append('%s_day%s'%(treatment, days))
    ax.set_xticklabels(y_axis, minor=True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xticklabels(x_labels)
    # ax.set_yticklabels(list(index_important))
    x_major_locator=MultipleLocator(1)
    ax.xaxis.set_major_locator(x_major_locator) 
    y_major_locator=MultipleLocator(1)
    ax.yaxis.set_major_locator(y_major_locator)
    # ax.set_adjustable('box')
    # plt.tight_layout()
    plt.savefig('../3.figures_files/Core_cluster_in_samples_w_FMT_coloring_20210614.pdf')
    # plt.show()

    ## print No. 2
    Z2 = []
    for index_name in indexes_FMT_top50:
        Z2.append(df_FMT_top50.loc[index_name, '626AC1'])
        for treatment in ['CTL', 'CCl4', 'JFH-CCl4', 'JFD-CCl4']:
            for days in ['1', '14', '28']:
                Z2.append(df_OTU_table.loc[str(index_name), 'D%s-%s_values'%(days, treatment)])
    X2 = np.arange(13)
    Y2 = np.arange(len(indexes_FMT_top50))
    Z2 = np.array(Z2)
    # 网格化坐标
    xx2, yy2 = np.meshgrid(X2, Y2)
    # 矩阵扁平化
    X2, Y2 = xx2.ravel(), yy2.ravel()
    # 
    height = np.zeros_like(Z2)
    width = depth = 0.6
    fig2 = plt.figure(figsize=(30,30))
    ax2 = fig2.gca(projection='3d')
    ax2.get_proj = lambda: np.dot(Axes3D.get_proj(ax2), np.diag([0.7, 2, 1, 1]))
    
    colors_for_plots2 = []
    
    for tax_names in taxonomy_FMT_top50:
        if tax_names in taxonomy_important:
            colors_for_plots2 = colors_for_plots2 + ['#F5474890']*13
            print('match')
        else:
            colors_for_plots2 = colors_for_plots2 + ['#f5e6ca90']*13
            print('NaN')
    # ax.bar3d(X, Y-0.6, height, width, depth, Z,  color='#F5474890', shade=True)
    ax2.bar3d(X2, Y2-0.6, height, width, depth, Z2,  color=colors_for_plots2, shade=True)
    y_axis = []
    y_axis.append('FMT')
    for treatment in ['CTL', 'CCl4', 'JFH-CCl4', 'JFD-CCl4']:
        for days in ['1', '14', '28']:
            y_axis.append('%s_day%s'%(treatment, days))
    
    ax2.set_xticklabels(y_axis, minor=True)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_xticklabels(x_labels)
    # ax.set_yticklabels(list(index_important))
    x_major_locator=MultipleLocator(1)
    ax2.xaxis.set_major_locator(x_major_locator) 
    y_major_locator=MultipleLocator(1)
    ax2.yaxis.set_major_locator(y_major_locator)
    # ax.set_adjustable('box')
    # plt.tight_layout()
    plt.savefig('../3.new_outputs/Top_20_FMT_OTUs_in_samples_w_core_cluster_coloring_20210614.pdf')
    plt.show()


if __name__ == '__main__':
    get_plot_index()
    Core_cluster_selection(cluster_select=[3,4,7])
    plot_core_clusters()