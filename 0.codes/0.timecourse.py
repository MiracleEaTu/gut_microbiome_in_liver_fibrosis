import pandas as pd
import scipy.stats
import numpy as np


def test():
    df_in = pd.read_excel('../1.rawdata/OTU_table_20210611.xlsx', index_col=0)
    correlation_list = find_correlation_by_pearson_test_with_PCC(df_in)
    get_clusters(correlation_list)
    file_out = open('timecourse_clusters.txt', 'w')
    file_out.write('node1\tnode2\n')
    for relation_list in correlation_list:
        if relation_list[0] < relation_list[1]:
            file_out.write(str(relation_list[0])+'\t'+str(relation_list[1])+'\n')
    file_out.close()
    print_clusters(correlation_list)
    pass


def find_correlation_by_pearson_test_with_PCC(matrix_pandas_dataframe):
    index_list = list(matrix_pandas_dataframe.index.values)
    correlation_list = []
    data_num = len(matrix_pandas_dataframe.iloc[:, 0])
    print(data_num)
    # print(scipy.stats.pearsonr(matrix_pandas_dataframe.iloc[69, :], matrix_pandas_dataframe.iloc[60, :]))
    for i in range(data_num):
        for j in range(data_num):
            pearsonr_result = scipy.stats.pearsonr(matrix_pandas_dataframe.iloc[i, :], matrix_pandas_dataframe.iloc[j, :])
            if pearsonr_result[0] >= 0.6 and i != j:
                correlation_list.append([index_list[i], index_list[j], pearsonr_result[0]])
    # print(correlation_list) 
    return correlation_list


def get_main_clusters(correlation_list):
    max_num = 0
    cluster_main_component = ''
    target_num = 0
    for i in range(len(correlation_list)-1):
        if correlation_list[i][0] == correlation_list[i+1][0]:
            target_num += 1
            if target_num >= max_num:
                cluster_main_component = correlation_list[i][0]
                max_num = target_num
            else:
                pass
        else:
            target_num = 0
    # print(max_num+1)
    # print(cluster_main_component)
    main_cluster = [cluster_main_component]
    for connections in correlation_list:
        if connections[0] == cluster_main_component:
            main_cluster.append(connections[1])
    # print(main_cluster)
    return main_cluster


def get_clusters(correlation_list):
    all_clusters = []
    cluster = get_main_clusters(correlation_list)
    all_clusters.append(cluster)
    while len(cluster) >= 2:
        correlation_list = [connections for connections in correlation_list if connections[0] not in cluster and connections[1] not in cluster]
        cluster = get_main_clusters(correlation_list)
        if cluster != ['']:
            all_clusters.append(cluster)
    return all_clusters


def print_clusters(correlation_list):
    # to print clusters, delete # in this function
    all_clusters = []
    fOUT = open('ctimecourse_clusters_group.txt', 'w')
    fOUT.write('nodeID\tClusters\n')
    cluster_names = ['Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'ClusterOthers']
    n = 0
    cluster = get_main_clusters(correlation_list)
    all_clusters.append(cluster)
    for elements in cluster:
        fOUT.write(str(elements)+'\t'+cluster_names[n]+'\n')
    n += 1
    while len(cluster) >= 2:
        correlation_list = [connections for connections in correlation_list if connections[0] not in cluster and connections[1] not in cluster]
        cluster = get_main_clusters(correlation_list)
        if cluster != ['']:
            all_clusters.append(cluster)
            if n == 11:
                n = 10
            for elements in cluster:
                fOUT.write(str(elements)+'\t'+cluster_names[n]+'\n')
            n += 1
    # print(all_clusters)
    fOUT.close()
    return all_clusters


if __name__ == '__main__':
    test()