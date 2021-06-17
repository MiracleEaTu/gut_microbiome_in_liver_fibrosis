from matplotlib.cbook import index_of
import pandas as pd 
import numpy as np 

def get_keys(d, value):
    return [k for k,v in d.items() if v == value]


def data_filter_by_std(df_in_path=''):
    # 将sample-id对应的样本信息对应上
    metadata = pd.read_table('../1.data/metadata.txt', index_col=0)
    dict_sampleid2sampletype = {}
    for index_id in metadata.index.values:
        dict_sampleid2sampletype[index_id] = metadata.loc[index_id, 'SampleType']
    print(dict_sampleid2sampletype)
    # 开始对原始数据进行处理
    df = pd.read_excel(df_in_path, index_col=0)
    df_new = pd.DataFrame(None, index=df.index.values)
    sample_types_all = list(set(dict_sampleid2sampletype.values()))
    for sample_types in sample_types_all:
        columns_using = [column_names for column_names in df.columns.values if sample_types in column_names]
        df4keys = df[columns_using].T
        mean_values = np.mean(df4keys)
        std_values = np.std(df4keys)
        df_new['%s_values'%sample_types] = mean_values
        df_new['%s_std'%sample_types] = std_values
    df_new.to_csv('../1.data/summarized_OTUs_20210612.csv', index=True)
    return df_new


def data_rename_rows(df_in_path=''):
    # 将sample-id对应的样本信息对应上
    metadata = pd.read_table('../1.data/metadata.txt', index_col=0)
    dict_sampleid2sampletype = {}
    for index_id in metadata.index.values:
        dict_sampleid2sampletype[index_id] = metadata.loc[index_id, 'SampleType']
    print(dict_sampleid2sampletype)
    # 开始对原始数据进行处理
    df = pd.read_csv(df_in_path, index_col=0)
    print(df.columns.values)
    column_names_replace = []
    columns_selected = []
    # 最后一列为物种注释信息
    for column_name in df.columns.values[:-1]:
        if int(column_name) in dict_sampleid2sampletype.keys():
            print(column_name)
            columns_selected.append(column_name)
            sample_type = dict_sampleid2sampletype[int(column_name)]
            column_names_replace.append(sample_type)
    # df_new.columns = column_names_replace
    df_new = df[columns_selected]
    df_new.columns = column_names_replace
    df_new.to_excel('../1.data/OTU_table_20210611.xlsx')
    return df_new


def filter_data_based_on_OTU_values(df_in, df_out_name=''):
    # 每个样本OTU数据存放在列
    # 每个OTU在各个样本中的数据存放在行
    # df_out = pd.DataFrame(None, columns=df_in.columns.values)
    index_using = []
    for index_name in df_in.index.values:
        # 最后一列为taxonomy
        sum_OTU_in_samples = np.sum(df_in.loc[index_name][:len(df_in.loc[index_name])-1])
        if sum_OTU_in_samples >= 100:
            index_using.append(index_name)
    df_out = df_in.loc[index_using]
    df_out.to_csv('../2.output_files/%s.csv'%df_out_name)
    return 0


def test():
    df_in = pd.read_table('../1.data/table.from_biom_w_taxonomy_bak.txt', index_col=0)
    filter_data_based_on_OTU_values(df_in=df_in, df_out_name='table.from_biom_w_taxonomy_filtered')
    data_rename_rows(df_in_path='../1.data/OTU_table.from_biom_w_taxonomy_filtered.csv')
    df_filtered = data_filter_by_std(df_in_path='../1.data/OTU_table_20210611.xlsx')


if __name__ == '__main__':
    test()