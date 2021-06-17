from numpy.core.numeric import True_
import pandas as pd 
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler


 
def get_keys(d, value):
    return [k for k,v in d.items() if v == value]


def forest_selection(df_in_path='../2.output_files/OTU_table_20210611.xlsx'):
    metadata = pd.read_table('../1.data/metadata.txt', index_col=0)
    dict_sampleid2sampletype = {}
    for index_id in metadata.index.values:
        dict_sampleid2sampletype[index_id] = metadata.loc[index_id, 'SampleType']
    print(dict_sampleid2sampletype)
    # 开始对原始数据进行处理
    df = pd.read_excel(df_in_path, index_col=0)
    sample_types_all = list(set(dict_sampleid2sampletype.values()))
    for sample_types in sample_types_all:
        columns_using = [column_names for column_names in df.columns.values if sample_types in column_names]
        df4keys = df[columns_using]
        df4keys.to_excel('../2.output_files/%s_20210610.xlsx'%sample_types, index=True)
    pass


#StandardScaler按列z-score
def optimization(df_control, df_treatment, name_output=''):
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split
    # 数据表格样本为行
    X_control = df_control.values
    X_treatment = df_treatment.values
    X = np.concatenate((X_control, X_treatment), axis=0)
    Y = []
    for i in range(len(df_control)):
        Y.append(0)
    for i in range(len(df_treatment)):
        Y.append(1)
    print(len(X))
    print(len(Y))
    s1 = StandardScaler()
    X = s1.fit_transform(X)
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(X,Y,test_size=0.2)
    rfc = RandomForestClassifier(random_state=30)
    rfc = rfc.fit(Xtrain,Ytrain)
    score_r = rfc.score(Xtest,Ytest)
    print("Random Forest:{}".format(score_r))
    # 特征重要性提取
    features = df_control.columns.values
    feature_importances = rfc.feature_importances_
    features_df = pd.DataFrame({'Features':features,'Importance':feature_importances})
    features_df.sort_values('Importance',inplace=True,ascending=False)
    features_df.to_excel('../3.figures_files/RandomForestTree_featureImportance_%s.xlsx'%name_output, index=False)
    import seaborn as sns
    
    sns.set(rc={"figure.figsize": (8, 8)})
    ax = sns.barplot(features_df['Features'][:15], features_df['Importance'][:15],)
    plt.ylabel('Importance')
    # 数据可视化：柱状图
    sns.despine(bottom=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.figure.tight_layout()
    plt.savefig('../3.figures_files/%s_20210612.pdf'%name_output)
    plt.close()
    # plt.show()
    return 0


def randomforest_modeling():
    df_D28_JFH_CCl4 = pd.read_excel('../2.output_files/D28-JFH-CCl4_20210610.xlsx', index_col=0).T
    df_D28_JFD_CCl4 = pd.read_excel('../2.output_files/D28-JFD-CCl4_20210610.xlsx', index_col=0).T
    df_D28_CTL = pd.read_excel('../2.output_files/D28-CTL_20210610.xlsx', index_col=0).T
    df_D28_CCl4 = pd.read_excel('../2.output_files/D28-CCl4_20210610.xlsx', index_col=0).T
    df_D14_JFH_CCl4 = pd.read_excel('../2.output_files/D14-JFH-CCl4_20210610.xlsx', index_col=0).T
    df_D14_JFD_CCl4 = pd.read_excel('../2.output_files/D14-JFD-CCl4_20210610.xlsx', index_col=0).T
    df_D14_CTL = pd.read_excel('../2.output_files/D14-CTL_20210610.xlsx', index_col=0).T
    df_D14_CCl4 = pd.read_excel('../2.output_files/D14-CCl4_20210610.xlsx', index_col=0).T
    df_D1_JFH_CCl4 = pd.read_excel('../2.output_files/D1-JFH-CCl4_20210610.xlsx', index_col=0).T
    df_D1_JFD_CCl4 = pd.read_excel('../2.output_files/D1-JFD-CCl4_20210610.xlsx', index_col=0).T
    df_D1_CTL = pd.read_excel('../2.output_files/D1-CTL_20210610.xlsx', index_col=0).T
    df_D1_CCl4 = pd.read_excel('../2.output_files/D1-CCl4_20210610.xlsx', index_col=0).T
    df_JFH = df_D14_JFH_CCl4.append(df_D28_JFH_CCl4)
    df_JFD = df_D28_JFD_CCl4.append(df_D14_JFD_CCl4)
    df_CCl4 = df_D14_CCl4.append(df_D28_CCl4)
    df_CTL = df_D28_CTL.append(df_D14_CTL)
    optimization(df_control=df_CCl4, df_treatment=df_JFD, name_output='JFD2CCl4')
    optimization(df_control=df_CCl4, df_treatment=df_JFH, name_output='JFH2CCl4')
    optimization(df_control=df_CCl4, df_treatment=df_CTL, name_output='CTL2CCl4')
    return 0
   

if __name__ == '__main__':
    forest_selection()
    randomforest_modeling()