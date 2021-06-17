import pandas as pd 

df_used_index = pd.read_excel('../3.figures_files/labels.xlsx', index_col=1)
index_core_cluster = df_used_index.index.values
df_FMT_top20 = pd.read_excel('../1.rawdata/FMT_OTUs_20210612.xlsx', index_col=0)
index_FMT_top20 = df_FMT_top20.index.values[:20]
index_FMT_top20_STR = []
for index_name in index_FMT_top20:
    index_FMT_top20_STR.append(str(index_name))

index_selected = list(set(list(index_core_cluster)+list(index_FMT_top20_STR)))


df_OTU_replicates = pd.read_excel('../2.output_files/OTU_table_20210611.xlsx', index_col=0)
df_OTU_selected = df_OTU_replicates.loc[index_selected]
df_OTU_selected.to_excel('../2.output_files/OTUs4cooccurance.xlsx', index=True)