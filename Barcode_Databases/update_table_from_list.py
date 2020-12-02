import pandas as pd
import sys

input_list = sys.argv[1]
input_table = sys.argv[2]

ls_df = pd.read_table(input_list, header=None)
table_df = pd.read_csv(input_table)

print(ls_df.shape)
print(table_df.shape)
table_df2 = table_df[table_df.SeqID.isin(ls_df[0])]
print(table_df2.shape)
table_df2.to_csv(input_table.replace('.csv','_update.csv'),index=False)
