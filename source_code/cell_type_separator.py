import pandas as pd
import numpy as np
import os

file_path = './data/'
file = 'Normal_lung__I_log2TPM_matrix.txt'
file_annotation = 'Normal_lung__I_annotation.txt'

separator = 'Cell_type'

annotation = pd.read_csv(file_path + file_annotation, sep=',')
print(annotation)
separated_index = {annotation['Index'][i]:annotation[separator][i] for i in range(len(annotation))}
sep_values = list(separated_index.values())

#getting columns in data
columns = ''
with open(file_path + file) as data:
  columns = data.readline()
  data.close()
file_list = {}
for val in sep_values:
  if val not in file_list:
    val_formated = val.replace(' ', '_')
    if not os.path.exists(f'{file_path}{val_formated}'):
      os.mkdir(f'{file_path}{val_formated}')
    file_list[val] = open(f'{file_path}{val_formated}/{file}', 'w')
    file_list[val].writelines(columns)

print('ok')
with open(file_path + file) as data:
  i = 0
  batch = {val:[] for val in sep_values}
  columns = data.readline()
  line = data.readline()
  while line:
    values = line.replace('\n', '').split('\t')
    index = values[0]
    sep_val = separated_index[index]
    batch[sep_val] += [line]
    if i % 100 == 0:
      for val in sep_values:
        file_list[val].writelines(batch[val])
        batch[val] = []
    line = data.readline()
    i += 1

# writting last batch
for val in sep_values:
  file_list[val].writelines(batch[val])



