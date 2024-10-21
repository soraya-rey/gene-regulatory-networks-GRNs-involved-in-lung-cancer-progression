import pandas as pd
import numpy as np
import os

"""Processing script: separate data by cell type, sample source and tumor stage"""

columns_anotation_file = 'GSE131907_Lung_Cancer_cell_annotation.txt'
cell_data_file = 'cell_data_simplified.csv'
data_file = 'GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt'

# 1. Getting each cell destination
# 1.1 Fiding samples cell type
annotation = pd.read_csv(columns_anotation_file, sep='\t')
cells_type = {annotation['Index'][i]:annotation['Cell_type'][i].replace(' ', '_') for i in range(len(annotation))}

# 1.2 Fiding cells corresponding sample
cells_sample = {annotation['Index'][i]:annotation['Sample'][i] for i in range(len(annotation))}

# 1.3 Fiding samples stages
sample_data = pd.read_csv(cell_data_file, sep=',')
samples_stage = {sample_data['Sample_title'][i]:sample_data['tumor stage'][i].replace(' ', '') for i in range(len(sample_data))}

# 1.4 Fiding samples sources
samples_source = {sample_data['Sample_title'][i]:sample_data['Sample_source_name_ch1'][i].replace(' ', '_') for i in range(len(sample_data))}

# 1.5 Writing cells destination file's title
cells_file = {}
for cell in cells_sample.keys():
  sample = cells_sample[cell]
  cells_file[cell] = f'./data/{cells_type[cell]}/{samples_source[sample]}_{samples_stage[sample]}'
destination_files = np.unique(list(cells_file.values()))


# 2. Getting indexes of cells in RNAseq-file's header
# 2.1 Reading column line
columns = ''
with open(data_file) as data:
  columns = data.readline()
  data.close()

# 2.2 Parsing columns
parsed_columns = columns.replace('\n','').split('\t')

# 2.4 Accounting for Index column at 0, calculating indexes for each cell properties combination
indexes = {dest:[0] for dest in destination_files}
for i in range(1, len(parsed_columns)):
  cell = parsed_columns[i]
  dest = cells_file[cell]
  indexes[dest] += [i]

# 3. Writting data in files
# 3.1 Creating folders
if not os.path.exists('./data'):
  os.mkdir('./data')

for type in np.unique(list(cells_type.values())):
  folder = f'./data/{type}'
  if not os.path.exists(folder):
    os.mkdir(folder)

# 3.2 Opening files and writting header
opened_files = {
  dest: open(f'{dest}_matrix.csv', 'w') for dest in destination_files
}

# 3.3 Going through data and writting each cell ensemble in defined file
with open(data_file) as data:
  # 3.3.a line initialisation
  line = 42

  # 3.3.b Batch initialisation to limit in/out operations
  batch = {dest: [] for dest in destination_files}
  i = 0

  # 3.3.c reading each line and for each subset, selecting values and writting in batch
  while line:
    line = data.readline()
    if len(line) == 0:
      break
    parsed_line = np.array(line.replace('\n', '').split('\t'))
    for dest in indexes.keys():
      subset = parsed_line[indexes[dest]]
      batch[dest] += [','.join(subset) + '\n']
    i += 1

    # When batch is full, write in files and empty batch
    if i % 100 == 0:
      print(f'Lines processed: {i}')
      for dest in batch.keys():
        opened_files[dest].writelines(batch[dest])
        batch[dest] = []

  # 3.3.d Writting last batched lines
  print(f'Total lines processed: {i}')
  for dest in batch.keys():
    opened_files[dest].writelines(batch[dest])

# 3.4 (GOOD PRACTICE) Closing files
for file in opened_files.keys():
  opened_files[file].close()








