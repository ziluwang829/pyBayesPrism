This is a Python implementation of [BayesPrism](https://github.com/Danko-Lab/BayesPrism).


Usage
```python

import os
import pandas as pd
from pybayesprism import process_input, prism

os.system("curl -L -O https://github.com/ziluwang829/pyBayesPrism/raw/main/data/data.tar.gz")
os.system("mkdir -p BP_data")
os.system("tar -xzvf data.tar.gz -C BP_data")

bk_dat = pd.read_csv("BP_data/bk_dat.csv", sep=",", index_col=0)
sc_dat = pd.read_csv("BP_data/sc_dat.csv", sep=",", index_col=0)


cell_state_labels = pd.read_csv("BP_data/cell_state_labels.csv", header=None).iloc[:,0].tolist()

cell_type_labels = pd.read_csv("BP_data/cell_type_labels.csv", header=None).iloc[:,0].tolist()

sc_dat_filtered = process_input.cleanup_genes(sc_dat, "count.matrix", "hs", \
                  ["Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"], 5)
                  
sc_dat_filtered_pc = process_input.select_gene_type(sc_dat_filtered, ["protein_coding"])

my_prism = prism.Prism.new(reference = sc_dat_filtered_pc, 
                          mixture = bk_dat, input_type = "count.matrix", 
                          cell_type_labels = cell_type_labels, 
                          cell_state_labels = cell_state_labels, 
                          key = "tumor", 
                          outlier_cut = 0.01, 
                          outlier_fraction = 0.1)

bp_res = my_prism.run(n_cores = 36, update_gibbs = True)      
```
