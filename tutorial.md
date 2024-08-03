This tutorial is based on the [BayesPrism Tutorial](https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.pdf). It is tested with Python 3.8, but newer versions should work.

### Create and Use Python Virtual Environment venv

```bash
python3 -m venv venv
source venv/bin/activate
```

### Install pybayesprism with pip
```bash
python3 -m pip install pybayesprism==0.1.0
``` 

You should see this output if installation is successful:
```
...
Successfully installed numpy-1.24.4 packaging-24.1 pandas-2.0.3 pybayesprism-0.1.0 python-dateutil-2.9.0.post0 pytz-2024.1 scipy-1.10.1 six-1.16.0 tqdm-4.66.4 tzdata-2024.1 xarray-2023.1.0
```

### Import Packages and Download Tutorial Data

```python
# import necessary packages
import os
import pandas as pd
import numpy as np
from pybayesprism import process_input, prism, extract

# download tutorial file
os.system("curl -L -O https://github.com/ziluwang829/pyBayesPrism/raw/main/data/tutorial.tar.gz")
os.system("mkdir -p BP_data")
os.system("tar -xzvf tutorial.tar.gz -C BP_data")
```

If untar is successful,
```
x bk_dat.csv
x sc_dat.csv
x cell_state_labels.csv
x cell_type_labels.csv
0
```

### Load Data and Check Data

Note that this step will take a while.
```python
bk_dat = pd.read_csv("BP_data/bk_dat.csv", sep=",", index_col=0).astype(np.int32)
sc_dat = pd.read_csv("BP_data/sc_dat.csv", sep=",", index_col=0).astype(np.int32)

cell_state_labels = pd.read_csv("BP_data/cell_state_labels.csv", header=None).iloc[:,0].tolist()
cell_type_labels = pd.read_csv("BP_data/cell_type_labels.csv", header=None).iloc[:,0].tolist()
```

Now, check the dataframes we just loaded.
#### bk_dat

`bk_dat` : The sample-by-gene raw count pandas dataframe of bulk RNA-seq expression. Index are bulk sample IDs, while columns are gene names/IDs.

```python
bk_dat.shape
```
```
(169, 60483)
```
<br/>

```python
bk_dat.index
```
```
Index(['TCGA-06-2563-01A-01R-1849-01', 'TCGA-06-0749-01A-01R-1849-01',
       'TCGA-06-5418-01A-01R-1849-01', 'TCGA-06-0211-01B-01R-1849-01',
       'TCGA-19-2625-01A-01R-1850-01', 'TCGA-19-4065-02A-11R-2005-01',
       'TCGA-06-0211-02A-02R-2005-01', 'TCGA-15-0742-01A-01R-1850-01',
       'TCGA-06-0138-01A-02R-1849-01', 'TCGA-06-0645-01A-01R-1849-01',
       ...
       'TCGA-12-3650-01A-01R-1849-01', 'TCGA-06-2564-01A-01R-1849-01',
       'TCGA-06-0644-01A-02R-1849-01', 'TCGA-27-2523-01A-01R-1850-01',
       'TCGA-14-0817-01A-01R-1849-01', 'TCGA-06-0168-01A-01R-1849-01',
       'TCGA-02-0047-01A-01R-1849-01', 'TCGA-32-2638-01A-01R-1850-01',
       'TCGA-19-4065-01A-01R-2005-01', 'TCGA-76-4926-01B-01R-1850-01'],
      dtype='object', length=169)
```
<br/>

```python
bk_dat.columns
```
```
Index(['ENSG00000000003.13', 'ENSG00000000005.5', 'ENSG00000000419.11',
       'ENSG00000000457.12', 'ENSG00000000460.15', 'ENSG00000000938.11',
       'ENSG00000000971.14', 'ENSG00000001036.12', 'ENSG00000001084.9',
       'ENSG00000001167.13',
       ...
       'ENSGR0000263980.4', 'ENSGR0000264510.4', 'ENSGR0000264819.4',
       'ENSGR0000265658.4', 'ENSGR0000270726.4', 'ENSGR0000275287.3',
       'ENSGR0000276543.3', 'ENSGR0000277120.3', 'ENSGR0000280767.1',
       'ENSGR0000281849.1'],
      dtype='object', length=60483)
```
<br/>

#### sc_dat

`sc_dat` : The cell-by-gene raw count pandas dataframe of single cell RNA-seq expression. Index are cell IDs, while columns are gene names/IDs. `sc_dat` should be a dense matrix.

```python
sc_dat.shape
```
```
(23793, 60294)
```
<br/>

```python
sc_dat.index
```
```
Index(['PJ016.V3', 'PJ016.V4', 'PJ016.V5', 'PJ016.V6', 'PJ016.V7', 'PJ016.V8',
       'PJ016.V9', 'PJ016.V10', 'PJ016.V11', 'PJ016.V12',
       ...
       'PJ048.V3077', 'PJ048.V3078', 'PJ048.V3079', 'PJ048.V3080',
       'PJ048.V3081', 'PJ048.V3082', 'PJ048.V3083', 'PJ048.V3084',
       'PJ048.V3085', 'PJ048.V3086'],
      dtype='object', length=23793)
```
<br/>

```python
sc_dat.columns
```
```
Index(['ENSG00000130876.10', 'ENSG00000134438.9', 'ENSG00000269696.1',
       'ENSG00000230393.1', 'ENSG00000266744.1', 'ENSG00000202281.1',
       'ENSG00000165244.6', 'ENSG00000275916.1', 'ENSG00000173597.7',
       'ENSG00000158022.6',
       ...
       'ENSG00000182798.9', 'ENSG00000138829.9', 'ENSG00000204614.7',
       'ENSG00000173214.5', 'ENSG00000273451.1', 'ENSG00000280003.1',
       'ENSG00000222979.1', 'ENSG00000242860.3', 'ENSG00000249321.1',
       'ENSG00000174796.11'],
      dtype='object', length=60294)
```
<br/>

Note that BayesPrism will perform deconvolution over genes shared between the mixture and scRNA- seq reference, i.e., by intersecting index(mixture) and columns(reference). Therefore, please make sure to use consistent gene annotations (TCGA RNA-seq uses GENCODE v22).

We recommend the use of unnormalized and untransformed count data. When raw count is not available, linear normalization, such as TPM, RPM, RPKM, FPKM, is also acceptable, as BayesPrism is robust to linear multiplicative difference between the reference and mixture. Ideally, if using normalized data, it is the best to supply reference and bulk transformed by the same method. log transformation should be avoided.

#### cell_state_labels

`cell_type_labels` : a list of the same length as row count of `sc_dat` to denote the cell type of each cell in the reference.

```python
len(cell_state_labels)
```
```
23793
```
<br/>

```python
pd.Series(cell_state_labels).value_counts()
```
```
PJ025-tumor-0    971
PJ025-tumor-1    941
PJ025-tumor-2    630
PJ016-tumor-0    621
PJ016-tumor-1    619
                ... 
PJ032-tumor-3     62
PJ032-tumor-4     57
myeloid_8         49
PJ032-tumor-5     41
PJ017-tumor-6     22
Name: count, Length: 73, dtype: int64
```
<br/>

#### cell_type_labels

`cell.state_labels`  a list of the same length as row count of `sc_dat` to denote the cell state of each cell in the reference. In our example, cell states of malignant cells were obtained by sub-clustering the malignant cells from each patient, and cell states of myeloid cells were obtained by clustering myeloid cells from all patients. We de"ne multiple cell states for these two cell types, as they contain substantial heterogeneity while also having su$cient number of cells for sub-clustering.


```python
len(cell_type_labels)
```
```
23793
```
<br/>

```python
pd.Series(cell_type_labels).value_counts()
```
```
tumor          20080
myeloid         2505
endothelial      492
pericyte         489
oligo            160
tcell             67
Name: count, dtype: int64
```
<br/>


Please make sure that all cell states contain a reasonable number of cells, e.g. >20 or >50, so that their profile can be represented accurately.

What to supply for `cell_state_labels` and `cell_type_labels`? The definition of cell type and cell state can be somewhat arbitrary (similar to the issue of assigning cell types for scRNA-seq) and depends on the question of interest. Their definitions depend on the granularity we aim at and the confidence of the `cell_type_labels` in scRNA-seq data. Usually, a good rule of thumb is as follows. 

1) Define cell types as the cluster of cells having a sufficient number of significantly differentially expressed genes than other cell types, e.g., greater than 50 or even 100. For clusters that are too similar in transcription, we recommend treating them as cell states, which will be summed up before the final Gibbs sampling. Therefore, cell states are often suitable for cells that form a continuum on the phenotypic manifold rather than distinct clusters. 
2) Define multiple cell states for cell types of significant heterogeneity, such as malignant cells, and of interest to deconvolve their transcription.

### Process Input

Next, we remove the genes from selected groups. Note that when sex is not identical between the reference and mixture, we recommend excluding genes from chrX and chrY. We also remove lowly transcribed genes, as the measurement of transcription of these genes tend to be noise-prone. Removal of these genes can also speed up computation.

```python
sc_dat_filtered = process_input.cleanup_genes(sc_dat, "count.matrix", "hs", \
                  ["Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"], 5)
```
Expected output is:
```
EMSEMBLE IDs detected.
number of genes filtered in each category: 
Rb            89
Mrp           78
other_Rb    1011
chrM          37
MALAT1         1
chrX        2464
chrY         594
dtype: int64
A total of 4214 genes from ['Rb', 'Mrp', 'other_Rb', 'chrM', 'MALAT1', 'chrX', 'chrY'] have been excluded
A total of 24343 gene expressed in fewer than 5 cells have been excluded
```
<br/>

Check the dimensions of `sc_dat_filtered`:
```python
sc_dat_filtered.shape
```
```
(23793, 31737)
```
<br/>


We observe that protein coding genes are the most concordant group between two assays. To reduce batch effects and speed up computation, we perform deconvolution on protein coding genes. We have also tried to run BayesPrism on all genes. The results were similar.

```python
sc_dat_filtered_pc = process_input.select_gene_type(sc_dat_filtered, ["protein_coding"])
```
Expected output is:
```
EMSEMBLE IDs detected.
number of genes retained in each category: 
category
protein_coding    16148
Name: count, dtype: int64
```

### Construct a Prism Object

A prism object contains all data required for running BayesPrism, namely, a scRNA-seq reference matrix, the cell type and cell state labels of each row of reference, and the mixture matrix for bulk RNA-seq.

When using scRNA-seq count matrix as the input (recommended), user needs to specify `input_type = "count.matrix"`. The other option for `input_pe` is `"GEP"` (gene expression profile) which is a cell state by gene matrix. This option is used when using reference derived from other assays, such as sorted bulk data.

The parameter key is a character in `cell_type_labels` that corresponds to the malignant cell type. Set to `None` if there are no malignant cells or the malignant cells between reference and mixture are from matched samples, in which case all cell types will be treated equally.


```python
my_prism = prism.Prism.new(reference = sc_dat_filtered_pc, 
                          mixture = bk_dat, input_type = "count.matrix", 
                          cell_type_labels = cell_type_labels, 
                          cell_state_labels = cell_state_labels, 
                          key = "tumor", 
                          outlier_cut = 0.01, 
                          outlier_fraction = 0.1)
```

Expected output is:
```
PJ025-tumor-0    971
PJ025-tumor-1    941
PJ025-tumor-2    630
PJ016-tumor-0    621
PJ016-tumor-1    619
                ... 
PJ032-tumor-3     62
PJ032-tumor-4     57
myeloid_8         49
PJ032-tumor-5     41
PJ017-tumor-6     22
Name: count, Length: 73, dtype: int64
Number of outlier genes filtered from mixture = 6
Aligning reference and mixture...
Normalizing reference...
```

### Run BayesPrism

Next, we start the run of BayesPrism.

```python
bp_res = my_prism.run(n_cores = 10, update_gibbs = True)
```
Your output should look similar to:
```
Run Gibbs sampling...
Current time:  2024-08-01 22:48:05.832530
Estimated time to complete:  59mins
Estimated finishing time:  2024-08-01 23:46:16.605730
Start run...
100%|█████████████████████████████████████████| 169/169 [27:31<00:00,  9.77s/it]
Now Merging...
Update the reference matrix ...
running with 10 cores!
Run Gibbs sampling using updated reference ...
Current time:  2024-08-01 23:22:23.613335
Estimated time to complete:  24mins
Estimated finishing time:  2024-08-01 23:46:00.420535
Start run...
```

### Extract the Results
```python
theta = extract.get_fraction(bp_res, "final", "type")
theta_cv = bp_res.posterior_theta_f.theta_cv
Z_tumor = extract.get_exp(bp_res, "type", "tumor")
```

theta:
```
>>> theta
                                 tumor   myeloid  ...         tcell     oligo
TCGA-06-2563-01A-01R-1849-01  0.835802  0.043125  ...  6.473308e-04  0.011552
TCGA-06-0749-01A-01R-1849-01  0.706223  0.169320  ...  1.359410e-06  0.107721
TCGA-06-5418-01A-01R-1849-01  0.859105  0.098019  ...  5.512243e-07  0.005164
TCGA-06-0211-01B-01R-1849-01  0.885878  0.044639  ...  3.478573e-06  0.000131
TCGA-19-2625-01A-01R-1850-01  0.936858  0.035302  ...  2.422870e-06  0.008850
...                                ...       ...  ...           ...       ...
TCGA-06-0168-01A-01R-1849-01  0.760788  0.136781  ...  2.817239e-04  0.031485
TCGA-02-0047-01A-01R-1849-01  0.784189  0.100653  ...  1.997586e-06  0.026175
TCGA-32-2638-01A-01R-1850-01  0.781563  0.106922  ...  2.376922e-04  0.004507
TCGA-19-4065-01A-01R-2005-01  0.726260  0.193135  ...  4.569846e-05  0.015416
TCGA-76-4926-01B-01R-1850-01  0.892356  0.041993  ...  1.174681e-05  0.030072

[169 rows x 6 columns]
```

theta_cv:
```
>>> theta_cv
                                 tumor   myeloid  pericyte  endothelial     tcell     oligo
TCGA-06-2563-01A-01R-1849-01  0.063373  0.063394  0.063446     0.063392  0.077282  0.063601
TCGA-06-0749-01A-01R-1849-01  0.063373  0.063375  0.763787     0.063562  1.045820  0.063383
TCGA-06-5418-01A-01R-1849-01  0.063373  0.063380  0.063744     0.063457  0.837628  0.064056
TCGA-06-0211-01B-01R-1849-01  0.063373  0.063391  0.063600     0.063395  0.574891  0.285020
TCGA-19-2625-01A-01R-1850-01  0.063373  0.063408  0.068472     0.063611  0.685201  0.063851
...                                ...       ...       ...          ...       ...       ...
TCGA-06-0168-01A-01R-1849-01  0.063374  0.063391  0.063715     0.063518  0.119633  0.063630
TCGA-02-0047-01A-01R-1849-01  0.063373  0.063381  0.063453     0.063429  0.735219  0.063553
TCGA-32-2638-01A-01R-1850-01  0.063373  0.063377  0.063419     0.063392  0.089633  0.063845
TCGA-19-4065-01A-01R-2005-01  0.063373  0.063376  0.063421     0.063429  0.280551  0.063582
TCGA-76-4926-01B-01R-1850-01  0.063373  0.063394  0.065483     0.063431  0.386786  0.063443

[169 rows x 6 columns]
```

Z_tumor:
```
>>> Z_tumor
<xarray.DataArray (bulk_id: 169, gene_id: 16145)>
array([[  55.776,  205.412,   17.316, ...,  957.188,  131.228,  503.008],
       [ 443.148,   38.036,    6.544, ...,   26.972,   84.056,  302.932],
       [   8.964,  291.592,   21.704, ...,   89.404,  185.356,  581.952],
       ...,
       [  20.908,  350.384,    6.216, ..., 1259.216,  233.484,  565.92 ],
       [  87.62 ,  428.596,   25.804, ...,  329.572,  162.708,  557.768],
       [ 225.092,  384.988,    9.136, ..., 1821.112,  244.676,  499.212]])
Coordinates:
  * bulk_id           (bulk_id) object 'TCGA-06-2563-01A-01R-1849-01' ... 'TC...
  * gene_id           (gene_id) object 'ENSG00000130876.10' ... 'ENSG00000174...
    cell_type_merged  <U11 'tumor'
```

### Clean Up Environment
Exit Python, deactivate venv, and remove the files.
```shell
deactivate
rm -rf venv
rm -rf BP_data
rm tutorial.tar.gz
```
