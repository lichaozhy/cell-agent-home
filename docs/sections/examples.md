### CellAgent surpasses GPT-4 in its advanced automation capabilities on single-cell analysis.

<v-container class="mt-8 mb-16 py-0 px-8">
  <v-row>
    <v-col cols="6">
      <h3 class="mb-2">Traditional Method</h3>
      <v-card
        border="surface-variant sm opacity-100"
        title=""
        variant="text"
        class="rounded-0"
      >
        <v-card-text>
1.Researchers are required to have programming skills ⌨️<br>
2.Expert knowledge and searches of relevant literature databases📚<br>
3.Alternatively, use automated annotation tools📐<br>
        </v-card-text>
      </v-card>
      <v-card
        border="surface-variant sm opacity-100"
        variant="text"
        class="rounded-0 border-t-0"
      >
        <v-card-text>
        <chat-bubble>
            <template v-slot:subject><v-img src="/user.png" /></template>
<div>

```python
import scanpy as sc
import pandas as pd
adata = sc.datasets.pbmc3k()
print(adata)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['leiden'])
```
</div>
        </chat-bubble>

<chat-bubble>
  <template v-slot:subject><v-img src="/user.png" /></template>
  <span class="bg-grey-lighten-2">
    "Based on current single-cell RNA sequencing research and literature data, we can use the following marker genes to annotate."
  </span>
<div>

```python
marker_genes = {
    'T_cells': ['IL7R', 'CD3D'],
    'B_cells': ['CD79A', 'MS4A1'],
    'NK_cells': ['GNLY', 'NKG7'],
    'Monocytes': ['LYZ', 'CD14'],
    'Dendritic_cells': ['FCER1A', 'CST3'],
    'Megakaryocytes': ['PPBP']
}

sc.tl.score_genes_cell_cycle(adata, s_genes=marker_genes['T_cells'], g2m_genes=marker_genes['B_cells'])
sc.pl.umap(adata, color=['CD3D', 'CD79A', 'GNLY', 'LYZ', 'FCER1A', 'PPBP'])

adata.obs['cell_type'] = adata.obs['leiden'].astype(str)

for cell_type, genes in marker_genes.items():
    adata.obs.loc[adata[:, genes].X.mean(1) > 0.5, 'cell_type'] = cell_type

sc.pl.umap(adata, color=['cell_type'])
```
</div>
  <span class="bg-grey-lighten-2">
    "📄 Finally，we get the annotation result."
  </span>
</chat-bubble>
        </v-card-text>
      </v-card>
    </v-col>
    <v-col cols="6">
      <h3 class="mb-2 text-success">CellAgent</h3>
      <v-card
        border="success sm opacity-100"
        title=""
        variant="text"
        class="rounded-0 text-success"
      >
        <v-card-text>
1. Allows users to perform various single-cell data analyses through dialogue, completely code-free. <br>
2. Multi-Agent architecture, leveraging collaboration among multiple experts to accomplish data analysis tasks.<br>
3. Self-iterative evaluation effectively enhances the quality of data analysis.<br>
        </v-card-text>
      </v-card>
      <v-card
        border="success sm opacity-100"
        title=""
        variant="text"
        class="rounded-0 border-t-0 text-success"
      >
        <v-card-text>
          <chat-bubble>
            <span class="bg-grey-lighten-2">
              This is raw PBMC dataset. Please help me complete cell type annotation.<br>
            </span>
          </chat-bubble>
          <chat-bubble>
            <template v-slot:subject><v-img src="/logo.png" /></template>
            <span class="bg-light-green-lighten-2">
              After an initial anlysis, 6 steps are required to complete your request:<br>
              1.Quality Control<br>
              2. Normalization<br>
              3. Identification of Highly<br>
              4. Dimensionality Reduction<br>
              5. Clustering<br>
              6. Clustering<br>
            </span>
            <span class="bg-light-green-lighten-2">
              ...
            </span>
          </chat-bubble>
          <chat-bubble>
            <span class="bg-grey-lighten-2">
              Please using the Leiden algorithm to cluster.
            </span>
          </chat-bubble>
          <chat-bubble>
            <template v-slot:subject><v-img src="/logo.png" /></template>
            <span class="bg-light-green-lighten-2">
              I will first compute the neighborhood graph and then apply the Leiden algorithm to identify clusters. Finally, visualize the clustering results using UMAP.
            </span>
<div>

![An image](/cellexample1.png)
</div>
          </chat-bubble>
          <chat-bubble>
            <span class="bg-grey-lighten-2">
              Great, help me complete the cell type annotation finally.
            </span>
          </chat-bubble>
          <chat-bubble>
            <template v-slot:subject><v-img src="/logo.png" /></template>
            <span class="bg-light-green-lighten-2">
              Combining differential gene expression across clusters and results from cell type annotation tools like Celltypist, the cell labels for these clusters were finally confirmed and saved as `.obs['final_type']` after evaluation.
            </span>
<div>

![An image](/cellexample2.png)
</div>
          </chat-bubble>
        </v-card-text>
      </v-card>
    </v-col>
  </v-row>
</v-container>
