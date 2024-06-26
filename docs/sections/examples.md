## CellAgent surpasses traditional analysing process in its advanced automation capabilities on single-cell data.

<div class="mb-4"></div>

<v-container class="mb-16 py-0 px-0">
  <v-row>
    <v-col cols="6">
      <div class="mb-2 text-h5">Traditional scRNA-seq analysis</div>
      <v-card
        border="surface-variant sm opacity-50"
        variant="text"
        class="rounded-0"
      >
        <v-card-text>
● Require programming skills.<br>
● Require biological expertise.<br>
● Adjust hyperparameters manually.<br>
        </v-card-text>
      </v-card>
      <v-card
        border="surface-variant sm opacity-50"
        variant="text"
        class="rounded-0 border-t-0"
      >
        <v-card-text>
        <chat-bubble>
            <template v-slot:subject><v-img src="/user.png" /></template>
<div>

```python
# Filter based on the properties of current dataset and empirical thresholds
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
# Calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells with less than 2500 genes and mitochondrial gene percentage less than 5%
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# Other steps omitted

sc.tl.pca(adata, svd_solver='arpack')
# For this dataset, n_neighbors=10 and n_pcs=40 show clear results
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['leiden'])
```
</div>
        </chat-bubble>

<!-- <chat-bubble> -->
  <!-- <template v-slot:subject><v-img src="/teacher.png" /></template> -->
  <img src="/teacher.png" alt="" style="height: 25px; display: inline-block; ">

  <span class="bg-grey-lighten-2">
    "Based on the current data analysis results, and in biomedical research, these genes can typically serve as markers to help us identify different cell types..."
  </span><br>
<!-- </chat-bubble> -->
<chat-bubble>
            <template v-slot:subject><v-img src="/user.png" /></template>
<div>

```python
# Find differentially expressed genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

#After comparing the differentially expressed genes in each cluster, and under the guidance of experts, cell type labels were assigned to each cluster.
cluster_annotations = {
    '0': 'T cells','1': 'B cells','2': 'NK cells','3': 'Monocytes','4': 'Dendritic cells'...
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)
```
</div>
  <span class="bg-grey-lighten-2">
    📄 "Finally，we get the annotation result."
  </span>
</chat-bubble>
        </v-card-text>
      </v-card>
    </v-col>
    <v-col cols="6">
      <div class="mb-2 text-h5 text-success">CellAgent</div>
      <v-card
        border="success sm opacity-50"
        variant="text"
        class="rounded-0 text-success"
      >
        <v-card-text>
● Code less, dialogue only.<br>
● Collaborate with LLM biological experts.<br>
● Optimize result automatically.<br>
        </v-card-text>
      </v-card>
      <v-card
        border="success sm opacity-50"
        variant="text"
        class="rounded-0 border-t-0 text-success"
      >
        <v-card-text>
          <chat-bubble>
            <template v-slot:subject><v-img src="/user.png" /></template>
            <span class="bg-grey-lighten-2">
              "This is raw PBMC dataset. Please help me complete cell type annotation."<br>
            </span>
          </chat-bubble>
          <chat-bubble>
            <template v-slot:subject><v-img src="/logo.png" /></template>
            <span class="bg-light-green-lighten-2">
              After an initial anlysis, 6 steps are required to complete your request:<br><br>
              1.Quality Control<br>
              2. Normalization<br>
              3. Identification of Highly<br>
              4. Dimensionality Reduction<br>
              5. Clustering<br><br>
            </span>
            <span class="bg-light-green-lighten-2">
              ...
            </span>
          </chat-bubble>
          <chat-bubble>
            <template v-slot:subject><v-img src="/user.png" /></template>
            <span class="bg-grey-lighten-2">
              "Please using the Leiden algorithm to cluster."
            </span>
          </chat-bubble>
          <!-- <chat-bubble> -->
            <span>
              ...<br><br>
            </span>
          <!-- </chat-bubble> -->
          <chat-bubble>
            <template v-slot:subject><v-img src="/user.png" /></template>
            <span class="bg-grey-lighten-2">
              "Great, help me complete the cell type annotation finally."
            </span>
          </chat-bubble>
          <chat-bubble>
            <template v-slot:subject><v-img src="/logo.png" /></template>
            <span class="bg-light-green-lighten-2">
              Combining differential gene expression across clusters and results from cell type annotation tools like Celltypist, the cell labels for these clusters were finally confirmed and saved as `.obs['final_type']` after evaluation.<br><br>
            </span>
<div>

<img src="/cellexample2.png" alt="" style="height: 200px; display: inline-block;"><br>
</div>
          </chat-bubble>
        </v-card-text>
      </v-card>
    </v-col>
  </v-row>
</v-container>
