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
        title="Output"
        variant="text"
        class="rounded-0 border-t-0"
      >
        <v-card-text style="text-wrap-mode:wrap;white-space:pre-wrap;white-space-collapse:preserves">
          
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

![An image](/user.png)
<span class="bg-grey-lighten-2">"Based on current single-cell RNA sequencing research and literature data, we can use the following marker genes to annotate..."
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
<span class="bg-grey-lighten-2">"📄 Finally，we get the annotation result."
</span>
<!-- <div>

::: danger
GPT-4 provided code for a basic solution.
:::
</div> -->
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
          <v-expansion-panels tile elevation="0">
            <v-expansion-panel title="Thought and code generation">
              <v-expansion-panel-text>
                <span class="bg-light-green-lighten-2">
                  This is raw PBMC dataset. Please help me complete cell type annotation.<br>
                </span>
<div>

::: code-group

```python [Basic]
# Cell Type Annotation using AnnotatorCellmarkerACT
annotator = AnnotatorCellmarkerACT()
adata = annotator.run(species='Human', tissue_type='Blood', adata=adata, obs_cluster='leiden')
...
```

```python [Optimized]
# To optimize the cell type annotation step, let's use the
# `AnnotatorCelltypist` tool and the "Immune_All_Low.pkl" model,
# which is suitable for immune sub-populations.
annotator = AnnotatorCelltypist()
adata = annotator.run(model_name='Immune_All_Low.pkl', adata=adata, obs_cluster='leiden')
...
```

```python [Further Optimized]
# To further optimize the cell type annotation step, we can try
# using the `AnnotatorSCType` tool.
annotator = AnnotatorSCType()
adata = annotator.run(adata=adata, obs_cluster='leiden', path=cfg['output_dir'], tissue_type='Immune system')
...
```

:::
</div>
                <span class="bg-light-green-lighten-2">
                  After being evaluated by GPT-4, the labels for these categories were finally confirmed and saved as ob1.obs['final_type']:
                </span>
              </v-expansion-panel-text>
            </v-expansion-panel>
          </v-expansion-panels>
        
<div>

![An image](/final_annotation.png)

::: tip
CellAgent tried various solutions, evaluated their results, and ultimately produced a higher-quality result.
:::
</div>
        </v-card-text>
      </v-card>
    </v-col>
  </v-row>
</v-container>
