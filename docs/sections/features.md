## By coordinating several LLM-driven biological experts, CellAgent automatically conducts step-by-step execution and iterative optimization for various tasks, substantially reducing the workload for science data analyses, bringing us into the “Agent for Science” era.

<div class="mb-4"></div>

<v-tabs class="border-b-sm" fixed-tabs v-model="tab">
  <v-tab text="Task decomposition" value="0"></v-tab>
  <v-tab text="Execution for sub-step" value="1"></v-tab>
  <v-tab text="Self-optimization" value="2"></v-tab>
  <v-tab text="Handle user feedback" value="3"></v-tab>
</v-tabs>

<v-tabs-window class="mt-8" v-model="tab">
  <v-tabs-window-item value="0">
    <demo-frame no-thought>
      <template v-slot:description>
        For complex scRNA-seq data analysis tasks, CellAgent can emulate the process of a human expert by first breaking down the task into sub-steps and then executing them sequentially.
      </template>
      <template v-slot:input>
        This is a classic single-cell data analysis dataset, consisting of 3k human peripheral blood mononuclear cells from a healthy donor. Please perform the cell type annotation task.
      </template>
      <template v-slot:output>
After an initial analysis, 6 steps are required to complete your request:
<div>

```
1. Quality Control
2. Normalization
3. Identification of Highly Variable Genes
4. Dimensionality Reduction
5. Clustering
6. Cell Type Annotation
```
</div>
      </template>
    </demo-frame>
  </v-tabs-window-item>

  <v-tabs-window-item value="1">
    <demo-frame no-input>
      <template v-slot:description>
        During each step, CellAgent can automatically generate and execute code, providing the resulting plots and data.
      </template>
      <template v-slot:input>
        your input
      </template>
      <template v-slot:thought>
<div>
Data Format: 

```python
AnnData object with n_obs × n_vars = 2638 × 32738
    obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt'
    var: 'gene_ids', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'mt', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'log1p', 'hvg'
    layers: 'counts'
```
The code has been ran: `[history code]`
The current step is: Dimensionality Reduction. Generate code to complete the current step.

</div>
      </template>
      <template v-slot:output>
<div>
To perform dimensionality reduction, we will use Principal Component Analysis (PCA) to reduce the dimensionality of the dataset. This step is crucial for visualizing the data and for subsequent steps like clustering and trajectory inference. Here is the code to perform PCA and visualize the results

```python
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True, show=True, save=False)
plt.savefig(f"{cfg['img_save_dir']}/pca_variance_ratio.png")

sc.pl.pca(adata, color='total_counts', show=True, save=False)

```

<img src="/image_2.png" alt="" style="height: 200px;">
</div>
      </template>
    </demo-frame>
  </v-tabs-window-item>

  <v-tabs-window-item value="2">
    <demo-frame no-input>
      <template v-slot:description>
        CellAgent can automatically optimize the generated solutions and ultimately select the best one to present.
      </template>
      <template v-slot:thought>
      The current result doesn't look good enough yet. You can optimize this step by adjusting parameters, changing tools, and so on.
      </template>
      <template v-slot:output>
<div>
Previously, CellAgent used the `AnnotatorCellmarkerACT` tool to annotate cell types. To optimize the cell type annotation step, we can try using a different tool or adjusting the parameters. Given the options available, let's use the `AnnotatorCelltypist` tool, which is known to work well if it can find a good matching model. We will use the "Immune_All_Low.pkl" model, which is suitable for immune sub-populations.

Here is the complete code to perform cell type annotation using the `AnnotatorCelltypist` tool:

```python
from tools.annotator_celltypist import AnnotatorCelltypist

annotator = AnnotatorCelltypist()

adata = annotator.run(model_name='Immune_All_Low.pkl', adata=adata, obs_cluster='leiden')

print("Cell type annotation complete. The annotated cell types are stored in adata.obs['cell_typist_label'].")
```

This code initializes the `AnnotatorCelltypist` tool, runs the annotation process using the "Immune_All_Low.pkl" model, and stores the annotated cell types in the `adata.obs['cell_typist_label']` column. This model is specifically designed for immune sub-populations, which should be suitable for the given dataset of human peripheral blood mononuclear cells.

![An image](/output.png)

</div>
      </template>
    </demo-frame>
  </v-tabs-window-item>

  <v-tabs-window-item value="3">
    <demo-frame>
      <template v-slot:description>
        CellAgent can understand your natural language requirements and, at any point during the execution, respond to your natural language feedback by generating new code to meet your needs.
      </template>
      <template v-slot:input>
        Set the npca parameter to 30 and display the Louvain clustering results.
      </template>
      <template v-slot:thought>
<div>
The current step is: Dimensionality Reduction.

In this step, CellAgent executed the code `sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)` to compute the neighborhood graph.

The current step has been completed, and the user has provided new feedback: Set the npca parameter to 30 and display the Louvain clustering results.

Please continue to generate code and natural language explanations to better meet the user's needs.
</div>
      </template>
      <template v-slot:output>
<div>
To meet the user's requirements, we will set the `npca` parameter to 30 for the PCA and then perform Louvain clustering. We will visualize the clustering results using UMAP. Here is the code to complete the current step:

```python
...

# Perform PCA with npca=30
sc.tl.pca(adata, svd_solver='arpack', n_comps=30)

sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

sc.tl.umap(adata)

# Perform Louvain clustering
sc.tl.louvain(adata)

# Visualize the UMAP with Louvain clustering results
sc.pl.umap(adata, color=['louvain'])
...
```

Explanation:
1. **PCA with npca=30**: We perform PCA with the number of principal components set to 30.
2. **Louvain Clustering**: We perform Louvain clustering to identify clusters in the data.
...


<img src="/output1.png" alt="" style="height: 200px;">
<img src="/output2.png" alt="" style="height: 200px;">

</div>
      </template>
    </demo-frame>
  </v-tabs-window-item>
</v-tabs-window>
