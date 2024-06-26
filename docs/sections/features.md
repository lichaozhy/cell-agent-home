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
# Perform PCA
sc.tl.pca(adata, svd_solver='arpack')

# Plot the explained variance ratio to understand how much variance is captured by each principal component
sc.pl.pca_variance_ratio(adata, log=True, show=True, save=False)
plt.savefig(f"{cfg['img_save_dir']}/pca_variance_ratio.png")

# Plot the PCA result
sc.pl.pca(adata, color='total_counts', show=True, save=False)
plt.savefig(f"{cfg['img_save_dir']}/pca_result.png")

# Print description of the figures
print("The first figure shows the explained variance ratio of the principal components. The second figure shows the PCA result, with cells colored by their total counts.")
```
![](/image_2.png)
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
# Cell Type Annotation using AnnotatorCelltypist
from tools.annotator_celltypist import AnnotatorCelltypist

# Initialize the AnnotatorCelltypist tool
annotator = AnnotatorCelltypist()

# Run the annotation using the "Immune_All_Low.pkl" model
adata = annotator.run(model_name='Immune_All_Low.pkl', adata=adata, obs_cluster='leiden')

# Print a message indicating that cell type annotation is complete
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
# Dimensionality reduction using PCA with npca=30 and Louvain clustering
import scanpy as sc
import matplotlib.pyplot as plt

# Perform PCA with npca=30
sc.tl.pca(adata, svd_solver='arpack', n_comps=30)

# Visualize the explained variance ratio to determine the number of principal components to use
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig(f"{cfg['img_save_dir']}/pca_variance_ratio_30.png")
print("PCA variance ratio plot with npca=30 saved as 'pca_variance_ratio_30.png'.")

# Compute the neighborhood graph using 30 principal components
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

# Perform UMAP
sc.tl.umap(adata)

# Perform Louvain clustering
sc.tl.louvain(adata)

# Visualize the UMAP with Louvain clustering results
sc.pl.umap(adata, color=['louvain'])
plt.savefig(f"{cfg['img_save_dir']}/umap_louvain.png")
print("UMAP plot with Louvain clustering results saved as 'umap_louvain.png'.")

# Save the data with PCA, UMAP, and Louvain clustering results
adata.write(f"{cfg['output_dir']}/dimensionality_reduction_louvain_data.h5ad")
print("Dimensionality reduction with Louvain clustering complete and data saved as 'dimensionality_reduction_louvain_data.h5ad'.")
```

Explanation:
1. **PCA with npca=30**: We perform PCA with the number of principal components set to 30.
2. **Visualization**: We visualize the explained variance ratio for the PCA with 30 components.
3. **Neighborhood Graph**: We compute the neighborhood graph using 30 principal components.
4. **UMAP**: We perform UMAP for further dimensionality reduction and visualization.
5. **Louvain Clustering**: We perform Louvain clustering to identify clusters in the data.
6. **Visualization**: We visualize the UMAP results with the Louvain clustering results.
7. **Save Data**: We save the AnnData object with the PCA, UMAP, and Louvain clustering results.

![An image](/output1.png)
![An image](/output2.png)
</div>
      </template>
    </demo-frame>
  </v-tabs-window-item>
</v-tabs-window>
