---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "CellAgent"
  text: "An LLM-driven agent for single-cell data analysis, ensuring high-quality results with minimal effort."
  # tagline: My great project tagline
  actions:
    - theme: brand
      text: Try on CellAgent 👉
      link: http://cell.agent4science.cn/
    - theme: alt
      text: View CellAgent research >
      link: https://www.biorxiv.org/content/10.1101/2024.05.13.593861v1
---

<script setup>
import { ref, onMounted } from 'vue'

const tab = ref('0')
const isCN = ref(false)

onMounted(async function assertInCN() {
  try {
    const response = await fetch('//ipinfo.io/json');
    const address = await response.json();

    isCN.value = address.country === 'CN';
  } catch {
    isCN.value = false;
  }

  console.log(isCN.value)
})
</script>

<v-responsive
  :aspect-ratio="16 / 9"
  class="border-0 px-md-16 px-sm-0 py-0"
>
  <iframe
    v-if="isCN"
    src="//player.bilibili.com/player.html?isOutside=true&aid=112613522411165&bvid=BV1dVGoeCEQ4&cid=500001581492325&p=1"
    scrolling="no"
    allowfullscreen="true"
    class="h-100 w-100 border-0"
  ></iframe>
  <iframe
    v-else
    class="h-100 w-100 border-0"
    src="https://www.youtube.com/embed/7a4M3ymp5ng?si=Hp-jAv9KkYHy-4-w"
    title="YouTube video player"
    allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share"
    referrerpolicy="strict-origin-when-cross-origin"
    allowfullscreen
    ></iframe>
</v-responsive>

<br><br>By constructing and coordinating several LLM-driven biological expert roles,
CellAgent conducts step-by-step execution and iterative optimization for various
tasks, substantially reducing the workload for science data analyses, bringing
us into the “Agent for Science” era.<br><br>

<v-tabs fixed-tabs v-model="tab">
  <v-tab text="Task decomposition" value="0"></v-tab>
  <v-tab text="Execution for sub-step" value="1"></v-tab>
  <v-tab text="Self-optimization" value="2"></v-tab>
  <v-tab text="Handle user feedback" value="3"></v-tab>
</v-tabs>

<v-tabs-window v-model="tab">
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
Data Format: AnnData object with n_obs × n_vars = 2638 × 32738
    obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt'
    var: 'gene_ids', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'mt', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'log1p', 'hvg'
    layers: 'counts'

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


### CellAgent surpasses traditional data analysing approaches in its advanced automation capabilities on single-cell analysis.

<v-container class="my-8 py-0 px-8">
  <v-row>
    <v-col cols="6">
      <h3 class="mb-2">Traditional Method</h3>
      <v-card
        border="surface-variant sm opacity-100"
        title="Input"
        variant="text"
        class="rounded-0"
      >
        <v-card-text>
          Data file: xxx.h5ad<br>
          Task Requirement: Annotate cell types
        </v-card-text>
      </v-card>
      <v-card
        border="surface-variant sm opacity-100"
        title="Output"
        variant="text"
        class="rounded-0 border-t-0"
      >
        <v-card-text style="text-wrap-mode:wrap;white-space:pre-wrap;white-space-collapse:preserves">
          <span class="bg-grey-lighten-2">Sure, I can help with that. We'll use the Scanpy library in Python for single-cell RNA-seq data analysis. Here's the general workflow we'll follow:<br>
1. Load the dataset.<br>
2. Preprocess the data (filtering, normalization, etc.).<br>
3. Perform clustering.<br>
4. Annotate the cell types.<br>
Let's start by providing the complete Python code for this task:<br>
</span>
<div>

```python
...
cell_types = pd.DataFrame(index=adata.obs_names)
for cell_type, genes in marker_genes.items():
    adata.obs[cell_type] = adata[:, genes].X.mean(1)
...
```
</div>

<div>

::: danger
GPT-4 provided code for a basic solution.
:::
</div>
        </v-card-text>
      </v-card>
    </v-col>
    <v-col cols="6">
      <h3 class="mb-2 text-success">CellAgent</h3>
      <v-card
        border="success sm opacity-100"
        title="Input"
        variant="text"
        class="rounded-0 text-success"
      >
        <v-card-text>
          This is a classic single-cell data analysis dataset, consisting of 3k human peripheral blood mononuclear cells from a healthy donor. Please perform the cell type annotation task.
        </v-card-text>
      </v-card>
      <v-card
        border="success sm opacity-100"
        title="Output"
        variant="text"
        class="rounded-0 border-t-0 text-success"
      >
        <v-card-text>
        <span class="bg-light-green-lighten-2">
          CellAgent is starting the process. Please wait for a while.<br>
          ... [the process of task decomposition and the execution of previous steps]<br>
        </span>
<div>

```python
# Cell Type Annotation using AnnotatorCellmarkerACT
annotator = AnnotatorCellmarkerACT()
adata = annotator.run(species='Human', tissue_type='Blood', adata=adata, obs_cluster='leiden')
...
```

```python
# To optimize the cell type annotation step, let's use the
# `AnnotatorCelltypist` tool and the "Immune_All_Low.pkl" model,
# which is suitable for immune sub-populations.
annotator = AnnotatorCelltypist()
adata = annotator.run(model_name='Immune_All_Low.pkl', adata=adata, obs_cluster='leiden')
...
```

```python
# To further optimize the cell type annotation step, we can try
# using the `AnnotatorSCType` tool.
annotator = AnnotatorSCType()
adata = annotator.run(adata=adata, obs_cluster='leiden', path=cfg['output_dir'], tissue_type='Immune system')
...
```

<span class="bg-light-green-lighten-2">
After being evaluated by GPT-4, the labels for these categories were finally confirmed and saved as ob1.obs['final_type']:
</span>

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

### CellAgent consistently adapts appropriate tools and hyperparameters to achieve superior outcomes.

<v-container class="my-16 py-0 px-8">
  <v-row>
    <v-col>
      <v-card
        title="Task Completion Rate"
        variant="plain"
      >
        <v-card-text class="py-0">
          <v-container class="pa-0">
            <v-row no-gutters>
              <v-col>
                <div style="font-size:84px">47%</div>
                <div>GPT-4</div>
              </v-col>
              <v-col class="text-light-green">
                <div style="font-size:84px">92%</div>
                <div>CellAgent</div>
              </v-col>
            </v-row>
          </v-container>
        </v-card-text>
      </v-card>
    </v-col>
    <v-col>
      <v-card
        title="Task Performance**"
        variant="plain"
      >
        <v-card-text class="py-0">
          <v-container class="pa-0">
            <v-row no-gutters>
              <!-- <v-col>
                <div style="font-size:84px"></div>
                <div>CellAgent</div>
              </v-col> -->
              <v-col class="text-light-green">
                <div style="font-size:84px">107.23%</div>
                <div>CellAgent*</div>
              </v-col>
            </v-row>
          </v-container>
        </v-card-text>
      </v-card>
    </v-col>
  </v-row>
</v-container>

\* In typical scRNA-seq data analysis tasks, CellAgent's performance can reach 107.23% compared to the widely used and effective existing algorithms.

** The tasks referred to here mainly include batch effect correction, cell type annotation, and trajectory inference, corresponding to the existing algorithms Scanorama, GPT-4 annotation, and Slingshot, respectively.<br><br><br>


CellAgent can streamline your single-cell data analysis workflow, ensuring
high-quality results with minimal effort. Our intuitive interface and robust
algorithms make it easy to process and interpret your data, regardless of your
level of expertise. With CellAgent website, you can:

<div class="mx-8 px-8">

* **Upload and Analyze:** Click to upload your single-cell data file, then chat
with CellAgent, providing necessary descriptions. The more detailed your input,
the better the results. CellAgent will finish executing through each step,
displaying the results as it goes.

* **Explore with Demos:** Click on our provided examples to quickly experience the impressive capabilities of CellAgent.
* **Interactive Requests:** Through ongoing dialogue, you can continuously submit new requests, then CellAgent will try to meet your needs at all times
</div>


## More on CellAgent

### Research

CellAgent is publicly accessible on BiorXiv. 

<!-- ### Meet the team

<div class="mx-16 px-16">

- **Prof.** [Jiajie Peng](https://github.com) Northwestern Polytechnical University
- **Prof.** [Jianye Hao](https://github.com) Tianjin University
</div> -->

### Contact us

If you have any suggestions or concerns during use, please feel free to contact
the developer (email: TBD) or the corresponding authors (Professor Peng and
Professor Hao). In your email, please specify the time the issue occurred and
include your usage record of CellAgent (screenshots, etc.) to help us identify
the problem. Thank you.

<v-sheet class="mt-8 d-flex align-center justify-center flex-wrap text-center mx-auto pa-16 bg-grey-darken-4" elevation="4" max-width="800" width="100%">
  <div class="text-h5 font-weight-medium mb-8">
We are excited to see the potential of CellAgent to greatly enhance productivity,
foster new discoveries, and deepen our understanding of biological systems.
  </div>
  <v-btn rounded href="http://cell.agent4science.cn/">Try on CellAgent 👉</v-btn>
  <v-btn variant="text" href="https://www.biorxiv.org/content/10.1101/2024.05.13.593861v1" target="_blank">View CellAgent research ></v-btn>
</v-sheet>
