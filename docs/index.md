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
<v-responsive
  :aspect-ratio="16 / 9"
  class="border px-md-16 px-sm-0 py-0"
>
  <iframe
    src="//player.bilibili.com/player.html?isOutside=true&aid=112613522411165&bvid=BV1dVGoeCEQ4&cid=500001581492325&p=1"
    scrolling="no"
    allowfullscreen="true"
    class="h-100 w-100 border-0"
  ></iframe>
</v-responsive>

By constructing and coordinating several LLM-driven biological expert roles,
CellAgent conducts step-by-step execution and iterative optimization for various
tasks, substantially reducing the workload for science data analyses, bringing
us into the “Agent for Science” era.

<v-tabs
  bg-color="indigo-darken-2"
  fixed-tabs
>
  <v-tab text="Option"></v-tab>

  <v-tab text="Another Option"></v-tab>
</v-tabs>
[[TABS]]

### CellAgent consistently adapts appropriate tools and hyperparameters to achieve superior outcomes

<v-container class="my-16 py-0 px-8">
  <v-row>
    <v-col>
      <v-card
        title="Uniform Bar Exam"
        variant="plain"
      >
        <v-card-text class="py-0">
          <v-container class="pa-0">
            <v-row no-gutters>
              <v-col>
                <div style="font-size:84px">25th</div>
                <div>CellAgent</div>
              </v-col>
              <v-col class="text-light-green">
                <div style="font-size:84px">90th</div>
                <div>GPT-4</div>
              </v-col>
            </v-row>
          </v-container>
        </v-card-text>
      </v-card>
    </v-col>
    <v-col>
      <v-card
        title="Biology Olympiad"
        variant="plain"
      >
        <v-card-text class="py-0">
          <v-container class="pa-0">
            <v-row no-gutters>
              <v-col>
                <div style="font-size:84px">33rd</div>
                <div>CellAgent</div>
              </v-col>
              <v-col class="text-light-green">
                <div style="font-size:84px">99th</div>
                <div>GPT-4</div>
              </v-col>
            </v-row>
          </v-container>
        </v-card-text>
      </v-card>
    </v-col>
  </v-row>
</v-container>

### CellAgent surpasses GPT-4 and GPT-3.5 in its advanced automation capabilities on single-cell analysis.

<v-container class="my-8 py-0 px-8">
  <v-row>
    <v-col cols="6">
      <h3 class="mb-2">GPT-4.0</h3>
      <v-card
        border="surface sm opacity-100"
        title="Input"
        variant="text"
        class="rounded-0"
      >
        <v-card-text>
          This is a classic single-cell data analysis dataset, consisting of 3k human peripheral blood mononuclear cells from a healthy donor. Please perform the cell type annotation task and provide the complete Python code.
        </v-card-text>
      </v-card>
      <v-card
        border="surface sm opacity-100"
        title="Output"
        variant="text"
        class="rounded-0 border-t-0"
      >
        <v-card-text style="text-wrap-mode:wrap;white-space:pre-wrap;white-space-collapse:preserves">
          <span class="bg-grey-lighten-2">
Here is a list of their availability:<br>
- Andrew: 11 am to 3 pm<br>
- Joanne: noon to 2 pm, and 3:30 pm to 5 pm<br>
- Hannah: noon to 12:30 pm, and 4 pm to 6 pm<br><br>
Based on their availability, there is a 30-minute window where all three of them are available, which is from 4 pm to 4:30 pm. So, the meeting can be scheduled at 4 pm.
          </span>
        </v-card-text>
      </v-card>
    </v-col>
    <v-col>
      <h3 class="mb-2 text-success">CellAgent</h3>
      <v-card
        border="success sm opacity-100"
        title="Input"
        variant="text"
        class="rounded-0 text-success"
      >
        <v-card-text>
          Andrew is free from 11 am to 3 pm, Joanne is free from noon to 2 pm
          and then 3:30 pm to 5 pm. Hannah is available at noon for half an
          hour, and then 4 pm to 6 pm. What are some options for start times
          for a 30 minute meeting for Andrew, Hannah, and Joanne?
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
Andrew: 11 am - 3 pm<br>
Joanne: 12 pm - 2 pm, 3:30 pm - 5 pm<br>
Hannah: 12 pm - 12:30 pm, 4 pm - 6 pm<br><br>
Common availability for a 30-minute meeting: 12 pm - 12:30 pm
          </span>
        </v-card-text>
      </v-card>
    </v-col>
  </v-row>
</v-container>

CellAgent can streamline your single-cell data analysis workflow, ensuring
high-quality results with minimal effort. Our intuitive interface and robust
algorithms make it easy to process and interpret your data, regardless of your
level of expertise. With CellAgent website, you can:

<div class="mx-16 px-16">

* **Upload and Analyze:** Click to upload your single-cell data file, then chat
with CellAgent, providing necessary descriptions. The more detailed your input,
the better the results. CellAgent will finish executing through each step,
displaying the results as it goes.

* **Explore with Demos:** Click on our provided examples to quickly experience the impressive capabilities of CellAgent.
* **Interactive Requests:** Through ongoing dialogue, you can continuously submit new requests, then CellAgent will try to meet your needs at all times
</div>


## More on CellAgent

### Research

CellAgent is publicly accessible on BiorXiv. View CellAgent research.

### Meet the team

<div class="mx-16 px-16">

- **Prof.** [Jiajie Peng](https://github.com) Northwestern Polytechnical University
- **Prof.** [Jianye Hao](https://github.com) Tianjin University
</div>

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
  <v-btn rounded>Try on CellAgent 👉</v-btn>
  <v-btn variant="plain">View CellAgent research ></v-btn>
</v-sheet>
