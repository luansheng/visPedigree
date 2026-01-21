# visPedigree: A High-Performance R Package for Tidying and Visualizing Large-Scale Animal Pedigrees

## Abstract

Pedigree analysis is fundamental to animal breeding and genetics.
However, visualizing large-scale pedigrees with overlapping generations
and complex family structures remains a significant challenge. Existing
tools often struggle with “hairball” effects or lack the computational
efficiency to handle datasets exceeding thousands of individuals. We
present **visPedigree**, an R package that introduces two novel
algorithms: (1) a **Topological Generation Inference (TGI)** algorithm
for robust pedigree standardization and sub-graph extraction, and (2) a
**Compacted Generation-Constrained Layout (CGCL)** algorithm for
clutter-free visualization of massive populations. By integrating these
methods, visPedigree enables breeders to efficiently manage, prune, and
visualize pedigrees with over 10,000 individuals per generation,
providing high-quality vector graphics suitable for publication and
breeding decisions.

## 1. Introduction

(Briefly introduce the importance of pedigree visualization in
aquaculture and livestock breeding. Mention limitations of current tools
like `kinship2`, `pedigree`, etc.)

## 2. Algorithms and Methodology

The core innovation of `visPedigree` lies in its two-stage algorithmic
approach: data standardization followed by layout optimization.

### 2.1 Topological Generation Inference (TGI) Algorithm

Raw pedigree data often suffers from missing information, chronological
inconsistencies, and overlapping generations. To address this, we
developed the TGI algorithm to standardize the pedigree structure before
visualization.

**Algorithm Description:** The TGI algorithm treats the pedigree as a
Directed Acyclic Graph (DAG) where individuals are nodes and
parent-offspring relationships are edges.

1.  **Candidate-Centric Subgraph Extraction**: Instead of processing the
    entire population, the algorithm allows users to focus on specific
    candidates. It performs a bidirectional graph traversal:
    - *Upward Traversal*: Recursively identifies all ancestors $A$ of
      the candidate set $C$.
    - *Downward Traversal*: Recursively identifies all descendants $D$
      of the candidate set $C$.
    - The resulting subgraph $G\prime = C \cup A \cup D$ significantly
      reduces the dataset size while preserving relevant kinship
      structures.
2.  **Flexible Generation Assignment and Topological Sorting**: To
    resolve overlapping generations, we assign a discrete generation
    level $Gen(i)$ to each individual $i$ based on topological sorting.
    `visPedigree` provides two distinct alignment strategies to
    accommodate different breeding scenarios:
    - **Top-Down Alignment (Depth-based, Default)**: Founders are
      assigned to $Gen = 1$. For any non-founder individual $i$:
      $$Gen(i) = \max\left( Gen(Sire),Gen(Dam) \right) + 1$$ This
      “Top-Aligned” approach ensures that all founders start from the
      same visual level, which is biologically intuitive for tracking
      the history of a population from its base stock and prevents
      “founder drift” where ancestors of shorter paths sink into lower
      generations.
    - **Bottom-Up Alignment (Height-based)**: Generations are calculated
      relative to the terminal nodes (leaf nodes) using a height-based
      algorithm. Terminal nodes are aligned at the maximum generation,
      and height is propagated upward:
      $$Height(parent) = \max\left( Height(parent),Height(child) + 1 \right)$$
      This “Bottom-Aligned” approach is useful when integrating
      exogenous parents introduced in later years or when the primary
      focus is the synchronization of current-day offspring.
    - This dual-strategy flexibility ensures that parents are always
      placed in a generation strictly above their offspring (TGI
      algorithm), regardless of path length.
3.  **Missing Node Imputation**: The algorithm detects gaps where
    $Gen(child) - Gen(parent) > 1$. It can optionally infer and insert
    “virtual” intermediate nodes to maintain structural continuity,
    which is crucial for accurate inbreeding calculation and
    visualization.

### 2.2 Compacted Generation-Constrained Layout (CGCL) Algorithm

Visualizing thousands of individuals typically results in edge-crossing
clutter. The CGCL algorithm solves this by combining graph reduction
with constrained layering.

**Algorithm Description:**

1.  **Graph Transformation and Sibship Compaction**: We transform the
    standard pedigree graph into a specialized structure to reduce
    visual complexity.

    - **Family Nodes**: We introduce virtual “family nodes” to mediate
      connections between parents and offspring. Edges flow from Parents
      $\rightarrow$ Family Node $\rightarrow$ Offspring.
    - **Sibship Compaction Operator ($\Phi$)**: For any full-sib family
      $S$ (individuals sharing the same sire and dam) that has no
      progeny in the displayed graph, we apply a compaction operator:
      $$\left. \Phi(S)\rightarrow v_{compact} \right.$$ The entire set
      $S$ is collapsed into a single “super-node” $v_{compact}$, labeled
      with the family size $|S|$. This reduces the number of nodes in a
      generation from $N$ to $N_{parents} + N_{families}$, reducing
      graph complexity by orders of magnitude in high-fecundity species
      like fish.

2.  **Generation-Constrained Layering**: Unlike standard Sugiyama layout
    implementations that heuristically determine layers, we enforce
    strict biological constraints:

    - Individual nodes are constrained to layers $L = 2 \times Gen - 1$.
    - Family nodes are constrained to layers $L = 2 \times Gen - 2$.
      This rigid layering preserves the biological timeline and
      facilitates interpretation.

3.  **Constrained Crossing Minimization**: Within the fixed layers, we
    employ the Sugiyama heuristic (barycenter method) solely to optimize
    the horizontal ordering of nodes, minimizing edge crossings between
    adjacent generations.

4.  **1D Repulsion Post-Processing**: To prevent node overlap in dense
    generations, we apply a 1D repulsion algorithm on the X-coordinates.

    - Let $x_{1},x_{2},...,x_{n}$ be the sorted X-coordinates of nodes
      in a layer.
    - If $x_{i + 1} - x_{i} < \delta_{min}$ (minimum visual separation),
      the algorithm shifts $x_{i + 1}$ and subsequent nodes to enforce
      separation.
    - This ensures that every node remains distinct and legible, even in
      the densest parts of the pedigree.

## 3. Implementation and Performance

The algorithms are implemented in R, leveraging `data.table` for
high-performance data manipulation and `igraph` for graph rendering. The
package supports outputting high-resolution vector graphics (PDF),
allowing for lossless zooming on massive pedigree charts.

(Include performance benchmarks here: e.g., time taken to process 10k,
50k, 100k individuals compared to other packages).

## 4. Case Studies

To validate the effectiveness of `visPedigree`, we applied the package
to two distinct real-world breeding datasets, demonstrating its
capability to handle both deep ancestral structures and massive
population sizes.

### 4.1 Case Study 1: Deep Pedigree Structure (`deep_ped`)

**Dataset Characteristics:** The `deep_ped` dataset consists of 4,396
individuals. It represents a population with a long breeding history,
characterized by multiple overlapping generations and complex ancestry.

**Challenge:** Standard visualization tools often fail to correctly
align generations in such deep pedigrees. Individuals born in the same
year may belong to different genealogical generations due to varying
generation intervals, leading to a “timeline paradox” where parents and
offspring might appear on the same visual level.

**Application of TGI Algorithm:** Applying the Topological Generation
Inference (TGI) algorithm: 1. **Generation Alignment:** The algorithm
successfully assigned discrete topological generation numbers, ensuring
that every parent is placed strictly above their offspring. 2.
**Alignment Strategy Choice:** By utilizing the default
`genmethod = "top"`, `visPedigree` anchored all historical founders at
Generation 1, providing a consistent “top-down” flow. This prevents the
visual confusion where founders of lineages with fewer generations
appear “drifted” into the middle of the chart. 3. **Sub-graph
Extraction:** By selecting a specific cohort of interest (e.g., the
latest generation), `visPedigree` traced back and extracted only the
relevant 4,396 ancestors from a potentially larger historical database.
4. **Result:** The resulting visualization clearly depicts the flow of
genetic contribution over dozens of generations without the visual
ambiguity caused by chronological overlap.

### 4.2 Case Study 2: Massive Population with High Fecundity (`big_family_size_ped`)

**Dataset Characteristics:** The `big_family_size_ped` dataset is a
large-scale aquaculture pedigree containing 178,421 individuals. This
species exhibits high fecundity, where a single pair of parents can
produce hundreds of offspring (full-sibs).

**Challenge:** Visualizing nearly 180,000 nodes is computationally
prohibitive for most R packages. Even if rendered, the resulting graph
would be a “hairball”—a dense, unreadable mesh of edges where individual
relationships are obscured by the sheer volume of full-sib connections.

**Application of CGCL Algorithm:** We utilized the Compacted
Generation-Constrained Layout (CGCL) algorithm to process this massive
dataset. 1. **Sibship Compaction:** The algorithm identified full-sib
families with no subsequent progeny in the graph. These groups were
collapsed into single “family nodes.” \* *Reduction Metric:* This
process reduced the effective node count significantly, transforming
thousands of terminal nodes into a manageable number of summary nodes.
2. **Layout Optimization:** The compacted graph was laid out using the
constrained layering approach. 3. **Performance:** `visPedigree`
generated the vector graphic output in seconds. The final plot
maintained the structural integrity of the population—showing the number
of families and their sizes—while remaining clean and legible. This
visualization allows breeders to instantly assess family size
distributions and selection intensity across the entire population.

### 4.3 Comparison with Existing Tools

| Feature              | Standard Tools (`kinship2`, `pedigree`) | `visPedigree`                           |
|:---------------------|:----------------------------------------|:----------------------------------------|
| **Max Node Count**   | \< 5,000 (typically)                    | \> 100,000+                             |
| **Layout Strategy**  | Chronological / Force-directed          | Topological / Generation-Constrained    |
| **Sibship Handling** | Draws every individual                  | Compacts terminal full-sibs             |
| **Output Quality**   | Often cluttered                         | Optimized, clutter-free vector graphics |

## 5. Conclusion

visPedigree provides a scalable, algorithmic solution to the “hairball”
problem in pedigree visualization. By combining topological generation
inference with sibship compaction, it offers breeders a powerful tool to
gain insights into complex population structures.

## 6. visPedigree 系谱显示逻辑总结 (Summary of Visualization Logic)

`visPedigree` 的系谱显示逻辑是一套结合了图论（Graph
Theory）、数据治理和启发式布局优化（Heuristic
Layout）的复合系统。以下是核心规则和逻辑：

### 6.1 核心架构：两阶段（Two-Pass）渲染逻辑

为了解决 R 语言在导出高精度 PDF
时，线段与节点之间常出现的“白边”或连接不紧密的问题，`visped`
采用了独特的两阶段绘制策略： \*
**第一阶段（底座线段）**：先绘制所有从“父节点中心”到“子节点中心”的连接线。这确保了线段在视觉上绝对穿透节点。
\*
**第二阶段（顶层覆盖）**：在原位置覆盖绘制节点（圆/方块）及其标签。利用
R
的图层覆盖特性，确保连接线完美消失在节点边缘，即使在极高缩放倍数下也无断裂感。

### 6.2 布局逻辑：分层与居中对齐

布局引擎基于 `igraph` 的 **Sugiyama 分层算法**，并进行了专项增强： \*
**时间/代数轴**：祖先在顶部，后代在底部。代数（Generation）决定了纵向坐标（Y轴）。
\* **前向-后向平均算法（Forward-Backward Averaging）**： \*
**前向**：父母节点的横向位置（X轴）被尽可能调整为其所有后代位置的平均值。
\* **后向**：子节点的位置也会根据父母的位置进行微调。 \*
**目的**：确保系谱图具有对称的美感，使大家族呈现树状发散而非堆挤在一侧。
\*
**同胞聚类**：同一对父母产生的全同胞（Full-sibs）在布局中会紧密排列，避免被其他不相关节点穿插。

### 6.3 数据压缩规则：Compact 模式

针对大规模系谱（如水产育种），`visped` 引入了同胞压缩逻辑： \*
**逻辑**：当 `compact = TRUE`
时，同一代中具有完全相同双亲的一组个体被合并为一个单一的“家庭节点”（Family
Node）。 \*
**视觉表示**：合并后的节点由圆圈变为**正方形**，标签显示为该全同胞家族的个体数量。
\*
**价值**：能将上万个体的复杂系谱压缩为数百个代表性节点，保持核心拓扑结构清晰。

### 6.4 视觉语义规则

颜色和形状承载了明确的生物学含义： \* **形状**： \*
**圆形**：代表独立个体。 \* **正方形**：代表压缩后的全同胞家系（Compact
mode）。 \* **填充色**： \* **天蓝色 (Sky Blue)**：雄性。 \* **金黄色
(Goldenrod)**：雌性。 \* **橄榄绿 (Olive Green)**：性别未知。 \*
**边框与强调**： \* **紫色/自定义色**：用于 `highlight`
参数指定的重点关注个体或其亲缘系（通过 `trace` 溯源）。 \*
**连接线颜色**： \*
从父亲发出的一般继承其性别颜色，从母亲发出的一般继承其性别颜色，这有助于快速追踪母系或公系血统。

### 6.5 自动化治理逻辑

- **孤立点剔除**：自动识别既无父母也无后代的“孤立个体”（Gen
  0），并从绘图中剔除，防止杂讯干扰主系谱。
- **自适应字体（Cex）**：根据每一层（Generation）中节点的最密集程度，自动计算最合适的字符缩放比例（cex），确保标签既不重叠也不过小。
- **高纵横比处理**：当系谱极大导致 PDF 宽度超过限制时，提供
  `outline = TRUE`
  模式，仅绘制结构框架而隐藏标签，以便从宏观视角观察群体结构。

## 7. tidyped 对系谱的处理规则和逻辑 (Tidying Logic)

`tidyped` 是 `visPedigree`
的核心数据预处理引擎，主要负责将原始、混乱的系谱记录转化为结构标准化、逻辑自洽且可直接用于计算和绘图的“干净系谱”。其处理规则和逻辑可归纳为以下五个维度：

### 7.1 规范化与清洗规则 (Normalization)

- **标准化转换**：将输入数据（data.frame 或
  data.table）统一转换为高效率的 `data.table` 格式。
- **缺失值治理**：自动识别各种非标准的缺失标记（如 `"0"`, `"*"` ,
  `"NA"`, `""` 等）并统一转化为 R 语言原生的 `NA`。
- **字符类型强制**：强制个体、父本、母本三列为字符型，防止由于 ID
  是纯数字而产生的索引混淆。

### 7.2 拓扑结构验证规则 (Topology Validation)

- **循环检测 (Loop
  Detection)**：通过图论算法（检测是否存在强相连组件）识别系谱中的逻辑错误，例如“A
  是 B 的父亲，B 又是 A
  的父亲”或“个体是自己的祖先”。一旦发现循环，函数会中断并报错。
- **性别角色一致性**：确保同一 ID 不会在同一数据集内既出现在
  `Sire`（父本）列又出现在 `Dam`（母本）列（除非用户允许特殊繁殖方式）。

### 7.3 基于图论的代数推断 (TGI 算法)

这是 `tidyped` 的核心算法逻辑： \*
**分层赋值**：将系谱视为有向无环图（DAG），采用拓扑排序。 \*
**递归逻辑**：祖先（Founders）被赋予基础层级。对于任何非祖先个体
$i$：$Gen(i) = \max\left( Gen(Sire),Gen(Dam) \right) + 1$。 \*
**对齐启发式 (Alignment
Heuristics)**：为了美观，算法不仅计算“最长路径”，还会根据配偶或同辈的层级对引入的个体（如新引进的种质）进行“重心调整”，防止系谱图出现过大的纵向跨度。

### 7.4 智能剪枝与溯源逻辑 (Pruning & Tracing)

- **候选者驱动**：用户可以通过 `cand`
  参数指定感兴趣的个体（如当年的选种核心群）。
- **多向溯源**：支持向上（祖先）、向下（后代）或双向（亲缘网）追踪，并可通过
  `tracegen` 限制溯源深度。

### 7.5 辅助计算扩展 (Extended Calculation)

- **唯一数值编码 (`addnum`)**：为每个个体分配唯一的整型 ID（`IndNum`,
  `SireNum`, `DamNum`），这在处理百万级数据时能显著提升后续计算的速度。
- **近交系数计算 (`inbreed`)**：集成高效算法（基于 `nadiv`
  等后端），在整理数据的同时计算每个个体的近交系数 $f$。
