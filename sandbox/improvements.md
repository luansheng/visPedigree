# visPedigree 优化与验证计划

> 创建日期：2026-03-08
> 分支：`pedigreeAnalysis`

---

## 一、文档与示例修复（高优先级）

### 1.1 `pedstats`

- [X] 补充缺失的 `@param` 文档：`calc_ecg`、`calc_genint` 两个参数在函数签名中存在，但 roxygen 中完全未文档化
- [X] 完善 `@return`：在 `$summary`、`$ecg`、`$gen_intervals` 条目中明确列出所有列名
  - `$summary`：`N`, `N_Sire`, `N_Dam`, `N_Founder`, `Max_Gen`
  - `$ecg`：`Ind`, `ECG`, `FullGen`, `MaxGen`
  - `$gen_intervals`：`Pathway`, `N`, `Mean`, `SD`
- [X] 将 `@examples` 中的 `\dontrun{}` 改为 `\donttest{}`，并改用包内数据集（`simple_ped`、`big_family_size_ped`）实现可运行示例
- [X] 验证：带 `timevar` 和不带 `timevar` 两种场景下的完整运行测试

### 1.2 `pedinbreed_class`

- [X] 补全 `@return`：明确列出返回列名 `F_Class`（有序 factor）、`Count`（integer）、`Percentage`（numeric）及 5 个固定等级
- [X] 新增 `@examples`：使用 `inbred_ped` 或 `simple_ped` 包内数据集，至少一个可直接运行的 `\donttest{}` 示例
- [X] 验证：各近交等级阈值（0, 0.0625, 0.125, 0.25）的分类逻辑是否正确

### 1.3 `pedrel`

- [X] **修复文档 Bug**：`@return` 中声称返回列名为 `Group`，但代码实际上执行了 `setnames(result, "Group", by)`，真实列名是 `by` 参数的值（如 `Gen`、`Year`）
- [X] 新增 `@examples`：覆盖 `by = "Gen"` 和 `by = "Year"` 两种用法，以及 `reference` 参数的使用
- [X] 验证：`reference` 参数在不同 `by` 列下的过滤行为

### 1.5 补充缺失 `@examples` 的函数

以下函数目前完全没有 `@examples`，需新增：

- [X] `pedne`：使用 `simple_ped` / `inbred_ped`，覆盖 `method = "coancestry"`（默认）、`"inbreeding"`、`"demographic"` 三种
- [X] `pedecg`：使用 `simple_ped`，展示 ECG / 系谱完整度计算
- [X] `pedsubpop`：使用包内数据集，展示亚群体统计输出

---

## 二、代码验证（中优先级）

- [x] `pedstats`：在 `big_family_size_ped` 上执行带 `timevar = "Year"` 的完整测试，检查 `$gen_intervals` 输出
- [x] `pedrel`：大型系谱（>25,000 个体）下 `compact = TRUE` 与 `compact = FALSE` 的结果一致性检验
- [x] `pedinbreed_class`：构造包含所有 5 个等级个体的人工系谱，验证分类统计正确
- [x] `pedancestry`：多品系混合系谱下各标签比例之和是否等于 1
- [x] `pediv` → `pedcontrib` → `pedne` 调用链：确保 `reference` 参数在所有三层函数中传递一致

---

## 三、Vignette 覆盖补充（中优先级）

`vignettes/pedigree-analysis.Rmd` 目前缺少以下函数的独立示例章节：

- [x] `pedrel`：新增"平均亲缘系数趋势"章节
- [ ] `pedecg`：新增"系谱完整度（ECG）"章节
- [ ] `pedgenint`：新增"世代间隔"章节（目前仅通过 `pedstats` 间接涉及）
- [ ] `pedsubpop`：新增"亚群体分析"章节

---

## 四、育种群体扩展功能提案（低优先级 / 讨论中）

以下功能在现有包中**完全缺失**，适合育种实践场景，供评估是否纳入：

### 4.1 选择强度（Selection Intensity）

- 基于截断选择理论，计算不同选择比例下的理论选择强度 $i$
- 可结合 `pedgenint` 与 Ne 估计，推导遗传进展 $\Delta G = i \cdot r_{IH} \cdot \sigma_A / L$
- 候选函数名：`pedselint()`

### 4.2 近交率与遗传多样性趋势

- 按年份/世代绘制 $\Delta F$、$\Delta c$（协亲缘变化率）趋势曲线
- 可在 `pedne` 基础上扩展，返回逐期 Ne 的时间序列
- 候选函数名：`pedtrend()`（或扩展 `vispstat()`）

### 4.3 亲本贡献均衡性评估

- 计算亲本实际贡献与理想均匀贡献的偏差（Gini 系数或变异系数）
- 帮助识别"超级亲本"对群体 $f_e$ / $f_a$ 的影响
- 可作为 `pedcontrib` 的扩展输出

### 4.4 迁移率 / 外血引入分析

- 基于 `pedancestry` 结果，量化每代外来品系血统比例的变化
- 候选函数名：`pedmigration()` 或扩展 `pedancestry()`

### 4.5 近交衰退分析（需外部数据）

- 需要用户提供表型数据，分析近交系数与性状表现的回归关系
- 可考虑对接 `purgeR` 包，本包仅提供系谱侧的辅助输出
- **建议**：本包内不实现，在文档中注明推荐搭配使用的包

### 4.6 遗传趋势（需 EBV）

- 需要用户提供估计育种值（EBV），绘制逐年/逐代遗传趋势图
- **建议**：本包内不实现，在文档或 vignette 中提供与 `optiSel` 对接的示例

---

## 五、命名规范一致性检查

按 [Positron.md](../Positron.md) 中约定的命名规范：

- [ ] 确认所有导出函数名均为全小写无下划线格式（`pedXxx`/`visXxx`）
- [ ] 确认所有返回 data.table 的列名均为 PascalCase
- [ ] 确认函数参数名均为全小写（`ped`, `by`, `reference`, `compact` 等）
- [ ] 检查 `pedrel` 返回列名 `NTotal`, `NUsed`, `MeanRel` 是否符合 PascalCase 规范 ✅
- [ ] 检查 `pedinbreed_class` 返回列名 `F_Class`：`_` 在列名中的使用是否符合规范（目前规范要求 PascalCase，`FClass` 更合规）
