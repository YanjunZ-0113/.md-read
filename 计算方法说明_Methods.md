# SaCas9 编辑效率计算方法说明

## 1. 用的是什么方法？

本项目的编辑效率 **不是** 简单比对 `.seq` 文本序列，而是基于 **Sanger 测序色谱图（.ab1）** 的 **indel 反卷积分析**。

实际使用的是：

| 项目 | 内容 |
|------|------|
| **工具名称** | **ICE**（Inference of CRISPR Edits） |
| **开发方** | Synthego |
| **开源地址** | https://github.com/synthego-open/ice |
| **Python 包** | `synthego-ice`（本地通过 `pip install synthego-ice` 安装） |
| **算法来源** | ICE 在 Python 中 **重新实现了 TIDE 算法**，并做了扩展 |

ICE 与 TIDE 的关系：

- **TIDE**（Tracking of Indels by Decomposition）是更早的 Sanger CRISPR 编辑分析网页工具
- **ICE** 参考 TIDE 思路，用 **非负最小二乘（NNLS）回归** 对混合色谱进行分解，估计各种 indel 的比例
- 本项目命令行调用 ICE，结果中的 **`ice` 字段 = 编辑效率（%）**

---

## 2. 参考文献与网站

### 2.1 TIDE（原始方法）

| 项目 | 信息 |
|------|------|
| **文章** | Brinkman EK, Chen T, Amendola M, van Steensel B. **Easy quantification of template-directed CRISPR/Cas9 editing**. *Nucleic Acids Research*, 2014, 42(22):e168. |
| **DOI** | https://doi.org/10.1093/nar/gku936 |
| **网页工具** | https://tide.nki.nl （荷兰癌症研究所 NKI） |
| **原理** | 将编辑后样本的 Sanger trace 与未编辑对照 trace 对齐，在切割位点后把混合信号分解为 WT + 多种 indel 模式的叠加 |

### 2.2 ICE（本项目实际使用）

| 项目 | 信息 |
|------|------|
| **文章** | Shen MW, et al. **Predictable and precise template-free CRISPR repair of repetitive elements**. 文中描述了 ICE 算法；预印本见 bioRxiv: https://www.biorxiv.org/content/10.1101/251082v1 |
| **GitHub** | https://github.com/synthego-open/ice |
| **在线工具** | https://ice.synthego.com |
| **与 TIDE 比较** | ICE 作者报告在多数样本上与 TIDE 编辑效率高度相关，且能处理更大 indel 等情况 |

### 2.3 其他相关

- **Biopython**：读取 `.ab1` 文件中的碱基序列与四色峰信号（`Bio.SeqIO.read(..., "abi")`）
- **SciPy NNLS**：ICE 内部用于分解 indel 贡献

---

## 3. 本项目具体计算流程

### 3.1 输入文件

**Upward（U1–U5）**

| 角色 | 文件 |
|------|------|
| Reference / Control | `controls for U1/Ruoyan_NT_F1_U_768982_A06.ab1` |
| 编辑样本 | `U1-5/Ruoyan_U1~U5_U2F_771103_*.ab1` |
| gRNA | 每个样本对应一条（见下表） |

**Downward（D1–D3）**

| 角色 | 文件 |
|------|------|
| Reference / Control | `D1-3/Ruoyan_NT_D3F_770920_A10.ab1` |
| 编辑样本 | `D1-3/Ruoyan_D1~D3_D3F_771103_*.ab1` |
| gRNA | 每个样本对应一条 |

### 3.2 gRNA 序列

**Upward**

| 靶点 | gRNA（5'→3'） |
|------|--------------|
| U1 | CCAAGTCCCAGACATGTCAG |
| U2 | CTAGGATTCTTAAAACTAGT |
| U3 | CAAAAACAATAAATAAAACG |
| U4 | TTTGGGGTAATCAGTGGGTA |
| U5 | TTAAAACTTTAAAGAGACAA |

**Downward**

| 靶点 | gRNA（5'→3'） |
|------|--------------|
| D1 | GAGAGAGAGGGAGTGAAAGA |
| D2 | GAGTGAGAGTGAGAGAGAGA |
| D3 | CGAAGAAGAGGGAAACCAAT |

### 3.3 计算步骤（ICE 算法）

对每个样本执行以下步骤：

```
步骤 1  读取色谱
        ├─ 读取 NT 对照 .ab1（wild-type trace）
        └─ 读取编辑样本 .ab1（edited trace）

步骤 2  定位 gRNA
        ├─ 在对照序列中搜索 gRNA（正向或反向互补）
        └─ 根据 SpCas9 默认规则估计切割位点（cut site ≈ PAM 上游 3 bp）

步骤 3  对齐（Alignment）
        ├─ 在切割位点上游选取高质量窗口
        └─ 将编辑样本 trace 与对照 trace 逐碱基对齐

步骤 4  生成 indel 假设（Edit proposals）
        ├─ 在切割位点附近枚举 -20 至 +N bp 的插入/缺失
        └─ 每种 indel 预测对应的理论色谱

步骤 5  分解（Decomposition）★ 核心
        ├─ 观测色谱 ≈ Σ (各 indel 模式贡献 × 理论色谱)
        └─ 用非负最小二乘（NNLS）拟合各 indel 比例

步骤 6  输出
        ├─ ice（%）= 编辑效率
        ├─ rsq（R²）= 拟合优度
        ├─ mean_discord_before / after = 切割位点前后信号不一致度
        └─ ice_d = 主导 indel 的碱基偏移
```

### 3.4 编辑效率的定义

在 ICE 输出中：

$$\text{Editing efficiency (\%)} = \text{ICE score} = (1 - f_{\text{WT}}) \times 100$$

其中 $f_{\text{WT}}$ 为 **零 indel（野生型）** 在分解中的相对比例。

等价理解：**切割位点处带有 indel 的等位基因所占比例**。

### 3.5 质量判据

| 指标 | 含义 | 本项目阈值 |
|------|------|-----------|
| **R²（rsq）** | 模型对色谱的解释度 | ≥ 0.9 认为可靠 |
| **discord_after** | 切割位点后 trace 与对照差异 | 越高越可能有编辑 |
| **discord_before** | 切割位点前差异 | 应较低；过高提示对齐问题 |

**本项目判读：**

- U4：25%，R²=0.99 → **可靠阳性**
- D3：0%，R²=0.99 → **可靠阴性**
- D1/D2：5%/7%，R²=0.05/0.07 → **不可靠**（重复区导致测序质量差）

---

## 4. 本地运行脚本

| 脚本 | 功能 |
|------|------|
| `analyze_u_editing.py` | U1–U5 ICE 分析 |
| `analyze_d_editing.py` | D1–D3 ICE 分析（含低质量样本放宽 QC） |
| `plot_editing_results.py` | U 柱状图 |
| `plot_d_editing_results.py` | D 柱状图 |
| `generate_ppt.py` | 生成汇报 PPT |

命令示例：

```bash
python analyze_u_editing.py
python analyze_d_editing.py
```

等价的 ICE 命令行（单个样本）：

```bash
synthego_ice \
  --control "controls for U1/Ruoyan_NT_F1_U_768982_A06.ab1" \
  --edited "U1-5/Ruoyan_U4_U2F_771103_D03.ab1" \
  --target "TTTGGGGTAATCAGTGGGTA" \
  --out results/U4_test
```

---

## 5. 注意事项（汇报/写文章时可引用）

1. **SaCas9 PAM**：ICE 默认按 **SpCas9 NGG PAM** 定位切割位点；SaCas9 实际 PAM 为 **NNGRRT**，因此可能出现 “No PAM upstream/downstream” 警告，但编辑效率估计在高质量 trace 上仍可用。

2. **不能用 .seq 直接算效率**：`.seq` 是单条 consensus，丢失套峰比例信息；必须用 `.ab1` 色谱。

3. **重复序列区**：D1/D2 位于 CT/GA 重复区，Sanger 易滑链，ICE 的 R² 极低，结果不应采信。

4. **与 TIDE 网页版等价**：若需交叉验证，可将相同 control ab1、sample ab1、gRNA 上传至 https://tide.nki.nl 或 https://ice.synthego.com 对比。

---

## 6. 汇报/论文 Methods 段落（可直接粘贴）

> Editing efficiencies were quantified from Sanger sequencing traces (.ab1) using the ICE (Inference of CRISPR Edits) algorithm (Synthego), which implements a TIDE-like decomposition of mixed chromatograms into wild-type and indel components by non-negative least squares regression. Non-targeting (NT) controls served as wild-type references: Ruoyan_NT_F1_U_768982_A06 for upward targets (U1–U5) and Ruoyan_NT_D3F_770920_A10 for downward targets (D1–D3). The guide RNA sequence for each target was provided to ICE to localize the cut site. Editing efficiency was reported as the ICE score (% indels). Results with R² < 0.9 were considered unreliable. The original TIDE method is described by Brinkman et al. (Nucleic Acids Res, 2014; doi:10.1093/nar/gku936). ICE is described at https://github.com/synthego-open/ice.

中文版本：

> 编辑效率通过 Sanger 测序色谱图（.ab1）使用 ICE（Inference of CRISPR Edits）算法计算。ICE 采用与 TIDE 类似的原理，利用非负最小二乘回归将切割位点后的混合色谱分解为野生型与各 indel 类型的叠加。Upward 靶点以 Ruoyan_NT_F1_U_768982_A06 为野生型对照，Downward 靶点以 Ruoyan_NT_D3F_770920_A10 为对照。每个靶点输入对应 gRNA 序列以定位切割位点。编辑效率以 ICE score（% indels）报告；R² < 0.9 的结果视为不可靠。TIDE 原始方法见 Brinkman et al., Nucleic Acids Research, 2014 (doi:10.1093/nar/gku936)；ICE 见 https://github.com/synthego-open/ice。
