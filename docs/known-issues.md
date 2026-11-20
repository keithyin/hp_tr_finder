# 已知问题：ref 输出的冗余嵌套与边界重叠注释

- **状态**：问题一已修复（0.6.0，`--dedupNested` 启用，默认关闭）；问题二未修复
- **发现日期**：2026-11-20
- **发现版本**：[0.5.0]（commit `034a1c2`）
- **影响场景**：仅 `ref`（`called` 场景的区间会经 `merge_region` 合并，不受此问题影响）
- **坐标约定**：本文所有坐标均与 GFF 输出一致——第 4 列 start 为 1-based，第 5 列 end 为 0-based（沿用 tr-finder 风格）

## 现象总览

对 `hp_tr_finder ref str-seqs.fa --unitAndRepeats 1-3,2-3,3-3,4-3` 的输出
（`/home/keith/Datasets/str-seqs/str-seqs.gff`，70 条序列、3770 条注释）做区间扫描：

| 类别 | 定义 | 数量 | 状态 |
|---|---|---|---|
| 完全嵌套 | 一条注释的区间整体包含另一条注释 | 60 对，涉及 10 个外层 TR 位点 | ✅ 0.6.0 修复（`--dedupNested` 启用），见[修复说明](#修复说明060) |
| 部分重叠 | 两条注释区间相交但互不包含（共享 1–5 bp） | 196 对 | 未修复 |

两类都不在 0.5.0 去重（`dedup_rotation_hits`）的处理范围内，见
[为什么 0.5.0 的去重没有消除它们](#为什么-050-的去重没有消除它们)。

## 问题一：外层 TR 内部嵌套同聚物注释（已修复，0.6.0）

### 机理

`--unitAndRepeats 1-3,2-3,3-3,4-3` 同时启用 unit-1 同聚物正则与 unit-4 正则。
当 unit-4 motif **内部本身含 ≥3 连同种碱基**（`GAAA` 含 `AAA`、`CCCT` 含 `CCC`、
`CGCC` 跨单元拼出 `CCC`）时，每个重复单元都会再现这段同聚物，unit-1 正则
`(A){3,}` / `(C){3,}` 便在大 TR 内部独立命中一次，产出与外层注释完全嵌套的
同聚物注释——每个重复单元一条。

0.5.0 的 primitive motif 过滤（`generate_motif_core`，src/lib.rs:51）只拦截
"整体是短 unit 串联重复"的 motif（如 `ATAT`、`AAAA`），拦不住这种
"局部含同聚物"的 primitive motif（`GAAA` 不是任何更短单位的串联，按定义是最简的）。

本数据集中 60 对嵌套全部是 **unit-4 ⊃ unit-1**。unit-3 motif（如 `GAA`）内部
只产生 2 连碱基 run，达不到 min_repeats=3，因此不产生嵌套、只产生问题二的边界重叠。

### 全部 10 个嵌套位点

| 序列 | 外层注释（GFF 坐标） | 外层序列片段 | 嵌套内层 |
|---|---|---|---|
| 260805STR48-4 | 328–399 `(AGAA)18` | `AGAA`×18 | 17 × `(A)3` |
| STR33-5 | 51–102 `(GAAA)13` | `GAAA`×13 | 13 × `(A)3` |
| STR33-11 | 108–159 `(GAAA)13` | `GAAA`×13 | 13 × `(A)3` |
| STR33-6 | 1–20 `(AAGA)5` | `AAGA`×5 | 4 × `(A)3` |
| 260731STR37-1 | 209–220 `(TCCC)3` | `TCCCTCCCTCCC` | 3 × `(C)3` |
| 260827STR-93-1 | 265–276 `(GGGT)3` | `GGGTGGGTGGGT` | 2 × `(G)3` |
| STR-109-2 | 260–271 `(CCCT)3` | `CCCTCCCTCCCT` | 2 × `(C)3` |
| STR-114-1 | 807–818 `(CGCC)3` | `CGCCCGCCCGCC` | 2 × `(C)3` |
| STR-114-3 | 787–798 `(CGCC)3` | `CGCCCGCCCGCC` | 2 × `(C)3` |
| STR33-11 | 38–49 `(AAGA)3` | `AAGAAAGAAAGA` | 2 × `(A)3` |

（内层坐标例：STR33-5 的 13 条 `(A)3` 为 52–54、56–58、…、100–102，间隔 4bp 排满外层区间。）

### 典型例子

**STR33-5（一条大 TR 套出 13 条冗余）**，GFF 原文摘录：

```
STR33-5  .  .  51  102  .  .  .  Note=(GAAA)13,GAAAGAAAGAAA...GAAA
STR33-5  .  .  52   54  .  .  .  Note=(A)3,AAA      ← 每个单元的 AAA
STR33-5  .  .  56   58  .  .  .  Note=(A)3,AAA
...
STR33-5  .  .  100 102  .  .  .  Note=(A)3,AAA
```

序列 `GAAAGAAAGAAA…`：每个 `GAAA` 单元的尾部 `AAA` 都满足 `(A){3,}`，
于是 13 个单元产出 13 条与外层 `(GAAA)13` 完全嵌套的注释。

**STR-114-1（跨单元拼接的 run）**：`CGCCCGCCCGCC` 中 `CCC` 分别出现在
单元内部（第 3–5 位）与单元交界处（第 7–9 位），产出 2 条 `(C)3`，
嵌套于 `(CGCC)3` [807,818) 内。

## 修复说明（0.6.0）

**启用方式**：默认关闭，CLI 加 `--dedupNested` 启用；lib 侧
`all_seq_hp_tr_finder` / `single_seq_hp_tr_finder` 增加尾参 `dedup_nested: bool`。
实现位于 `dedup_rotation_hits` 之后的 `dedup_nested_hits`（src/lib.rs），
unit / motif 两种模式与 ref / called 两个 scenario 共用同一 pipeline。

**裁决规则**：对每个包含对（外层 ⊃ 内层），比较拷贝数
`copies = (end - start) / motif 长度`：

- 内层 copies > 外层 → 仅保留内层（删外层）；
- 内层 copies ≤ 外层 → 保留外层（删内层）。

当前数据 60 对全部为"内层 3 ≤ 外层 3–18"，即全部删内层、留外层。
">" 分支为规则完备性保留：假想的 `(TGGGG)3`（每单元含 `GGGG` run）⊃ 3×`(G)4`
会删外层、仅留 3 条 `(G)4`。理论保证：同一外层内所有嵌套同聚物的 copies 一致，
且不同周期的两条精确 TR 互相嵌套被 Fine–Wilf 定理排除（强制公共周期 →
motif 非最简），故逐对独立标记不会产生裁决冲突。

**修复后验证**（str-seqs 数据集，`--unitAndRepeats 1-3,2-3,3-3,4-3`）：

- 不加 flag：ref / called 输出与 0.5.0 逐字节一致；
- 加 `--dedupNested`：ref 恰删 60 条嵌套注释（3770 → 3710），删除行全部为
  `(A)3`/`(C)3`/`(G)3`，按 10 个位点的分布与本文表格逐一吻合；外层 TR 全部保留；
- 问题二不受影响：`(C)5`∩`(CCCT)3`、`(GGGT)3`∩`(GT)6` 等边界共享重叠原样保留；
- called 场景加 flag 输出亦与 0.5.0 一致（被删内层完全包含于保留外层，
  `merge_region` 合并结果不变）。

## 问题二：边界共享型部分重叠（未修复）

### 机理

每种 motif 的正则独立做贪婪最长匹配，各自的最优起点可以落在相邻命中的内部
——典型形态是"`N 连同种碱基` + 紧跟 `TR`"，TR 的首个重复单元"借用"同聚物
run 尾部的碱基。例如 `GGGAGAGAGAGA` 中 `(GA){3,}` 的最优起点是第 3 个 G，
与 `(G){3,}` 命中的 `GGG` 共享 1 bp。两个正则谁也不让谁，重叠部分被双方同时报告。

### 按 unit 组合分布（无序）

| 组合 | 对数 | 典型共享长度 |
|---|---|---|
| unit1 × unit2 | 115 | 1 bp（G-run 尾碱基） |
| unit1 × unit3 | 45 | 1–2 bp |
| unit1 × unit4 | 11 | 1–2 bp |
| unit2 × unit2 | 5 | 1 bp |
| unit2 × unit3 | 7 | 3 bp（整段） |
| unit3 × unit3 | 10 | 1–2 bp |
| unit3 × unit4 | 2 | 5 bp |
| **合计** | **196** | |

### 典型例子（覆盖全部组合类型，GFF 原文）

**unit1 × unit2** — `260623STR1-1`：`GGGAGAGAGAGA`
```
260623STR1-1  .  .  317 319  .  .  .  Note=(G)3,GGG
260623STR1-1  .  .  319 328  .  .  .  Note=(GA)5,GAGAGAGAGA   ← 共享第 319 位的 G
```

**unit1 × unit3** — `260623STR1-1`：`GGGGAGGAGGA`
```
260623STR1-1  .  .  122 125  .  .  .  Note=(G)4,GGGG
260623STR1-1  .  .  124 132  .  .  .  Note=(GGA)3,GGAGGAGGA   ← 共享 GG
```

**unit1 × unit4** — `STR33-11`：`AGGAAGGAAGGAAA`
```
STR33-11  .  .  231 242  .  .  .  Note=(AGGA)3,AGGAAGGAAGGA
STR33-11  .  .  242 244  .  .  .  Note=(A)3,AAA               ← 共享第 242 位的 A
```

**unit2 × unit2** — `STR27-1`：`CACACATATATATAT`
```
STR27-1  .  .  1  6  .  .  .  Note=(CA)3,CACACA
STR27-1  .  .  6 15  .  .  .  Note=(AT)5,ATATATATAT           ← 共享第 6 位的 A
```
（`STR27-1` 头部是连续化合物微卫星，`(CA)3|(AT)5|(AC)6|(AT)3|(AC)7|…`
每对相邻注释共享 1 bp。）

**unit2 × unit3** — `STR33-11`：`CTCTCTTCTTCT`（共享长度达 3 bp）
```
STR33-11  .  .  848 853  .  .  .  Note=(CT)3,CTCTCT
STR33-11  .  .  851 859  .  .  .  Note=(TCT)3,TCTTCTTCT       ← 共享 TCT
```

**unit3 × unit3** — `260805STR18-1`：`AAGAAGAAGAAGGAGGAGGAGG`
```
260805STR18-1  .  .  230 241  .  .  .  Note=(AAG)4,AAGAAGAAGAAG
260805STR18-1  .  .  240 251  .  .  .  Note=(AGG)4,AGGAGGAGGAGG  ← 共享 AG
```

**unit3 × unit4** — `STR-114-1`：`CGCCCGCCCGCCGCCGCCG`（共享 5 bp）
```
STR-114-1  .  .  807 818  .  .  .  Note=(CGCC)3,CGCCCGCCCGCC
STR-114-1  .  .  814 825  .  .  .  Note=(CCG)4,CCGCCGCCGCCG   ← 共享 CCGCC
```

### 嵌套与重叠可出现在同一位点

`STR-109-2`：序列 `AACCCCCTCCCTCCCTTCC` 中前导 `CC` 与 `(CCCT)3` 首单元的
`CCC` 连成 5 连 C——首单元的 `CCC` 与外层是**重叠**而非嵌套，后两个单元的
`CCC` 则是**嵌套**：
```
STR-109-2  .  .  258 262  .  .  .  Note=(C)5,CCCCC
STR-109-2  .  .  260 271  .  .  .  Note=(CCCT)3,CCCTCCCTCCCT   ← 与 (C)5 重叠共享 CCC
STR-109-2  .  .  264 266  .  .  .  Note=(C)3,CCC               ← 完全嵌套
STR-109-2  .  .  268 270  .  .  .  Note=(C)3,CCC               ← 完全嵌套
```

## 为什么 0.5.0 的去重没有消除它们

`dedup_rotation_hits`（src/lib.rs:225）只消除**同周期旋转冗余**：
同周期 p 的命中，起点错位 shift < p 视为同一 TR 的循环旋转（如 `(AT)6` 与
`(TA)7`），保留最靠左的一条。上述两类都在其保留规则之内：

1. **嵌套**是跨周期命中（周期 4 ⊃ 周期 1）。设计上有意保留跨周期重叠
   （src/lib.rs:209–214 注释："跨周期重叠（如 (AT) 与 (ATG)）是真实的不同注释"），
   由 `test_dedup_cross_period_overlap_kept` 覆盖。该设计对
   `(AT)`/`(ATG)` 这类**不同 motif 的真实不同注释**是对的，
   但没有考虑到"内层命中是外层 motif 内部同聚物"这种**派生冗余**。
2. **边界共享**是同周期或跨周期、shift ≥ p 的命中，被判定为"真实的不同 TR"而保留
   （`test_shift_ge_period_kept_even_if_overlapping`、
   `test_dedup_boundary_sharing_distinct_trs_kept`）。该判定对
   `AAT×2 0-6 与 TGG×2 5-11` 这类真实相邻 TR 是对的，但
   `(GA)5` 起点借走 `(G)3` 尾碱基这类**共享派生碱基**的情形同样是冗余。

即：0.5.0 的去重是"保守正确"的——宁可多报也不错杀；本文件记录的冗余
是它有意放行的两类。（0.6.0 起问题一可通过 `--dedupNested` 启用嵌套去重解决，
见[修复说明](#修复说明060)；问题二仍保持该默认行为。）

## 影响

- `ref` 输出中一个 locus 最多产生 18 条注释（1 条外层 + 17 条内层，见
  260805STR48-4）；按注释条数统计 TR/HP 数量会被显著高估。
  （问题一的影响可经 `--dedupNested` 消除，见[修复说明](#修复说明060)。）
- 内层注释与外层注释对同一段碱基给出两种周期解释（4-mer vs 1-mer），
  下游若直接取"每条注释一个 motif"，会得到互相矛盾的等位基因结构。
- `called` 场景不受影响（区间合并后输出），但 motif 信息在该场景本就不输出。

## 候选修复方向（待讨论）

1. **嵌套过滤**（✅ 已采用，0.6.0 `--dedupNested`，按拷贝数裁决而非仅按
   motif 包含关系，见[修复说明](#修复说明060)）：删除"完全包含于另一条注释、
   且自身 motif 是外层 motif 内部同聚物 run"的内层注释（如 `(A)3` ⊂
   `(GAAA)13` 且 `GAAA` 含 `AAA`）。注意与真实嵌套 STR 区分：`(A)3` ⊂
   `(AT)10` 的情形不存在（`AT` 内无 `AA` run），但更复杂的化合物 locus
   是否要保留内层注释需要定义。
2. **边界借碱基处理**：对 "同聚物 run + TR" 相接处，截去 TR 命中借用的
   共享碱基（改为 touching 不相交），或规定同聚物优先、TR 从 run 之后起匹配。
   风险：部分组合（如 `CTCTCTTCTTCT` 的 `(CT)3`/`(TCT)3` 共享 3 bp 整段）
   属于周期结构歧义，截断方向不唯一。
3. **与上游 tr-finder v2 对齐**：确认参考实现的期望行为，避免单方面
   引入不一致的注释语义。

## 复现方法

```bash
cd /home/keith/projects/hp_tr_finder
cargo build --release
./target/release/hp_tr_finder ref /home/keith/Datasets/str-seqs/str-seqs.fa \
    --unitAndRepeats 1-3,2-3,3-3,4-3 -o /tmp/out.gff
```

统计方法：解析 GFF 中每条注释的 (start, end)（start 减 1 转 0-based 后
为左闭右开区间），对每条序列内的注释按起点排序后做两两扫描，
判定完全包含（`b.start >= a.start && b.end <= a.end`）与部分相交，
并从 FASTA 中取对应片段核对序列。文中所有序列均已经过直接核对。
