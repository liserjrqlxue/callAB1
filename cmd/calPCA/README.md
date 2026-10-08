# calPCA — 基因合成 Sanger 测序批量质控分析

`calPCA` 是基于 [tracy](https://github.com/gear-genomics/tracy) 的批量 Sanger 测序（ab1）分析工具，用于基因合成（自合）批次的克隆质控：对每个分段并行执行碱基识别、参考比对与变异分解，在引物、片段、批次三个维度统计变异与准确率，判定正确克隆，并输出 Excel 报表；可选执行化学补充（化补2）引物设计。

## 功能特性

- 从设计 Excel（`自合.xlsx`）加载 **基因 → 分段 → 引物对** 三级序列结构
- 按分段并发调用 tracy（`basecall` → `align` → `decompose`），信号量控制线程数
- 支持两种测序文件组织方式：
  - **CY0130 模式**：按分段 ID 在测序目录（含一级子目录）glob `*.ab1`，支持 `rename.txt` 映射
  - **常规模式**：按 `BGA-<tag>-TY1-<克隆号>-T7(-Term).ab1` 固定命名、固定克隆数配对
- 变异三级汇总：同一克隆跨 Sanger 合并 → 同一变异跨克隆合并 → 按位置合并
- 统计 SNV / Insertion / Deletion / SV 数量与比率、变异比率、准确率、参考收率、参考单步准确率（几何平均）
- 判定有效克隆、正确克隆（及杂合正确克隆）、高频变异、序列性缺失，并给出片段/批次合格结论
- 输出多 Sheet 的 `<输出目录>.Sanger结果.xlsx` 及若干中间 txt
- `-fix` 对无正确克隆的分段重新做引物设计，产出化补2清单/引物订购单并打包 zip

## 目录结构与源码职责

| 文件 | 职责 |
| --- | --- |
| [main.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/main.go) | 入口；flag 定义、`Seq` 结构、并发调度、`processSeq` 单片段分析主流程、变异汇总 |
| [tool.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go) | tracy 批量封装（两种模式）、Excel 序列加载、引物/片段维度的变异与准确率统计、阈值常量 |
| [sanger.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/sanger.go) | 结果 Excel 各 Sheet 写入（Sanger结果、片段结果、板位图、批次统计、拼接引物板）、txt 输出、批次判定 |
| [simple.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/simple.go) | `-fix` 化学补充流程（调用 PrimerDesigner 设计引物、写清单/订购单、打包 zip） |
| [global.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/global.go) | 各报表表头、默认克隆数（32）、默认质量阈值（45）、96 孔板行列定义 |
| [xlsx.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/xlsx.go) | excelize 坐标转单元格名工具函数 |
| [pkg/tracy](file:///e:/liserjrqlxue/callAB1/pkg/tracy) | tracy 二进制封装与结果解析（`exec.go`、`result.go`、`global.go`） |

## 依赖

- **Go** ≥ 1.23（module `callAB1`）
- **tracy 二进制**：默认取 PATH 中的 `tracy`，可用 `-tracy` 指定路径。目录内附 Linux 版 `tracy_v0.7.8_linux_x86_64bit`
- **AZENTA.xlsx**：化补2引物订购单模板，须位于可执行文件同目录（`-fix` 时使用）
- 主要 Go 依赖：`excelize/v2`、`samber/lo`、`goUtil`、`PrimerDesigner/v2`、`liserjrqlxue/DNA`

## 输入

### 1. 自合.xlsx（`-i`，默认 `<输出目录>/自合.xlsx`）

| Sheet | 关键列 | 说明 |
| --- | --- | --- |
| `原始序列` | 基因名称、DNA序列 | 基因级序列 |
| `分段序列` | 片段名称、片段序列、起点、终点 | 分段序列；`RefID = 片段名称去掉末2字符`；固定跳过 `EGC_K` |
| `引物对序列` | 引物对名称、引物对序列、有效序列-起点/终点、左引物-起点/终点、右引物-起点/终点、片段-起点/终点 | 引物对及其左右引物坐标，加载时校验序列一致性 |

### 2. 引物订购单（`-io`，默认 `<输出目录>/引物订购单_BOM.xlsx`）

读取其中的 `拼接引物板` Sheet，在结果工作簿中重建并追加“平均准确率 / 参考单步准确率”行。

### 3. 测序 ab1（`-s`，默认 `<输出目录>/ab1`）

测序文件与片段/克隆的具体对应规则见下方 [序列层级与 Sanger 文件对应关系](#序列层级与-sanger-文件对应关系)。

### 4. rename.txt（`-r`，默认 `<输出目录>/rename.txt`）

tab 分隔两列、无表头：`片段名称<TAB>sanger文件前缀`。用于 CY0130 模式下把设计表中的片段名桥接到实际测序文件名；详细规则与不规则文件名示例见 [rename.txt 用法](#renametxt-用法不规则文件名)。

> 注意：只要 `-r` 非空即进入 CY0130 模式，而该 flag 有默认值，因此**程序默认运行在 CY0130 模式**；常规模式需显式 `-r ""`。

## 序列层级与 Sanger 文件对应关系

### 1. Excel 中的层级与 ID 命名

加载后内存中为四级树：**基因 → 片段 → 引物对 → 引物**，每一级靠 ID 的后缀规则反查父级（无需额外的 ID 列）：

| 层级 | Sheet · 列 | ID 示例 | 父级反查规则 |
| --- | --- | --- | --- |
| 基因 | `原始序列` · 基因名称 | `BGA_14A2` | — |
| 片段 | `分段序列` · 片段名称 | `BGA_14A2_1` | `RefID = 片段名去掉末 2 字符`（`BGA_14A2_1` → `BGA_14A2`）；固定跳过 `EGC_K` |
| 引物对 | `引物对序列` · 引物对名称 | `BGA_14A2_1_1` | 去掉末 1 字符后再去掉尾部 `_`（`BGA_14A2_1_1` → `BGA_14A2_1`） |
| 引物 | 程序派生，无 Sheet | 常规模式：`BGA_14A2_1_1_1`（左）、`..._2`（右）；CY0130 模式：仅 1 条，ID 同名引物对 | 由"左/右引物-起点/终点"列生成 |

挂载关系（[tool.go](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go)）：`gene.SubSeq = 该基因的全部片段`（[LoadSegmentSequence](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go#L637-L682)）；`segment.SubSeq = 引物对`、`pair.SubSeq = 左/右引物`（[LoadPrimerPairSequence](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go#L684-L739)）。

### 2. 序列坐标如何对应

- 片段行的 `起点/终点`：片段有效区在片段序列上的 0-based 半开区间 `[起点, 终点)`
- 引物对行的坐标：
  - `有效序列-起点/终点`：引物对**自身序列**上的有效区
  - `片段-起点/终点`：引物对有效区在**片段有效区内**的偏移
  - `左引物-起点/终点`、`右引物-起点/终点`：左/右引物在引物对序列上的坐标
- 加载时强校验（不一致直接 `log.Fatal`）：

```text
片段Seq[起点:终点][片段-起点:片段-终点] == 引物对Seq[有效序列-起点:有效序列-终点]
```

- 统计变异时（[RecordPrimer / RecordPair](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go#L171-L277)）把引物坐标换算回片段坐标：`offset = 片段起点 + 片段-起点 − 有效序列-起点`，只统计落在 `(引物Start+offset, 引物End+offset]` 内的变异。tracy 输出的 `Pos` 是相对该片段参考 fasta（即 `<id>.fa`）的位置。

### 3. Sanger 文件如何归属到片段和克隆

分析以**片段**为单位：每个片段写一个参考 fasta，tracy 以它为 reference 分析该片段的全部克隆；一个克隆可含 1～2 条 Sanger（T7 端 / term 端），两条中任一 `PASS` 该克隆即为有效克隆。

#### CY0130 模式（默认，[RunTracyBatchCy0130](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go#L87-L104)）

1. 在两处按前缀 glob 文件（递归深度仅一级子目录）：
   - `<sanger目录>/<前缀>*.ab1`
   - `<sanger目录>/*/<前缀>*.ab1`
2. 文件名匹配 `term.ab1`（小写）的文件不作为独立克隆；其余每个文件 = 一个克隆，**克隆 ID = 文件名去掉 `.ab1`**
3. term 端自动配对（[RunTracyCY0130](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go#L44-L71)）：同路径下把主文件名中的 `T7.ab1` 替换为 `T7term.ab1`，存在则作为该克隆第 2 条；不存在则按单端分析

```text
ab1/
├── BGA_14A2_1-1-T7.ab1        ── 克隆 "BGA_14A2_1-1-T7" 第1条
├── BGA_14A2_1-1-T7term.ab1    ── 同一克隆第2条（自动配对）
├── BGA_14A2_1-2-T7.ab1        ── 克隆 "BGA_14A2_1-2-T7"（单端）
└── plate1/
    └── BGA_14A2_2-1-T7.ab1    ── 一级子目录同样会被搜到
```

#### 常规模式（`-r ""`，[RunTracyBatch](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/tool.go#L74-L82)）

- tag 推导：片段 ID 按 `_` 切分取第 2 段前 4 字符、第 3 字符替换为 `B`。`BGA_14A2_1` → `14A2` → **`14B2`**
- 克隆号从 1 扫到 `-c`（默认 32），每个克隆号固定查找一对文件：

```text
BGA-14B2-TY1-1-T7.ab1
BGA-14B2-TY1-1-T7-Term.ab1     # 注意此模式是大写 Term
```

- 两个文件都存在才分析；**只找到一个会 `log.Fatal`**；都没有则跳过该克隆号。克隆 ID 即 `"1".."N"`

### 4. rename.txt 用法（不规则文件名）

格式（tab 分隔、无表头）：

```text
<分段序列中的片段名称>	<sanger文件实际前缀>
```

行为规则（[main.go L180-216](file:///e:/liserjrqlxue/callAB1/cmd/calPCA/main.go#L180-L216)）：

1. **文件存在时是白名单**：只有第一列列出的片段才会分析；未列出的片段打印 `Skip` 错误日志并跳过
2. **文件不存在时**：自动对全部片段建立恒等映射（前缀 = 片段名）
3. 第二列前缀只用于 glob 匹配（根目录 + 一级子目录），匹配到的每个文件仍以其文件 stem 作为克隆 ID
4. term 端必须是小写 `T7term.ab1` 且与主文件同目录才会自动配对

**示例**：测序目录结构不规则、文件放在不同板子目录：

```text
ab1/
├── plate1/
│   ├── CY0130-A01-1-T7.ab1
│   ├── CY0130-A01-1-T7term.ab1
│   └── CY0130-A01-2-T7.ab1
└── plate2/
    └── CY0130-A02-1-T7.ab1
```

`rename.txt`：

```text
BGA_14A2_1	CY0130-A01
BGA_14A2_2	CY0130-A02
```

则片段 `BGA_14A2_1` 会匹配到 `plate1/CY0130-A01-1(-T7/-T7term)` 与 `...-2-T7` 两个克隆；未在文件中出现的片段全部跳过。

**对结果显示的影响**：正确克隆在"片段结果/板位图"中会尝试把克隆 ID 缩写成 `片段ID-克隆号` 形式，要求文件名匹配 `^<片段ID>-(\d+)`，或片段 ID 首个 `_` 写作 `-` 的形式（如 `BGA-14A2_1-3-T7.ab1`）。完全不规则的文件名（如上例 `CY0130-A01-1-T7`）**不影响变异统计、正确克隆计数与合格判定**，仅列表中保留完整文件名并打印一条 `can not parse clone` 日志。

**前缀冲突提醒**：匹配规则是 `前缀*.ab1`，各前缀不能互为前缀（如 `CY0130-A1` 会误匹配到 `CY0130-A10-...`），请使用带分隔符、足够区分度的前缀。

## 参数说明

| flag | 类型 | 默认值 | 说明 |
| --- | --- | --- | --- |
| `-o` | string | — | **必填**，输出目录 |
| `-i` | string | `<o>/自合.xlsx` | 自合设计 Excel |
| `-io` | string | `<o>/引物订购单_BOM.xlsx` | 引物订购单 BOM |
| `-s` | string | `<o>/ab1` | ab1 测序目录 |
| `-r` | string | `<o>/rename.txt` | CY0130 重命名文件；置空可切换到常规模式 |
| `-tracy` | string | `tracy` | tracy 二进制路径 |
| `-c` | int | `0`（用内置值 32） | 常规模式每分段扫描的最大克隆数 |
| `-q` | int | `0`（用内置值 45） | 变异统计的最低 Qual，低于该值的变异被过滤 |
| `-w` | bool | `false` | 覆盖已有的 tracy 结果 JSON（否则断点复用） |
| `-fix` | bool | `false` | 分析完成后执行化学补充（化补2）引物设计 |
| `-t` | int | `0`（=分段总数） | 并发线程上限 |

## 处理流程

1. 加载 `原始序列 / 分段序列 / 引物对序列`，建立 `基因 → 分段 → 引物对 → 引物` 结构
2. 读取 `rename.txt`（或恒等映射）
3. 信号量 + WaitGroup 并发，每个分段执行 `processSeq`：
   1. 写出 `<片段ID>.fa`
   2. 按模式调用 tracy 批量分析，得到 `克隆号 → [正向(, 反向)] *tracy.Result`
   3. 汇总每条 Sanger 的变异，写 `<id>.variant.raw.txt`
   4. 仅保留 `PASS` 且 `Qual ≥ MaxQual` 的变异：
      - 同克隆跨 Sanger 合并 → `<id>.variant.clone.txt`
      - 同变异跨克隆合并（含 CloneCount / ClonePass / Ratio）→ `<id>.variant.set.txt`
      - 同位置合并 → `<id>.variant.pos.txt`
   5. 片段级统计 `RecordSeq` → `<id>.seq.result.txt`；引物级统计 `RecordSeqPrimer` → `<id>.result.txt`
4. 汇总各分段结果，计算批次统计与合格判定
5. 输出 `Result.txt`、`TracyResult.txt` 与 `<o>.Sanger结果.xlsx`
6. `-fix` 时对无正确克隆（`CloneHit == 0`）的分段追加 `C` 后缀重新设计引物，输出化补2文件并打包 `<批次名>.calPCA.zip`

### tracy 单条分析（pkg/tracy）

对每个 ab1 依次执行 `basecall`、`align`、`decompose`，结果缓存为 JSON，`-w` 强制重算：

- 比对边界：测序有效区域 `[60, 700]`，参考边界 `[40, 长度-40]`
- `BoundMatchRatio < 0.9` → 状态追加 `LowMatch`；basecall 长度 `< 100` → `TooShort`；无任何标记即为 `PASS`
- SV 判定：类型为 `Complex`，或 SNV 的 Alt 长度 ≥ 3，或插入/缺失长度差 ≥ 3

## 判定规则与指标

- **有效 Sanger / 有效克隆**：同一克隆正反向结果中至少一条 `PASS`
- **正确克隆**：有效且片段区域 `[起点, 终点]` 内变异数为 0；区域内变异全部为杂合（`het.`）记为**杂合正确克隆**
- **变异统计口径**：仅统计 `PASS`、`Qual ≥ MaxQual(45)` 且落在目标区域内的变异
- **变异位置比率**：位置变异克隆数 / 有效克隆数
- **参考收率 (%)**：`∏(1 - 各位置变异比率) × 100`
- **参考单步准确率 (%)**：`收率^(1/片段长度) × 100`（几何平均）
- **变异/准确率 (%)**：`变异个数×100/(长度×有效Sanger数)`，准确率 = 100 − 变异比率
- **SNV/INS/DEL/SV 比率 (%)**：各类型变异数 × 100 / (长度 × 有效Sanger数)
- **高频变异**：以变异位置为中心、窗口 k=1 求加权和，按和降序贪心选取互不重叠窗口，入选条件为窗口和 ≥ 4 或 窗口和/有效克隆数 ≥ 0.6
- **序列性缺失**：缺失碱基位置在有效克隆中缺失比率 ≥ 0.8，按碱基数计数
- **片段阈值判定**：无有效克隆，或 `SNV ≥ 7.22%`、`INS ≥ 0.40%`、`DEL ≥ 0.28%` 任一满足 → `不合格`
- **批次判定**：有效克隆为 0，或批次内有效片段的平均 SNV/INS/DEL 比率超过同上阈值 → `不合格`

## 输出

### 输出目录内文本

| 文件 | 内容 |
| --- | --- |
| `<id>.fa` | 分段参考序列 |
| `<id>_<克隆>.basecall.json / .align.json / .decompose.json` | tracy 原始结果（复用缓存） |
| `<id>_<克隆>.align.summary.txt` | 比对边界匹配率摘要 |
| `<id>_<克隆>.decompose.variants.txt` | 单条 Sanger 变异明细 |
| `<id>_<克隆>.stdout.txt / .stderr.txt` | tracy 运行日志 |
| `<id>.variant.raw.txt` | 全部克隆 × Sanger 原始变异 |
| `<id>.variant.clone.txt` | 克隆级合并变异 |
| `<id>.variant.set.txt` | 跨克隆变异集合统计 |
| `<id>.variant.pos.txt` | 位置级变异统计 |
| `<id>.result.txt` | 引物维度结果（`ResultTitle`） |
| `<id>.seq.result.txt` | 片段维度结果（`SeqTitle`） |
| `Result.txt` | 全批次引物维度汇总 |
| `TracyResult.txt` | 全部 Sanger 的 tracy 状态汇总 |

### `<输出目录>.Sanger结果.xlsx`

| Sheet | 粒度 | 内容 |
| --- | --- | --- |
| `Sanger统计` | 每条 Sanger | 状态、PASS、变异数、杂合数、比对率 |
| `Sanger结果` | 每条引物 | 有效 Sanger 数、各类变异数/比率、准确率/收率/单步准确率 |
| `片段结果` | 每个分段 | 长度、有效/无效 Sanger、正确克隆数、变异统计、阈值判定等（`SeqTitle`） |
| `测序结果板位图` | 基因 × 96 孔板 | 基因长度/节数及各分段正确克隆孔位标注 |
| `Clone变异结果` | 克隆 × 变异 | Qual/Filter/Genotype 合并结果 |
| `变异统计` | 变异集合 | 跨克隆的 CloneCount/ClonePass/CloneRatio |
| `拼接引物板` | 板位 | 订购单板位图 + 平均准确率 + 参考单步准确率（提供 `-io` 时） |
| `批次统计` | 批次 | 总/有效片段与克隆数、高频变异、序列性缺失、平均比率、合格状态 |

### `-fix` 追加产物

- `<批次名>-化补2清单.xlsx`（`化补2清单`、`化补2引物板位图`）
- `<批次名>-化补2引物订购单.xlsx`（基于可执行文件同目录 `AZENTA.xlsx`）
- `<批次名>.calPCA.zip`（Sanger结果 + 上述两个文件）

## 构建与运行

```bash
# 构建
go build -o calPCA ./cmd/calPCA

# CY0130 模式（默认，测序文件按 <片段ID>*.ab1 组织；-s 可省略，默认 <o>/ab1）
./calPCA -o path/to/batch -t 16

# 常规模式（固定 BGA 命名 + 成对 T7/T7-Term）
./calPCA -o path/to/batch -s path/to/ab1 -r "" -c 32

# 分析后追加化学补充
./calPCA -o path/to/batch -fix -w
```
