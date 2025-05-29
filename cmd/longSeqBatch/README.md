# 批量分析长片段的多段 sanger 测序

## 有效克隆

- `sanger` 内置有效条件
- 序列覆盖完全

## 输入输出结构

- `-i` `输入目录`
  - `-s` `测序目录`
    - `*.ab1` 测序文件
  - `-r` `目标序列.xlsx`
- `-d` `输出目录`
  - `-o` `result.xlsx`
  - `ref/[GeneID].ref.fa` 基因序列
  - `[GeneID]/[CloneID]/[Index]*` ab1 分析结构

## 参数说明

| flag | type     | default         | note                         |
| ---- | -------- | --------------- | ---------------------------- |
| `-i` | `string` |                 | **必填**，`输入目录`         |
| `-s` | `string` | `测序目录`      | `输入目录` 下的测序子目录     |
| `-r` | `string` | `目标序列.xlsx` | `输入目录` 下的输入信息 excel |
| `-d` | `string` | 同 `-i` 的值      | `输出目录`                   |
| `-o` | `string` | `result.xlsx`   | `输出目录` 的最终结果 excel   |
