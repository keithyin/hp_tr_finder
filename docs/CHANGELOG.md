# Changelog

本项目从 `0.5.0` 起记录变更。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。

## [0.6.0] - 2026-11-21

### Added

- 新增嵌套注释去重，默认关闭，由 flag 启用：CLI 增加 `--dedupNested`，
  lib 侧 `all_seq_hp_tr_finder` / `single_seq_hp_tr_finder` 增加尾参 `dedup_nested: bool`。
- 嵌套去重规则（`dedup_nested_hits`，在旋转去重之后运行）：外层命中完全包含内层命中时
  按拷贝数 `copies = (end - start) / motif 长度` 裁决——内层 copies > 外层 → 仅保留内层；
  内层 copies ≤ 外层 → 保留外层。详见 `docs/known-issues.md` 修复说明。
- 嵌套去重的单元测试与 CLI flag 端到端测试。

### Changed

- `all_seq_hp_tr_finder` / `single_seq_hp_tr_finder` 签名增加 `dedup_nested: bool` 参数
  （破坏性变更）。

## [0.5.0] - 2026-11-20

### Added

- `called` 场景新增区域扩展参数：`--extLeft` / `--extRight` 按固定长度扩展命中区间，
  `--extRatio` 按区间长度比例扩展（并提供 `--extMax` 上限），扩展后重叠或相邻的区间会合并输出。
- 公开 `expand_regions` 与 `merge_region` 工具函数，用于区间扩展与重叠合并。
- 公开 `simplify_motif` 与 `is_primitive_motif`：`motif` 模式会把非最简 motif 化简为最简形式并按拷贝数
  等比放大最小重复数（如 `GTGT-3` 等价于 `GT-6`）。
- `motif` / `unit` 参数校验：空 motif、`min repeats == 0` 等非法输入会给出明确错误信息。
- 补充 `build_finder_regs`、motif 化简、`called` 场景扩区合并等单元测试。

### Changed

- FASTA 解析改用 `noodles-fasta`；重复的序列 ID 按拼接方式合并，并保持首次出现顺序。
- 生成 motif 时过滤非最简 motif（含同聚物，如 `AAAA`、`GTGT`），避免同一串联重复被长短 unit 重复报告，
  同时消除重复正则导致的重复命中。
- 旋转等价命中去重（原 `dedup_overlaps` 内部实现改为 `dedup_rotation_hits`），相邻区间不再视为重叠。
- GFF3 输出的起始坐标改为 1-based（结束坐标保持与 tr-finder 一致的 0-based 风格）。
- `Cargo.toml` 中的项目描述更新为 “Find homopolymers and tandem repeats in a FASTA file.”。
