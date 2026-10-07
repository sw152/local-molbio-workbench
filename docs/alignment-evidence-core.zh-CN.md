# 受限 minimap2 比对证据核心

本模块实际执行 minimap2，读取 PAF，并依据输入序列逐碱基校验 CIGAR 后产生可追溯的局部比对证据。它是后续 FASTQ 分析适配器的内部核心，目前还没有接入登记输入快照、持久化分析队列、提交 API 或 UI。输入为显式提供的参考字符串与读段字符串列表，不接受随意的命令行参数，不在安装包中捆绑第三方可执行程序。

## 固定工具与运行范围

M2 上编译的版本为 `2.31-r1302`，官方 tag `v2.31` 的源提交 `3c28777e7e2dcc90f825de1b9f17a89cca7d4452`。源码/二进制留在被忽略的 `var/tools/minimap2-v2.31`，未提交到本仓库。可在 M2 项目目录复现：

```sh
git clone --depth 1 --branch v2.31 https://github.com/lh3/minimap2.git var/tools/minimap2-v2.31
make -C var/tools/minimap2-v2.31 -j4 aarch64=1
MOLBIO_MINIMAP2="$PWD/var/tools/minimap2-v2.31/minimap2" python -m pytest -q tests/test_alignment_evidence.py
```

版本字符串必须匹配已验收版本；运行前后核对二进制 SHA-256，并将其与参数、规模上限、超时和坐标约定写入结果。物化的参考/读段 FASTA 在执行前后也核对哈希，变化时不得以空 PAF 保存“没有匹配”的结果。版本字符串本身不是安全认证，binary 参数只能来自本机受信任配置。官方依据：[发布版](https://github.com/lh3/minimap2/releases/tag/v2.31)、[参数及 PAF 手册](https://lh3.github.io/minimap2/minimap2.html)、[固定版本计数实现](https://github.com/lh3/minimap2/blob/v2.31/align.c)。

调用 `localmolbio.alignment_evidence.align(reference, reads, topology=..., data_type=..., binary=..., output_dir=...)`。output_dir 必须为新目录，保留参考/读段 FASTA、原始 PAF、工具日志及成功时的 result.json。输出和日志含研究内容，应只留本地受控数据目录；失败保留诊断文件，不写成功结果，不自动删除既有路径。读段仅使用 q0、q1 等本次调用顺序编号；连接 FASTQ 登记时还需建立原始输入文件/记录序号映射。

数据类型必须显式选择：

| data_type | preset | 本模块范围 |
|---|---|---|
| ont-noisy | map-ont | Nanopore 长读预设 |
| ont-high-accuracy | lr:hq | 高准确度 Nanopore 长读预设 |
| pacbio-hifi | map-hifi | PacBio HiFi 预设 |
| short-single | sr | 短读，强制不使用配对 |

工具命令固定使用 `-c --eqx`，保留 primary/secondary，secondary 数限制为 20，分数比阈值 0.5，seed=11，线程参数 2、batch 5m。预设后显式 `--frag=no --secondary=yes`，不依赖短读预设隐含的配对或次级匹配关闭。并非完整候选穷举，输出一条命中也不能证明唯一定位。输入以 FASTA 传入，质量不参与比对或统计，明确 `base_quality_used: false`；不会凭质量摘要推断数据类型。

首版上限：参考 200,000 bp，最多 5,000 个读段，每段 100,000 bp，总读段 5,000,000 bp；仅支持大写 A/C/G/T/N，其他 IUPAC 字符明确拒绝，不默默转换。拓扑必须明确为 linear/circular。进程默认 60 秒，允许 1–300 秒；轮询时及退出后检查 PAF 64 MiB、日志 8 MiB，超限拒绝结果。轮询检查不是硬磁盘配额，输出可能短暂超过阈值；不承诺 OS 级内存上限。超时和 pulse 回调失败会终止并回收本次直接子进程。没有后台任务、自动重试或数据库发布权；接入队列时仍需 attempt 隔离和最终 token 校验。

## 科学边界

- PAF 必须包含受支持的 primary/secondary 类型和逐碱基 `= X I D` CIGAR。拒绝近似 M-only 输出、剪接/倒位类型、非法坐标、重复标签、错误参考或 query 身份及不一致的 CIGAR/计数；任何不一致拒绝整个结果。
- 所有坐标为 0 基半开区间；反向命中 query_start/query_end 仍使用原始读段方向，paired_blocks 带 query_step=-1。线性参考末端保持 length，不取模为 0。paired_blocks 只表达实际配对碱基，不将 deletion 视为读段覆盖。
- `local_exact_identity = A/C/G/T 精确配对数 / 全部比对列数`。分母包含错配、N 配对、插入和缺失。minimap2 2.31 的 PAF 第 11 列排除模糊碱基，含 gap 中的模糊碱基；模块用序列正文和 `nn` 标签交叉核对，另外保存 reported_paf_block_length，不能直接把 PAF 两列相除冒充本指标。N 不计精确匹配，也不是通配符。
- 圆形参考以双拷贝运行，归一到原长度；跨起点区间分为末段和起始段。query 区间、方向、起点模长、跨度和 CIGAR 全部相同的双拷贝命中合并，但保留 raw_hits。未做 gap 等价归一化、部分命中归并或跨 supplementary 链组装。
- 双拷贝会改变原始 MAPQ，模块不重新校准，不根据其设定置信度。MAPQ=255 保存为 null。圆形副本边缘可能报告额外局部命中，保留为多条证据；真正重复位置也保留。跨越超过一圈的整条命中 withheld，可能仍存在其他局部命中，不支持 concatemer 解释。
- `no_alignment_reported` 只代表此参数下工具未输出受支持命中。高局部 identity 与 aligned_query_fraction 分开保存；不把局部一致解释为完整参考验证。不产生 consensus、逐位点质量、变异结论、覆盖总览或通过判定，whole_reference_verified 始终为 false。

## 验收与下一步

专项测试分为格式/坐标校验和真正的 minimap2 集成。未设置 MOLBIO_MINIMAP2 时集成用例会明确跳过，不能据此声称工具已运行。M2 验收必须设置该变量并检查无跳过。所有输入为随机种子固定的合成序列，不是仪器性能或实际研究数据验收。

已包含四种显式预设、正反向/线性末端、未报告命中、圆形跨起点和副本去重、真实重复区域、多重命中、已知错配/插入/缺失/N、局部 overhang，以及独立解析的 N-gap 分母、非法输出、超时、输出超限和取消。进程控制失败测试使用短小替身程序，科学比对用例实际运行上述工具，两者不混称。

下一步把该核心接到已登记 FASTQ 的有限规模快照和受租约约束的队列执行：固定文件/记录序号、参考与参数，避免分析期间文件变化，保留 attempt 独立产物并通过 token 才发布。超出核心限制的输入必须明确拒绝，不截取前若干条后冒充完整分析。之后再扩展覆盖、差异与审阅界面。
