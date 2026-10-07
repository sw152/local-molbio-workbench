# 已登记 FASTQ 的队列身份复核

`registered-fastq-input-check-v1` 是受限的单机输入检查适配器。它复核已经登记的 FASTQ 原始字节与参考、质量编码声明等身份，不解析新的任意路径，不进行比对、配对、共识或变异分析，也不重新计算登记摘要。

## 显式提交和执行

应用与 worker 使用同一 `MOLBIO_DATA_DIR`，本项目命令在 M2 上执行。ID 来自该实例的 FASTQ 登记 API；不是文件路径。多文件重复提供 `--input-id`，要求 1–100 个不同 ID，全部属于同一参考修订。

```sh
python -m localmolbio.worker enqueue-fastq-check --revision REVISION_ID --input-id FASTQ_INPUT_ID --key UNIQUE_REQUEST_KEY
python -m localmolbio.worker once --adapter fastq --worker-id m2-fastq-check --lease-seconds 60
```

`once` 默认继续只领取 AB1 输入检查。`--adapter fastq` 仅领取 FASTQ 输入检查，最多一项后退出；没有常驻服务，也不会自动同时启动两类 worker。未知 adapter 拒绝执行。退出码和重试语义沿用 [单次 worker](input-check-worker.zh-CN.md)：成功/idle 为 0，本次失败或租约丢失为 1，非法请求为 2。`--max-attempts` 为 1–10，默认 3。

同一输入集合排序后生成清单，相同请求键返回原任务，不代表重跑；更改输入或声明后重用同键会冲突。新一轮检查应使用新键。本阶段提供 Python 接口和 CLI，尚未提供专用 FASTQ 检查提交/详情 UI；现有 AB1 检查 API/UI 仍限定 AB1。通用 `/api/jobs` 只可查看顶层状态。内部 analysis_jobs 沿用原有 job_kind，必须以 adapter、scope 和 false 标志判断实际执行范围。

## 被固定和复核的身份

- 参考修订 ID、实际序列 SHA-256 及拓扑；提交时和执行时检查数据库序列正文与哈希一致。
- FASTQ 登记 ID、原始字节 SHA-256、字节大小、none/gzip 压缩类型、Phred+33 声明，以及规范化登记摘要 JSON 的 SHA-256。JSON 空白与键顺序不影响摘要身份；统计值或声明变化会影响身份。
- 执行前及完成后重新读取登记，拒绝归属改变、哈希/编码/压缩/摘要变更。仅接受已支持的 `bounded-fastq-v1` 登记摘要及明确的单文件未配对范围。
- 实际文件大小须等于登记大小。受管理 FASTQ 目录内的普通文件按 1 MiB 分块计算 SHA-256；gzip 核对压缩原始字节，不对不同 gzip 包装做等价归并。
- 与 AB1 共用受管理文件读取检查：拒绝目录、FIFO、符号链接和目录外路径；检查打开前后的 inode、设备、大小、mtime、ctime。按块续租，结束后再次检查登记；当前 token 或租约失效则不得发布。

成功结果为 `scope: registered_fastq_integrity_only`，明确 `summary_recomputed: false`、`analysis_performed: false`、`whole_reference_verified: false`。仅包含参考身份和各输入 ID/哈希/大小/声明/摘要哈希，不复制文件、序列、读段标题或存储路径。登记摘要哈希只能证明所固定摘要没有改变，不能独立证明统计计算正确。

终止性失败不重试：文件缺失、大小/哈希改变、元数据不一致或不支持、错误归属/路径等。一般 I/O 错误可在预算内重新领取；未知异常只记录固定错误代码。任一文件失败时不发布整个任务的部分成功。取消或过期撤销结果发布权，不表示强制终止底层阻塞 I/O。

检查是逐文件的时间点身份复核，不是跨文件原子快照，也不能保证成功后文件仍未变化。未来分析适配器必须在执行时再次复核或使用独立不可变快照；不能只相信旧检查结果。当前打开方式依赖 POSIX 标志，已在 M2 macOS 验证，没有 Windows worker 支持。未对恶意并发文件系统操作者或大文件性能作出保证。

## 验收

26 项 FASTQ worker 专项测试覆盖普通/gzip 原始字节、完整元数据身份、多文件排序幂等和末项损坏无部分成功、字节/大小/参考/归属/声明/摘要变化、目录边界与非普通文件、读取中变化、JSON 规范化、重试耗尽、取消和过期接管、AB1/FASTQ 领取隔离、未知错误脱敏，以及独立 CLI 进程。输入都是本地合成数据；不代表真实仪器 FASTQ 或大规模测序分析验收。

## 后续比对适配器

输入身份检查的语义不变。真正的有界局部比对使用另一个 adapter，通过 `enqueue-alignment` / `once --adapter alignment` 显式选择，见 [已登记 FASTQ 比对](registered-fastq-alignment.zh-CN.md)。身份检查成功不会自动触发比对。
