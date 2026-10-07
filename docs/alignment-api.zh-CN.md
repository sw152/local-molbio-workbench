# 修订范围内的局部比对任务 API

此接口连接已登记 FASTQ、有界 minimap2 worker 与可分页的结果审阅。已连接 [比对提交/审阅 UI](alignment-ui.zh-CN.md)；调用提交接口不会启动 worker，仍需显式运行 `python -m localmolbio.worker once --adapter alignment`。应用与 worker 使用同一数据目录和固定的 MOLBIO_MINIMAP2 配置，部署、类型和规模边界见 [登记 FASTQ 比对任务](registered-fastq-alignment.zh-CN.md)。

## 路由与请求

前缀：`/api/revisions/{revision_id}/alignments`。

| 方法与后缀 | 行为 |
|---|---|
| POST 空后缀 | 提交固定比对适配器，返回 201 任务摘要；同键同输入返回原任务 |
| GET 空后缀 | 分页任务列表，默认 limit=5、offset=0，按创建时间和 ID 降序 |
| GET /{job_id} | 身份清单、当前文件名标签、结果摘要、分页尝试历史（默认 5 条，尝试序号降序） |
| GET /{job_id}/reads | 已正式发布的逐读段证据，可按参考区间/关系过滤，默认 20 条，按原查询索引升序 |
| GET /{job_id}/report | 全任务 HTML/JSON 离线审阅报告，只允许已发布结果 |
| POST /{job_id}/cancel | 取消 queued/running；重复取消返回相同状态；已成功/失败返回 409 |

提交 JSON 必须显式给出：

```json
{
  "input_ids": ["REGISTERED_FASTQ_ID"],
  "data_type": "ont-high-accuracy",
  "idempotency_key": "CALLER_GENERATED_KEY",
  "max_attempts": 3
}
```

input_ids 为同一修订的 1–100 个不同登记 ID。data_type 只能为 ont-noisy、ont-high-accuracy、pacbio-hifi、short-single。max_attempts 可省略，严格整数 1–10；不接受布尔值、任意 adapter、parameters、binary 或路径。相同键但输入/工具/类型/重试预算发生变化返回 409。响应丢失应以原键原请求重试，不应重新生成键。

分页统一返回 items、total、limit、offset、has_more；limit 范围 1–100，offset 非负。列表与详情读取同一个数据库事务快照。分页仍使用 offset，不承诺跨多次请求遇到新增任务时的快照一致性；前端刷新需重置页码。修订不存在或任务不属于该修订/比对适配器返回 404，非法请求返回 422。

## 正式结果与来源映射

只有任务状态为 succeeded 才提供 result 摘要和 reads。排队、运行、失败、取消的 result 为 null，reads 返回 409 alignment_evidence_not_published。接口不扫描尝试目录，不会把未成功提交的 attempt-result.json 当成可见结果。取消通过既有队列租约机制阻止迟到尝试发布。

详情给出参考身份、输入原始字节/登记摘要哈希、工具版本/二进制哈希和数据类型。文件名标签取当前登记，不作为不可变内容身份。公开尝试历史包含序号、状态、时间和稳定错误码；不包含 lease token、内部尝试 ID、worker 标识、请求键、存储路径或诊断目录。列表不携带结果数组；详情结果也不携带全部 reads/sources/query hashes。接口不提供原始文件下载。

每条 reads 项保留局部命中、原始 query 坐标、参考区间、方向、eqx CIGAR、paired_blocks、精确/错配/N/gap 计数、原始命中信息、withheld 与 requires_review。source 使用登记 input_id、文件内 **1 基** record_ordinal、本次 qN 和大写序列 SHA-256；不会以重复 FASTQ 标题合并记录。返回前核对数组长度、索引/查询名、输入归属、文件内连续序号和序列哈希；失配返回 409 alignment_source_mapping_invalid。未知 schema 版本也拒绝读取。

sources_sha256 是 worker 验证过的原始 sources.json **字节**哈希。数据库会规范化 JSON 键顺序，因此 API 不把重序列化 JSON 的哈希冒充原文件复核。API 不重读磁盘输入或重新执行比对；其内容是已提交历史快照，不表示登记文件当前仍未被外部修改。

## 科学与规模边界

成功的 scope 为 registered_fastq_local_alignments_only，analysis_performed=true；whole_reference_verified、consensus_performed、base_quality_used、pairing_used 均为 false。坐标为 0 基半开区间，方向相对原始读段。MAPQ 为原始工具值，圆形双拷贝未重新校准；单命中不表示唯一，候选非穷举，局部高 identity 不表示整质粒通过。未报告命中与被 withheld 的命中需审阅，不自动解释为生物学失败。删除不属于 paired_blocks 的覆盖。

分页限制返回的读段数，并未实现数据库增量读取：当前会解析完整的有界 JSON 结果并检查来源映射，再切页；单条复杂命中也可能较大。这不是海量 FASTQ 的内存/延迟或响应字节上限保证。新增 [配对位置总览与区间 API](alignment-coverage.zh-CN.md)，当前没有质量感知共识、变异报告、账户/项目权限；修订隔离不等同于多人授权。

## 验收

11 项 API 专项测试在 M2 配置真实固定 minimap2，覆盖普通/gzip 多文件、重复标题、正反向精确坐标、未报告命中、来源映射、任务/尝试/读段分页、幂等冲突、跨修订/适配器隔离、严格参数、凭证与路径不公开、取消后禁止发布和损坏来源/版本拒绝。没有实际仪器输入或大规模性能结论。

区间筛选的 start/end/relation、先筛选后分页和 region_evidence 语义见 [区间进入原始记录](alignment-coverage.zh-CN.md#从区间进入原始记录)。不筛选时原读段字段不变，响应增加 region=null 与 unfiltered_total。

完整任务导出及安全/科学边界见 [局部比对审阅报告](alignment-report.zh-CN.md)。
