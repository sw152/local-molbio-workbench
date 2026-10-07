# 局部比对审阅报告

报告让同事在不运行应用时复核一项已发布的 FASTQ 比对任务。它保留参考、输入、工具身份，覆盖统计、完整待复核区间和全部读段候选，不输出质粒通过结论。与 既有 Sanger 报告 是不同的结果类型，当前不合并为跨平台共识。

## 使用

在成功任务的 Alignment review 中选择 Download HTML report 或 Download JSON evidence。HTML 可独立离线打开并使用浏览器打印，包含所有记录的方向、跨起点坐标、paired-block 示意轨、计数、CIGAR、原始 MAPQ、withheld 和来源身份。JSON 额外保留结构化 paired_blocks、全部 raw_hits、限制参数和来源映射哈希，方便程序复核。

**导出整项任务，不仅是当前页或当前区间筛选。**页面下载按钮旁明确说明此范围；报告开头再次显示完整记录数。报告包含文件名标签、ID 和哈希，标签是在导出事务中读取的当前元数据；原始序列、质量数组、FASTQ 文件、私有存储路径和租约凭证不包含在报告中。

HTML 使用内置样式，无脚本、网络字体或外部图像；所有文字值转义，提供限制内容加载的 CSP。浅色样式用于阅读/打印，390px 宽度重排卡片和长哈希；报告配对轨与参考坐标对应，删除保留间隙。页首组成条只表示基数构成，不是空间深度图。浏览器打印样式已检查，但未提供独立 PDF 下载或纸张分页保证。

## API 与来源

`GET /api/revisions/{revision_id}/alignments/{job_id}/report?format=json`，format 可选 json（默认）或 html。响应为 attachment，文件名仅使用任务标识，Cache-Control 为 no-store。

JSON schema 为 `localmolbio.alignment-review-report`，schema_version=1，包含：

- exported_at：UTC 导出时间；job：原任务 ID、创建/完成时间、发布状态与 attempt_number。
- reference：修订 ID、参考 SHA-256、拓扑及长度；inputs：原文件/登记摘要 SHA-256、大小、压缩/质量声明和 label_at_export；tool：名称、版本、二进制哈希。
- evidence：原始结果 schema、scope、测序类型、坐标约定、工具规模/候选限制和各项未评估标记。
- coverage：与覆盖 API 相同的精确基数与记录深度统计，不含全部游程；review_regions：全部 unpaired、ambiguous_only、deletion 区间，不分页、不截断。
- reads：全部原记录来源、状态、候选/paired_blocks/raw_hits 和 withheld；selection 恒为 all_task_records。
- source_verification_at_export=not_performed，limitations 明示解释范围。

只读取修订/适配器相符且 queue.finish 正式发布的 succeeded 结果；未发布/取消返回 409，跨修订或其他适配器返回 404，格式非法返回 422。读取采用单一数据库事务，校验来源映射、参考/工具/输入/attempt 与提交清单一致性、结果 scope 和 false 标记，并重新检查全部候选的 CIGAR/paired_blocks 几何；损坏则拒绝整份导出，不输出部分成功报告。

导出不重跑比对、不选择较新的任务、不再读取原始磁盘文件；原文件事后变化不会改变这份已提交历史证据。标签可以变化，因此不作为不可变身份。当前实现仍读取和渲染整份有界结果，不是海量数据增量报告；较大任务的 HTML 可能很长，没有规模延迟承诺。

## 科学边界

坐标为 0 基半开，原文件内记录序号为 1 基。配对包含错配/N，删除不计配对，可以与其他记录的配对重叠；一个记录的重叠候选取并集，不冒充独立分子。单候选不证明唯一，withheld-only 不推断位置；候选非穷举、圆形 MAPQ 未校准、indel 未归一化。无质量感知共识、读对推断、正式变异调用或全质粒验证。

## 验收

报告专项使用真实固定 minimap2，覆盖完整多记录输出、反向/删除/无命中、身份与来源、未知 schema/scope/坐标和坏几何拒绝、参数/修订/适配器/取消边界、HTML 转义、下载安全头、私有路径/原始序列不公开，以及原文件事后变化不触发重新执行。

M2 Chrome 从筛选后第二页实际下载两种格式，验证仍包含全部 8 个合成读段和未分页区间；独立打开 HTML 检查 withheld、删除间隙、桌面/390px/打印媒体样式及零 HTTP(S) 资源请求。下载失败保留当前证据，关闭/切换修订后迟到导出不会触发下载。没有仪器数据、用户试用或省时对照结论。
