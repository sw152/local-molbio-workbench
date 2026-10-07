# 单机持久化队列基础

本模块只提供任务生命周期原语，不启动后台进程，不执行命令，也没有注册测序分析器。独立的输入完整性检查 worker 已接入，见 [使用说明](input-check-worker.zh-CN.md)。当前 UI 已支持同一修订多个已上传 AB1 的输入检查提交，见 [任务界面](input-check-ui.zh-CN.md)。阶段 3 仍在进行中，队列完成不表示全质粒测序功能已完成。

数据库以新增 `queued_jobs` 和 `job_attempts` 的方式扩展现有 `analysis_jobs`。旧的同步引物设计及 Sanger 任务不进入此队列，迁移不改写它们。使用现有 `MOLBIO_DATA_DIR`；调用前执行 `database.initialise()`。

## 当前接口

- `enqueue(revision_id, idempotency_key, adapter, input_manifest, parameters=None, max_attempts=3)`：将参考修订 ID、哈希、拓扑及调用方输入清单、适配器版本、参数固化保存。全局唯一请求键与规范 JSON 摘要绑定；相同请求返回原任务，包括已经结束的任务。不同内容复用同键抛出 `QueueConflict`。有意重新运行必须使用新请求键。
- `claim(worker_id, adapter, lease_seconds=60)`：原子恢复到期租约，再按创建时间领取一项兼容任务。返回新的尝试号、租约 token、到期时间和输入清单。没有匹配任务返回 `None`。每次失败/过期后领取都会产生新尝试与新 token。
- `renew(job_id, lease_token, lease_seconds=60)`：仅当前未到期租约可以续期。等于到期时间已经失效。
- `finish(job_id, lease_token, result)`：持有有效 token 才能把结果摘要与成功状态一起保存。任务成功仅表示适配器完成，结果不得隐含全质粒生物学通过。
- `fail(job_id, lease_token, error, retryable=False)`：保留尝试错误；明确可重试且尚未达到上限时重排，否则终止为失败。
- `recover_expired()`：将到期尝试标为 expired；预算内重排，耗尽则失败。恢复是幂等的；也由 claim 自动触发，不依赖应用启动时的一次扫描。
- `cancel(job_id)`：queued/running 可以取消，重复取消无副作用；已完成任务不改写。取消会撤销租约，但尚未接入进程终止功能。

## 一致性边界

每次状态转换使用 SQLite `BEGIN IMMEDIATE`，状态、尝试记录和审计在同一事务中提交。计算在事务外完成。结果提交同时检查任务 running、当前 token 和租约未到期，因此失联 worker 即使恢复也无法覆盖已接管的任务。

这是至少一次执行及有租约约束的结果发布，不是外部副作用恰好一次。未来 worker 必须把文件写入每次 attempt 独立目录，成功时再发布不可变产物引用；取消/过期后的迟到 worker 不能写共享的最终输出路径。输入文件哈希复核、参考哈希复核、产物存在性验证和进程控制由分析器集成层实施，队列目前只保存清单，不读取或保证清单指向的文件。

租约使用同一主机的 UTC 墙上时钟并持久化到期值，支持进程重启。系统时钟大幅调整可能提前或延后恢复。此实现限定单机 SQLite；没有宣称多主机租约、账户权限、公网 worker API 或 PostgreSQL 支持。重试没有退避调度；适配器未安装时任务保留 queued，不能假装已完成。调用方应只注册明确支持的适配器，错误文本不应包含凭据。

## 已验证

14 项队列专项回归：6 个并发提交者只产生一项任务和一次入队审计；6 个并发领取者只有一人获得租约；重复请求/冲突请求、续期和精确到期边界、旧 token 的成功/失败/续期拒绝、取消、尝试上限、永久失败、初始化幂等、非法 JSON 输入/结果回滚。单独新 Python 进程重新打开同一数据库，接管已到期任务，验证尝试历史保留与旧 worker 无法发布结果。

全量 238 项 Python、11 项 Node 测试通过。M2 Chrome 合成 AB1 上传、峰图、历史、整组导出、错误恢复和窄屏流程通过；本轮不重复使用公开 AB1，日志明确 `public_ab1_checked: false`。截图人工复核，未增加虚构的批量运行 UI。

后续进展：已加入受限的输入校验适配器和显式单次 worker 命令。任务状态/尝试历史 API 与输入检查批量界面也已接入。长短读比对、共识及变异分析另行实现与验收。
