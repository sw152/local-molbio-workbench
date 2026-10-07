# 显式开启前台 FASTQ 分析模式

此模式让使用者启动一次服务后，在网页登记文件、排队、刷新审阅结果，不再逐次另开终端执行 work。当前仅处理既有的受限 FASTQ 局部比对，未新增分析器、共识或整质粒结论。

## 开启和退出

```sh
molbio --data-dir "/absolute/path/to/molbio-data" \
  --minimap2 "/absolute/path/to/minimap2" \
  serve --port 8000 --with-alignment-worker
```

必须显式带 `--with-alignment-worker`；默认 serve 不启动执行进程。参数意味着授权本次前台服务处理**所选数据目录中全部 queued 的 FASTQ 比对任务，包括启动前已排队的任务**，不是仅限当前浏览器或当前修订。每次串行执行一个任务，所有输入/工具/规模和租约检查仍生效。不会领取 AB1 或 FASTQ 身份核验适配器任务。

先用 doctor --require-alignment 排查固定工具。启动检查依赖、可写路径、工具、端口；拿到同数据目录的托管进程锁，完成数据库初始化与修订准备之后，才启动执行进程。没有 HTTP 启用/命令执行接口，不安装系统服务，也不自动下载或更换比对器。当前托管模式使用 POSIX 文件锁和进程组，M2 macOS 已验收，Windows 不支持此模式。

Ctrl+C / 正常终止服务时停止领取新任务，给当前尝试最多 70 秒完成；超时终止本次 worker 及其比对器进程组。强制终止不伪造失败/成功结果：尚未发布的任务保留 running 租约，过期后由既有队列规则接管。进程异常终止时，父进程监视器清理残余比对器；父进程突然死亡时，worker 在轮询/续租点检测并中止、记录可重试的父进程丢失。操作系统挂起或不可中断 I/O 仍需依靠操作系统和租约恢复，不承诺瞬间结束。

同一数据目录最多一个由此入口托管的执行进程；锁会随父进程退出释放，保留的 `.alignment-worker.lock` 空文件不表示仍在运行，无需删除。文件锁不阻止既有显式手动 worker；这些 worker 仍受数据库任务租约约束。不要同时用多个方式执行同一批任务来做性能对照。

## 失败与页面状态

发生第一次 attempt_failed 或 lease_lost 后，本次托管执行停止，服务继续供审阅/取消/下载。后续排队保持 queued；不会自动反复消耗重试预算。先查看任务错误、处理原因，再停止并重新以分析模式启动服务。终止失败任务不会自动重跑；可重试任务按队列预算在后续启动时重新领取。

只读 `GET /api/runtime` 返回 alignment_worker：enabled、state、completed_tasks、adapter；开启时另有 failure_policy=stop_until_service_restart。不暴露 PID、worker 身份、文件系统路径、租约 token 或输入名称。状态可为 starting、idle、running、stopping、stopped、failed；未启用是 disabled。completed_tasks 为本次托管进程完成数，不是当前修订或生物学验证成功数。

网页 Local analysis mode 显示是否启用、正在执行、退出或因失败停止。点击 Refresh alignment tasks 同时刷新任务与执行状态；请求失败显示 unavailable，不能把旧的 idle 当作最新状态。没有自动状态轮询，页面不提供重启进程按钮。

## 验收和限制

生命周期专项包含真实 minimap2 正常发布/适配器隔离、重复托管锁、首错停止、当前任务正常完成后退出、强制杀组不发布、父进程丢失、worker 崩溃清理比对器、进程启动失败释放锁、数据库初始化失败不启动，以及缺工具拒绝。暂停/崩溃场景使用明确的进程控制替身，不冒充生物学准确性测试。

专用 M2 Chrome 流程从网页上传/排队出发，无任何单独 worker 命令，实际跑真实 minimap2、正反向/跨起点审阅和报告下载；检查状态请求失败恢复、真实不支持输入导致停机、后续任务保持 queued、桌面/390px 显示。默认手动模式原有完整浏览器流程另行回归。

这解决了每次手写 worker 命令的成本，但第一次安装 Python、选择存储和配置 minimap2 仍需按文档操作；自动化测试不能代替普通研究者首次独立完成，也没有证明比 Benchling 更省时或更愿意复用。接下来按 V1 冻结任务/数据和记录表，而不是扩展多人/AI 功能。
