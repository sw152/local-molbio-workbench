# 已登记 FASTQ 的有界比对任务

`registered-fastq-minimap2-v1` 将原始 FASTQ 登记、独立输入快照、真实 minimap2 和持久化队列连接起来。目前提供 Python 接口、CLI 和 [修订范围内分页 API](alignment-api.zh-CN.md)，并已连接 [局部比对审阅界面](alignment-ui.zh-CN.md)。成功表示产生了已校验的局部比对证据，不是生成了共识、变异报告或整质粒通过结论。

## 运行

应用与 worker 使用同一 MOLBIO_DATA_DIR。MOLBIO_MINIMAP2 只能由受信任的本机配置指定，不由请求提供路径或命令。本项目在 M2 执行，固定工具安装方式见 [比对核心](alignment-evidence-core.zh-CN.md)。

```sh
export MOLBIO_MINIMAP2="$PWD/var/tools/minimap2-v2.31/minimap2"
python -m localmolbio.worker enqueue-alignment --revision REVISION_ID --input-id FASTQ_INPUT_ID --data-type ont-high-accuracy --key UNIQUE_REQUEST_KEY
python -m localmolbio.worker once --adapter alignment --worker-id m2-alignment
```

`--data-type` 必填，只允许 ont-noisy / ont-high-accuracy / pacbio-hifi / short-single。多个文件可重复 `--input-id`，最多 100 个不同登记 ID、必须属于相同修订。按 ID 排序固定查询顺序，同一集合不同顺序复用请求键仍返回原任务；类型、输入、参考、工具哈希或重试预算变化后重用同键会冲突。默认 once 仍只领取 AB1 检查；alignment 必须显式选择，不会自动启动后台服务。

提交固定参考身份与哈希/拓扑、各 FASTQ 的原始哈希/大小/压缩/编码/摘要哈希、数据类型及 minimap2 版本/二进制哈希。执行时重新验证。工具缺失或版本不受支持明确拒绝，不切换主机、不改用其他算法。旧任务不会静默使用新工具二进制。

## 输入快照和限额

每次尝试在数据目录的 `alignment-attempts/<job UUID>/<独立 UUID>/` 保存：attempt.json、inputs/ 原始字节副本、sources.json、alignment/ 核心输入与输出、以及完成执行后的 attempt-result.json。文件从受管理普通文件描述符分块复制并同时校验，独立副本不使用硬链接；不覆盖其他尝试。不公开原始路径和 lease token，结果只提供数据目录内相对 artifact_directory。

快照再次完整解析校验，重新计算摘要并与登记摘要核对（仅限额配置字段不作内容比较）。`sources.json` 给出 q0/q1 等本次查询名到登记 input_id、1 基 record_ordinal、规范化序列哈希的映射。重复 FASTQ 标题不会合并记录，不依赖标题猜测配对。只将 DNA 字母大小写统一为大写，其他不支持 IUPAC 明确拒绝；原始字节和质量仍保存在快照内，但质量不参与当前比对。

本轮有意限制为小规模完整输入：所有原始文件合计 32 MiB，解压内容合计 32 MiB；合计 5,000 条、每条 100,000 bp、总计 5,000,000 bp；参考最多 200,000 bp。超限拒绝整项任务，不截前缀，不静默跳过记录。登记接口更大的接收限额不代表可以用此适配器直接分析同样大的数据。大规模 FASTQ 的分块和增量结果处理尚未实现。

快照扫描每约 64 KiB、记录转换每 100 条及复制/比对期间检查 pulse；丢失租约或取消则停止，不继续发布。每个独立目录是逻辑隔离的快照，并不是 OS 强制不可变存储；外部拥有者仍能修改文件，因此解析后和发布前都会核对哈希。

## 结果发布与失败

完成比对后重新核对原始登记文件实际字节、快照、来源映射、参考正文与数据库登记。即使快照已经成功分析，原始文件或参考随后在本次任务中变化，也按保守策略拒绝发布。任一输入失败不会发布部分成功。

attempt-result.json 仅为诊断产物，可能在最终提交前丢失租约；不能凭文件存在判断任务成功。只有 job_queue.finish 接受当前未过期 token 后，analysis_jobs 的 succeeded 与 result_summary_json 才是正式结果。取消不删除既有诊断文件；过期接管和 I/O 重试建立新的独立目录，结果保存 job_id 和 attempt_number 以追溯。尚未实现失败产物自动清理或磁盘配额管理。

科学/输入不一致和工具进程失败为终止性错误；一般输入/产物 I/O 错误在预算内重试，等待下一次显式 once。未知异常记为 unexpected_alignment_error，避免将异常中的文件路径或序列写入通用任务错误。任务成功设置 analysis_performed=true，scope=registered_fastq_local_alignments_only；consensus_performed、whole_reference_verified 仍为 false。

完整结果包含各读段局部匹配、工具与参数来源及快照索引。当前直接保存有界 JSON 结果，分页 API 读取已提交结果，并提供 [配对位置覆盖审阅](alignment-coverage.zh-CN.md)，但尚无质量感知覆盖/差异报告、共识或大规模性能承诺。圆形副本、非穷举多重匹配、MAPQ、N 与 gap 分母等边界沿用 [核心语义](alignment-evidence-core.zh-CN.md)。后续科学算法变更需更新适配器/schema 版本，不能静默改写历史结果。

## 验收

22 项专项测试配置固定 minimap2，成功比对和执行后身份变化等用例真实运行工具；涵盖普通/gzip 多文件和重复标题、独立副本/序号追溯、大小写规范化、圆形跨起点、坏摘要和不支持 IUPAC、整批限额拒绝、原始/快照/来源映射/参考变化、扫描取消、失效后接管、独立尝试目录、I/O 重试、工具身份冲突、错误脱敏以及独立 CLI 进程。需要设置 MOLBIO_MINIMAP2；缺少环境变量时相关集成测试明确跳过，不可声称实际执行。所有样例为合成数据，没有实际仪器或大文件性能结论。
