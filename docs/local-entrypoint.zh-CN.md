# 本地入口与环境检查

安装后使用 `molbio`，或等价的 `python -m localmolbio`。当前提供环境检查、前台网页服务和有界队列执行，不安装系统服务、不启动守护进程、不自动领取已有任务。此版本仍需要单独执行 work；尚未达到 V1 普通试用者无需手写 worker 命令的就绪门槛。

## 安装和选择数据位置

在仓库目录运行，示例适用于 macOS/Linux shell：

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
export MOLBIO_DATA_DIR="/absolute/path/to/molbio-data"
molbio doctor
molbio serve --port 8000
```

将数据目录替换为实际希望保存研究文件的位置。已有数据应继续选择原目录，不迁移或复制数据库。全局参数 `--data-dir` / `--minimap2` 应放在 doctor/serve/work 前，也可用对应环境变量配置。环境变量只作用于当前终端；另开终端执行任务时需使用相同数据位置与工具。

`serve` 在前台运行，启动时输出实际绝对数据路径和网址，固定绑定 `127.0.0.1`，没有对外监听选项。Ctrl+C 结束服务。远程工作机的连接方式仍见 [Tailscale 部署](workstation-tailscale.zh-CN.md)。启动前先占用所选端口，端口不可用会在数据库初始化/迁移之前退出。旧的 uvicorn/worker 命令保留兼容，但新入口不再默认为当前目录的 var，防止从不同目录启动时误建另一个库。

安装包包含 Python 应用和前端资源，不包含 minimap2；首次 pip 安装可能需要网络获取依赖。本次只验证 M2 macOS / Python 3.12 的新虚拟环境安装，未声称全新操作系统、Windows 或所有声明版本均已验收。

## doctor 的含义

`doctor` 不创建数据目录、不初始化数据库、不领取任务。检查依赖是否能导入、已安装版本、所选目录或最近存在父目录的写入权限、可用磁盘字节、是否已有数据库，并报告可选比对器状态。权限检查是当时状态，不是未来写入成功或数据正确性的保证，不做磁盘写测试或数据库完整性检测。

基础应用可用时 ready=true；未配置 minimap2 不影响图谱、引物和 Sanger，因此基础检查可成功但 alignment.ready=false。需要 FASTQ 比对时使用严格检查：

```sh
export MOLBIO_MINIMAP2="/absolute/path/to/minimap2"
molbio doctor --require-alignment
```

仅接受已验证的 minimap2 `2.31-r1302`，报告二进制 SHA-256。配置应指向本机受信任程序；版本检查会执行该程序的 `--version`，版本文字不是安全认证。工具缺失/版本不符给出可读说明，严格检查退出码为 2；安装/构建固定工具见 [比对器说明](alignment-evidence-core.zh-CN.md#固定工具与运行范围)。不会自动下载或改用其他工具。

## 执行已排队任务

先启动界面、登记输入并排队，再在另一个已激活相同环境的终端运行：

```sh
molbio --data-dir "/absolute/path/to/molbio-data" --minimap2 "/absolute/path/to/minimap2" work --adapter alignment --max-jobs 5
```

- adapter 必选：alignment 执行已登记 FASTQ 的局部比对；ab1 / fastq 仅核验登记输入的身份，不是 Sanger 分析或全质粒验证。
- max-jobs 默认 1，范围 1–100，计算领取尝试次数。每次输出一行 JSON；遇到 idle 或第一次失败/租约丢失立即停止，不自动耗尽一个失败任务的重试预算。默认只执行一次；运行后回网页刷新。
- 比对工具未通过预检不领取任务；目标目录未初始化数据库时也拒绝执行。成功/空闲返回 0，任务失败返回 1，配置问题返回 2；交互中断返回 130，未完成租约需等到过期再接管。
- 复用原队列的隔离、续租与正式发布规则，没有改变输入限制或科学结论。不得因命令返回 0 或 status=succeeded 把局部比对解释为质粒通过。

## 验收与仍未解决的成本

17 项入口测试检查明确存储、doctor 不写入、依赖错误、固定工具、端口冲突不建库、显式 adapter/范围、有限执行和首错停止。M2 新 Python 虚拟环境只安装本项目 wheel 及应用依赖，实际检查入口来自该环境 site-packages；浏览器测试驱动仍使用原开发环境的 Playwright，应用/服务/队列执行均使用新环境安装包。

通过新入口实际完成合成 GenBank 导入、FASTQ 上传、排队、真实 minimap2 分析、正反向/跨起点/删除审阅及 HTML/JSON 下载，继续检查窄屏与错误恢复。复用宿主已有 Python、Chrome、构建工具与 minimap2，不代表从空白机器完成安装。

仍需处理：用户不应在每次排队后手写 work 命令；固定工具获取/配置也有首次成本。下一步只补此闭环的可控执行与可见错误，不借此扩展多用户或更多分析器。开发者自动化运行不能算“首次用户独立完成”，时间/可靠性/回访仍按 V1 真实验收记录。
