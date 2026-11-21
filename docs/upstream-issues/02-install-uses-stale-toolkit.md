# Bug: `install.py --runtime-only` silently builds the wheel from the *current checkout's* toolkit

## 摘要
`scripts/install.py --runtime-only --toolkit-path <path>` 会用 `<path>` 里的 toolkit 打包
`cdxml-toolkit-community`。如果不显式给 `--toolkit-path`，它用的是**当前 checkout** 的
`toolkit/`。当本地 checkout 停留在**落后于 main** 的分支（例如 `claude/figure-coverage`）时，
装进 venv 的 toolkit 就缺少 main 上已有的能力（`fragment_colors` / `row_align` 等），
而这**不会给出任何警告**。

## 影响
- 用户看到的是 benchmark / 绘图报 `Unknown brief keys`，误以为是语料问题或版本问题
- 实际是"安装源比语料旧"
- 我们为此在三个错误方向上排查了两轮（venv 损坏、pypi SSL、运行时版本）

## 复现
1. `git checkout <落后 main 的分支>`
2. `python scripts/install.py --runtime-only --extra publication --apply`
3. 用该 venv 跑一个带新字段的 brief → `Unknown brief keys`

## 建议
1. 安装前检查：若 `HEAD` 不是 `origin/main` 的祖先（或 `HEAD` 落后），打印警告并询问是否改用 main
2. 或在输出的 JSON 里记录**实际使用的 toolkit commit**（`git -C toolkit rev-parse --short HEAD`），
   便于事后核对"装的和语料是不是同一代"
3. 或让 `--runtime-only` 默认用 `origin/main` 的 toolkit（而非当前 checkout）

## 顺带
`pyvenv.cfg` 的 `home` 指向的基础解释器被删除后，venv 的 `python.exe` 会报
`No Python at '...'`，此时 `install.py` 只返回 exit 103，没有提示"重建 venv"。建议在
`install_runtime.py` 里捕获这个错误并给出 `--force-venv` 之类的修复建议。
