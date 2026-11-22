# `D:\Work\forks\chemdraw-skill` —— 第三方 ChemDraw 技能的 fork（开发参考）

> **一句话**：fork 自 `ZiChenWang114514/chemdraw-skill`（MIT 开源）—— **一个把 ChemDraw 变成 AI agent 可调用技能的成熟项目**。
>
> **用途判断：它更可能是「参考/借鉴/上游」，不是「我自己在维护的分支」。**

**测量时间**：2026-11-08 · **111.5 MB / 1,732 文件 / 345 子目录** · 今天有改动

---

## 一、它是谁的、做什么的

**⇒ 原始项目**：`github.com/ZiChenWang114514/chemdraw-skill`（作者 ZiChenWang，MIT License）

**⇒ 它的定位（README 原文）：**
> "Editable ChemDraw chemistry for Claude Code, Codex and any AI coding agent"
> 画结构、重构论文图件、从表格建 SAR 面板、把活的 ChemDraw 对象塞进 Word/PowerPoint。
> 同一个 Skill 还能称量反应、写表征段落、在反应中映射原子。

**⇒ 支持的 agent**：Claude Code · Codex · Gemini CLI · GitHub Copilot · Cursor · OpenCode · **任何能跑 shell 的 agent**

---

## 二、它的核心理念（值得学的四点）

| 理念 | 说明 |
|---|---|
| **结构要有出处** | 结构来自数据库 / OPSIN / ChemScript / DECIMER —— **agent 不许凭记忆敲 SMILES** |
| **输出可编辑** | 结果是 CDXML，**在 ChemDraw 里还能改原子/键/文字** |
| **一条命令** | `cdx.py` 一个入口跑所有工具，MCP 是可选的 |
| **存完要回读校验** | ChemDraw 渲染保存后，**把化学读回来跟输入比对** |

**⇒ 这四条和本机 `web-standard` 的"先量后报 / 验收站在用户判据上"是同一套思路。**

---

## 三、和本机自有技能的关系

**⇒ 本机 `D:\Work\tools\skills\catalog\` 里有三个 ChemDraw 相关技能：**

| 本机技能 | 和这个 fork 的关系 |
|---|---|
| `chemdraw` | 🟡 **可能借鉴了它的思路** |
| `chemdraw-automation` | 🟡 同上 |
| `chemdraw-workbench` | 🟡 同上（"统一入口"这个思路很像它的 `cdx.py`） |

**⇒ 还有 `D:\Work\tools\live-chemdraw-mcp\`（2.4 MB）—— 本机自己的 ChemDraw MCP 驱动。**

**⇒ 所以本机是**两条腿**：**
```
自研：live-chemdraw-mcp（COM 驱动）+ 三个自有技能
参考：forks\chemdraw-skill（MIT 上游，学它的架构）
```

---

## 四、里面有什么

```
chemdraw-skill\
  ├── README.md / README.zh-cn.md     919 行 / 中文版
  ├── AGENTS.md                       给 agent 的规则
  ├── LICENSE                         MIT
  ├── .github\                        工作流（build / pages / release / site / validate）
  ├── skill\                          技能本体
  │     └── chemdraw\references\      参考资料（含 agent-integration.md）
  ├── assets\readme\                  919 行 README 用的图表
  │     ├── gallery\                  16 组对照图（原图 vs ChemDraw 复刻）
  │     ├── engines\                  引擎示例（forward/lcms/nmr/ocsr/pdf/retro）
  │     ├── markush\                  Markush 结构示例
  │     ├── mechanism\                反应机理示例
  │     └── paper-replica\            论文图件复刻示例
  └── docs\guide.md
```

**⇒ 111.5 MB 里大部分是 `assets\readme\` 的图（16 组 gallery 对照图 + 引擎示例）。**

---

## 五、能不能动

| 项 | 判定 |
|---|---|
| 能不能搬 | 🟡 可以（纯文件 + `.git`） |
| 能不能删 | 🟡 **看是否还需要参考** —— 它随时能重新 `git clone` |
| 有没有本地改动 | ❓ 待查（`git status` 看） |
| 为什么放在 `forks\` | 因为它是**第三方项目**，和自己的代码分开 |

**⇒ 建议：**
1. **如果本地有改动（自己加了东西）** → 留着，可能有用
2. **如果只是原样 fork** → 可删（`git clone` 能重得），或移到 `Archive\`

**⇒ 查本地改动：**
```
git -C "D:\Work\forks\chemdraw-skill" status --porcelain
git -C "D:\Work\forks\chemdraw-skill" log --oneline -5
```

---

## 六、这个目录存在的意义

**⇒ `D:\Work\forks\` 这个目录名本身就说明了用途：**

```
forks\  = 「别人写的、我 clone/fork 下来的」
projects\ = 「我自己做的」
tools\  = 「工具」（自研 + 第三方装的）
```

**⇒ 所以 `forks\` 是**工作资产的第四类：外部代码的本地副本**。**

**⇒ 目前只有 `chemdraw-skill` 一个。**

---

*本 README 由 2026-11-08 整理工程建立。*
