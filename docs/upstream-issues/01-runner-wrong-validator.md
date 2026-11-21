# Bug: `run_scope_benchmark.py --corpus journal` validates every brief with the scheme validator

## 摘要
`scripts/run_scope_benchmark.py` 的 journal 路径一律调用 `figure_draw.draw_figure(brief, ...)`，
而 journal 语料里的 brief 是 `kind: "scope"`，带 `row_align` / `fragment_colors` / `columns` 等
scope 专有键。`figure_draw._validate()` 用的是 scheme 的 `BRIEF_KEYS`，于是全部 26 例都判为
`Unknown brief keys`，journal 基准 0/26。

## 复现
```bash
python scripts/run_scope_benchmark.py --corpus journal --render --output-dir /tmp/j --limit 1
```
输出：
```
[01/01] FAIL  bjoc-articles-22-103-scheme1  pages=0  readback=None
    errors: ["Unknown brief keys: ['fragment_colors', 'row_align']",
             "header substrates[1] has unknown keys ['color']."]
```

## 根因
`figure_scope._validate()` 对同一个 brief 返回 `[]`（通过），`figure_draw._validate()` 返回上述错误。
runner 只调了后者。

## 建议修法
按 `brief["kind"]` 分派：

```python
kind = (brief or {}).get("kind") or "scheme"
if kind == "scope":
    from cdxml_toolkit.mcp_runtime.figure_scope import draw_scope
    report = draw_scope(brief, output_path=str(case_dir), render=bool(render))
else:
    report = draw_figure(brief, output_path=str(case_dir), render=bool(render))
```

## 验证
打完补丁后（并用**装好 runtime 的那个解释器**跑）：
```
26/26 passed, readback=matched x26, chk_err 0, warn 0
```

## 备注
这个坑容易被误诊成"runtime 版本落后"，因为报错文本里没有提示"用错入口"。
建议在 `_validate` 的错误信息里带上 `kind`，例如
`Unknown brief keys for kind 'scope': [...]`。
