# auxetic-transferability — 受控 QE 计算

论文《拉胀判定的不可迁移性》补算的 QE 输入与 GitHub Actions 工作流。
仓库只存输入，输出以 artifact 形式取回，不进版本库。

## 层级

| 层 | 内容 | job 数 | SCF |
|---|---|---|---|
| A | HfRh，展宽 0.010 Ry @ 20³ / 24³，PBEsol @ 24³ | 6 | 12 |
| B | InP(立方) / BN(六方) / SnS(正交) / HfAl(正交) / Zr₂Cu(四方)，k = 6³/12³/16³ | 78 | 156 |
| C | 同上 5 个材料，k = 24³ | 26 | 52 |

一个 job = 一个应变模式 = 两个幅度（±），对应 `m.in` / `p.in`。

## 计算约定

与 `03_程序/qe_k24_strain.ipynb` 完全一致，保证新旧数字可直接同表：

- `eps[j,j] = s*h`（j=0,1,2，单轴），`eps[k,l] = 0.5*s*h`（j=3,4,5，剪切），`cell' = cell @ (I + eps)`，`h = 0.01`
- `C[:, j] = (σ(+) − σ(−)) / (2h)`，`σ = −0.1 × QE 应力`（kbar → GPa，压为正取反）
- 立方只跑 j=0（C11/C12）与 j=3（C44），其余晶系跑 6 个 Voigt 模式
- 固定离子 `calculation='scf'`；`ecutwfc 60 / ecutrho 480`；`smearing='mv', degauss=0.020`（A 层 0.010 两档）；`conv_thr 1e-8`；`disk_io='none'`
- 赝势统一 PSLibrary 1.0.0 的 PBE 版，运行时从官方站点下载，不入库

## 运行

Actions → **QE SCF** → Run workflow：

- `plan` = `A` / `B` / `C` / `ALL`
- `only` = 可选子串过滤，如 `_k24_`
- `max_parallel` 默认 20（免费账户上限）

结果：run 页面底部 110 个 artifact，每个含 `m.out` / `p.out`。

## 本地脚本

```
python scripts/fetch_structures.py   # 从 MP 取常规胞（需 D:\aiwork\cail\.mp_key）
python scripts/gen_inputs.py         # 生成 runs/ 与 plan/jobs.json
python scripts/make_matrix.py --plan A
python scripts/parse_results.py      # 汇总 artifact → results/qe_summary.csv
```
