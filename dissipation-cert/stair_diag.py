#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""stair_diag.py —— 数学层误差沿 x 的可视化诊断：分段提取 + 抓异常 + 中心值复测。
不修改 dissipation_cert.py。跑法:
    python3 stair_diag.py

产出（当前目录）:
    stair_diag.csv   原始采样
    stair_diag.svg   双对数误差曲线（纯标准库生成，无第三方依赖）
    stdout           分段表 / 异常清单 / 中心值复测表

关键：级数段的界由本脚本自己递推 H_n 得到（精确有理，第二实现），
      与 dc 的 htab 互为对撞；直算段用两档 prec 对撞实测。
"""
from fractions import Fraction as Fr
from math import factorial, log10, sqrt
import dissipation_cert as dc

NK   = getattr(dc, "NK", 24)        # 必须与 dc 内部一致，不一致会算错台阶
SAFE = getattr(dc, "SAFE", Fr(1, 2))
PREC_LO, PREC_HI = 10, 20           # 对撞用的两档 prec
NOISE_LO, NOISE_HI = 1e-6, 1e-3     # 物理噪声预算
SEGS = [(1e-6, 1e-2), (1e-2, 0.1), (0.1, 0.25), (0.25, 0.5),
        (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0)]


# ---------- 第二实现：H_n 精确递推（与 dc.htab 对撞用） ----------
def h_coeffs(N):
    """F(x)=1/(sin x - x)=sum_{n>=0} H_n x^(2n-3)。
       sin x - x = sum_{k>=1} a_k x^(2k+1),  a_k = (-1)^k/(2k+1)!  """
    a = [None] + [Fr((-1) ** k, factorial(2 * k + 1)) for k in range(1, N + 2)]
    inv_a1, H = 1 / a[1], [Fr(0)] * (N + 1)
    H[0] = inv_a1
    for n in range(1, N + 1):
        s = Fr(0)
        for j in range(1, n + 1):
            s += a[j + 1] * H[n - j]
        H[n] = -inv_a1 * s
    return H


H = h_coeffs(NK + 2)
assert H[1] == Fr(-3, 10) and H[2] == Fr(-11, 1400) and H[3] == Fr(-17, 126000), \
    "H 递推与已知系数不符，先停：%s" % [H[1], H[2], H[3]]


# ---------- 两个分支的误差量级 ----------
def eps_series(x, K=NK):
    """级数分支截断上界；与 prec 无关。x>SAFE 时该分支不启用，返回 0。
    口径与 dc.eps_series 一致：R_series 保留 H_1..H_{K-1}，首丢项 H_K x^(2K-3)，
    余项 <= |T_K|(1+x^2+x^4+...) = |H_K| x^(2K-3)/(1-x^2) <= |H_K| x^(2K-3)·4/3。"""
    if x > SAFE:
        return Fr(0)
    return abs(H[K]) * x ** (2 * K - 3) * Fr(4, 3)


PI_LO, PI_HI = dc.pi_frac(PREC_LO), dc.pi_frac(PREC_HI)
def eps_direct(x):
    """直算分支误差量级 = |R(prec_hi) - R(prec_lo)|。级数分支下两者恒等 → 0。"""
    try:
        return abs(dc.R(x, PREC_HI, PI_HI) - dc.R(x, PREC_LO, PI_LO))
    except Exception:
        return None


def branch(x):
    return "series" if x <= SAFE else "direct"


def eps_used(x):
    """代码实际给出的界：走哪条分支就用哪条的。"""
    v = eps_direct(x)
    if branch(x) == "series":
        return eps_series(x)                     # 直算值在级数分支下恒为 0，不参与
    return v if v is not None else eps_series(x)


# ---------- 采样与分段 ----------
def xs_from(a, b, n):
    la, lb = log10(a), log10(b)
    return [Fr("%.10f" % (10.0 ** (la + (lb - la) * i / (n - 1)))) for i in range(n)]


CURVE_X = xs_from(1e-6, 8.0, 160)
rows = []
for x in CURVE_X:
    es, ed, eu = eps_series(x), eps_direct(x), eps_used(x)
    rows.append(dict(x=float(x), br=branch(x), eps_series=float(es),
                     eps_direct=(float(ed) if ed is not None else None),
                     eps_used=float(eu),
                     Rp=abs(float(dc.dR(x, PREC_HI, PI_HI)))))

# ---------- 分段提取 + 异常判定 ----------
print("=" * 96)
print("分段提取（每段对数均匀 9 点，中心取几何中点）")
print("%-16s %-10s %-12s %-12s %-12s %s"
      % ("段", "中心 x", "中心 eps", "段内 max/min", "中心|R'|", "标记"))
print("-" * 96)
anomalies = []
seg_report = []
for a, b in SEGS:
    xs = xs_from(a, b, 9)
    xc = Fr("%.10f" % sqrt(a * b))
    vals = [float(eps_used(x)) for x in xs]
    vmax, vmin = max(vals), max(min(vals), 1e-300)
    ratio = vmax / vmin
    ec = float(eps_used(xc))
    rp = abs(float(dc.dR(xc, PREC_HI, PI_HI)))
    flags = []
    if ratio > 1e4:
        flags.append("段内跨 %.0e 倍" % ratio)
    if a < 0.5 < b:
        flags.append("跨分支阈值")
    if ec > NOISE_LO:
        flags.append("仍在噪声之上")
    else:
        flags.append("已退到噪声下")
    seg_report.append(dict(a=a, b=b, xc=float(xc), ec=ec, ratio=ratio, Rp=rp))
    if ratio > 1e4:
        anomalies.append(("段内剧烈变化", a, b, ratio))
    if a < 0.5 < b:
        anomalies.append(("跨分支阈值", a, b, ratio))
    print("%-16s %-10.4f %-12.3e %-12.3e %-12.3e %s"
          % ("(%.0e,%.2f]" % (a, b), float(xc), ec, ratio, rp, "; ".join(flags)))

# ---------- 抓台阶：找 eps_used 的最大跳变 ----------
jump, jx = 0.0, None
for i in range(1, len(rows)):
    p, q = rows[i - 1], rows[i]
    if p["eps_used"] > 0 and q["eps_used"] > 0:
        r = abs(log10(p["eps_used"]) - log10(q["eps_used"]))
        if r > jump:
            jump, jx = r, q["x"]
print("-" * 96)
print("最大跳变: %.2f 个数量级  @ x ≈ %.6f   (分支阈值 SAFE = %.3f)"
      % (jump, jx if jx else -1, float(SAFE)))
print("级数分支首丢项 |H_%d| = %.6e" % (NK, float(abs(H[NK]))))
print("级数分支 @x=1/2 : %.3e   （含尾因子 4/3）" % float(eps_series(Fr(1, 2))))
print("级数分支 @x=1/2 : %.3e   （纯首丢项，口径对照）"
      % float(abs(H[NK]) * Fr(1, 2) ** (2 * NK - 3)))

# ---------- 中心值复测 ----------
print("=" * 96)
print("中心值复测：把每段中心当工作点，看数学项能否退到噪声预算之下")
print("%-14s %-10s %-12s %-14s %-12s %s"
      % ("段中心", "分支", "数学项", "等效物理误差", "噪声下界", "结论"))
print("-" * 96)
for s in seg_report:
    phys = s["ec"] * max(s["Rp"], 1e-300)
    ok = "✔ 已在噪声下" if phys < NOISE_LO else "✘ 仍在噪声上"
    print("%-14.4f %-10s %-12.3e %-14.3e %-12.1e %s"
          % (s["xc"], branch(Fr("%.10f" % s["xc"])), s["ec"], phys, NOISE_LO, ok))

# ---------- 风险段中心（(0.25,0.5] 的几何中点）单点细测 ----------
xr = Fr("%.10f" % sqrt(0.25 * 0.5))
print("-" * 96)
print("风险段中心 x = %.6f" % float(xr))
print("  级数截断       = %.3e" % float(eps_series(xr)))
print("  直算对撞       = %s" % ("%.3e" % float(eps_direct(xr)) if eps_direct(xr) is not None else "n/a"))
print("  |R'|           = %.3e" % abs(float(dc.dR(xr, PREC_HI, PI_HI))))
print("  数学项×|R'|    = %.3e" % (float(eps_series(xr)) * abs(float(dc.dR(xr, PREC_HI, PI_HI)))))
print("  噪声预算       = %.1e ~ %.1e" % (NOISE_LO, NOISE_HI))

# ---------- CSV + SVG ----------
with open("stair_diag.csv", "w",encoding="utf-8") as f:
    f.write("x,branch,eps_series,eps_direct,eps_used,Rp\n")
    for r in rows:
        f.write("%.10e,%s,%.10e,%s,%.10e,%.10e\n"
                % (r["x"], r["br"], r["eps_series"],
                   ("%.10e" % r["eps_direct"]) if r["eps_direct"] is not None else "",
                   r["eps_used"], r["Rp"]))


def _svg(curves, path, xr=(1e-6, 8.0), W=1060, H=650):
    ML, MR, MT, MB = 104, 158, 46, 74
    pw, ph = W - ML - MR, H - MT - MB
    lx0, lx1 = log10(xr[0]), log10(xr[1])
    ys = [y for _, _, _, p in curves for _, y in p if y and y > 0]
    ly0, ly1 = log10(min(ys)), log10(max(ys))
    pd = (ly1 - ly0) * 0.08
    ly0, ly1 = ly0 - pd, ly1 + pd
    X = lambda v: ML + (log10(v) - lx0) / (lx1 - lx0) * pw
    Y = lambda v: MT + (1 - (log10(max(v, 10 ** ly0)) - ly0) / (ly1 - ly0)) * ph
    s = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
         'viewBox="0 0 %d %d">' % (W, H, W, H),
         '<rect width="%d" height="%d" fill="#fff"/>' % (W, H)]
    # 噪声带
    s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#ffe9a8"/>'
             % (ML, Y(NOISE_HI), pw, Y(NOISE_LO) - Y(NOISE_HI)))
    s.append('<text x="%.1f" y="%.1f" font-size="12" fill="#8a6d1a">'
             '物理噪声预算 1e−6 ~ 1e−3</text>' % (ML + 8, Y(NOISE_HI) - 6))
    # 网格
    for e in range(int(ly0) - 1, int(ly1) + 2):
        yy = MT + (1 - (e - ly0) / (ly1 - ly0)) * ph    # Y(10^e)，避免 10.0**e 下溢
        if MT <= yy <= MT + ph:
            s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee"/>'
                     % (ML, yy, ML + pw, yy))
            s.append('<text x="%.1f" y="%.1f" font-size="11" text-anchor="end" fill="#666">'
                     '1e%d</text>' % (ML - 8, yy + 4, e))
    for e in range(int(lx0), int(lx1) + 1):
        xv = 10.0 ** e
        if xr[0] <= xv <= xr[1]:
            s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#eee"/>'
                     % (X(xv), MT, X(xv), MT + ph))
            s.append('<text x="%.1f" y="%.1f" font-size="11" text-anchor="middle" '
                     'fill="#666">%g</text>' % (X(xv), MT + ph + 18, xv))
    # 分支阈值
    s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#c00" '
             'stroke-dasharray="5,4"/>' % (X(0.5), MT, X(0.5), MT + ph))
    s.append('<text x="%.1f" y="%.1f" font-size="12" fill="#c00">x=1/2 分支切换</text>'
             % (X(0.5) + 6, MT + 14))
    # 曲线
    for label, color, wd, pts in curves:
        d = " ".join("%.1f,%.1f" % (X(x), Y(y)) for x, y in pts if y and y > 0)
        s.append('<polyline fill="none" stroke="%s" stroke-width="%s" '
                 'points="%s"/>' % (color, wd, d))
    # 图例
    for i, (label, color, wd, _) in enumerate(curves):
        yy = MT + 8 + i * 22
        s.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
                 'stroke-width="%s"/>' % (ML + pw + 14, yy, ML + pw + 44, yy, color, wd))
        s.append('<text x="%.1f" y="%.1f" font-size="12" fill="#333">%s</text>'
                 % (ML + pw + 50, yy + 4, label))
    s.append('<text x="%.1f" y="%.1f" font-size="13" text-anchor="middle" fill="#333">'
             '数学层误差沿 x（双对数，K=%d）</text>' % (ML + pw / 2, 26, NK))
    s.append('</svg>')
    open(path, "w",encoding="utf-8").write("\n".join(s))


_ser = [(r["x"], r["eps_series"]) for r in rows if r["x"] <= float(SAFE)]
_dir = [(r["x"], r["eps_direct"]) for r in rows
        if r["eps_direct"] and r["x"] > float(SAFE)]
_use = [(r["x"], r["eps_used"]) for r in rows]
_svg([("级数分支截断 |H_%d|·x^%d (与 prec 无关)" % (NK, 2 * NK - 3), "#1f6fb4", 2, _ser),
      ("直算分支两档 prec 对撞", "#d9480f", 2, _dir),
      ("代码实际给出的界", "#222", 3, _use)],
     "stair_diag.svg")
print("=" * 96)
print("已写出 stair_diag.csv / stair_diag.svg（浏览器直接打开看）")
for tag, a, b, v in anomalies:
    print("异常: %s  @ 段 (%.0e, %.2f]  比值 %.1e" % (tag, a, b, v))