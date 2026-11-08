import math
import sys
import time
from fractions import Fraction as Fr
from math import factorial


# ---------- 1. 高精度常数与三角函数（全程精确有理数） ----------

def atan_inv(n, prec):
    """arctan(1/n) = sum_{k>=0} (-1)^k / ((2k+1) n^(2k+1))
    交错级数，截断误差 < 首个被丢项 < 1e-prec。"""
    eps, s, k, sg = Fr(1, 10**prec), Fr(0), 0, 1
    while True:
        t = Fr(1, (2*k+1) * n**(2*k+1))
        if t < eps:
            break
        s += sg*t
        sg = -sg
        k += 1
    return s


def pi_frac(prec):
    """Machin: pi = 16*arctan(1/5) - 4*arctan(1/239)。
    两个分量各按 1e-(prec+10) 截断（误差 < 首个被丢项 < eps），
    故 |Δpi| < 16*eps + 4*eps = 20e-(prec+10) = 2e-(prec+9)。"""
    return 16*atan_inv(5, prec+10) - 4*atan_inv(239, prec+10)


def dpi_bound(prec):
    """pi_frac(prec) 的绝对误差上界：< 2e-(prec+9)。sin/cos 归约时要用到它。"""
    return Fr(2, 10**(prec+9))


def reduce_x(x, pi):
    """把 x 归约到 (-pi, pi]：返回 (xr, k)，使 x = xr + k*(2pi)（可能再减一次 2pi）。
    单独提出来给 sin_ex 和 (f) 共用，避免两处的归约逻辑各写一遍而漂移。"""
    tp = 2*pi
    k = int(x/tp)
    xr = x - k*tp
    if xr > pi:
        xr -= tp
        k += 1
    return xr, k


def _sin_series(xr, prec):
    """|xr| <= pi 时的 sin 交错 Taylor 级数：
    返回 (部分和, 首个被丢项的绝对值, 最后加进去的那一项的绝对值)。

    循环是「先加项、再判停」，退出时最后加进去的项 |T_N| <= eps = 1e-prec，
    但真余项 <= 首个被丢项 |T_(N+1)| = |T_N| * xr^2/((2N)(2N+1)) —— 比 eps 紧得多。
    |xr| <= pi 时自 m=1 起 |T_m| 单调下降（相邻比 <= pi^2/(20) < 1），
    交错级数余项不超过首个被丢项，所以这是真上界。

    第三个返回值给 (f) 用：m := 首丢项/末项，且由 m*eps < |T_N| <= eps 可得
    log10(eps/Δsin) ∈ [log10(1/m), 2*log10(1/m))。"""
    eps, t, s, n = Fr(1, 10**prec), xr, xr, 1
    while abs(t) > eps:
        t = -t*xr*xr / ((2*n)*(2*n+1))
        s += t
        n += 1
    last = abs(t)
    first_dropped = abs(-t*xr*xr / ((2*n)*(2*n+1)))
    return s, first_dropped, last


def sin_ex(x, prec, pi, dpi):
    """返回 (sin x, 误差上界)。
    上界 = 级数首个被丢项 + 2*pi 归约带来的自变量偏差（|d sin / d x| = |cos| <= 1）。
    归约 xr = x - k*(2pi) 里 pi 自带 dpi 的误差，使 xr 偏掉 (|k|+1)*2*dpi。"""
    xr, k = reduce_x(x, pi)
    s, tail, _last = _sin_series(xr, prec)
    return s, tail + (abs(k)+1)*2*dpi


def sin_frac_ex(x, prec, pi):
    return sin_ex(x, prec, pi, dpi_bound(prec))


def sin_frac(x, prec, pi):
    return sin_ex(x, prec, pi, dpi_bound(prec))[0]


def cos_frac_ex(x, prec, pi):
    """cos x = sin(x + pi/2)。级数按 prec+5 截断，但参数里的 pi 只按 prec 位，
    故归约项按 prec 记（dpi_bound(prec)），另加 pi/2 自身的偏差 dpi/2。
    本文件只在 x <= 10 用它（见 dR/d2R 的 assert）。"""
    s, tail = sin_ex(x + pi/2, prec+5, pi, dpi_bound(prec))
    return s, tail + dpi_bound(prec)/2


def cos_frac(x, prec, pi):
    return cos_frac_ex(x, prec, pi)[0]


def Si_frac_ex(x, prec):
    """Si(x) = sum_{k>=1} (-1)^(k-1) x^(2k-1) / ((2k-1)(2k-1)!)，返回 (部分和, 首个被丢项)。

    相邻项 T_(k+1)/T_k = -x^2 (2k-1)/((2k)(2k+1)^2)。
    精确有理运算：没有舍入，只有截断。循环在 |T| < eps 时立即返回，
    所以最后加进去的是 T_M（|T_M| < eps），首个被丢项再乘一个因子。
    代价随 x 走（项数 ~ e*x/2）；只在小 x 用它，再大需要渐近展开，这里直接拒绝。"""
    if x > 20:
        raise ValueError("Si_frac 只在小 x 使用（见函数注释）；收到 x=%s" % x)
    eps, t, s, k = Fr(1, 10**prec), x, x, 1
    while True:
        t = -t*x*x*(2*k-1) / (2*k*(2*k+1)**2)
        s += t
        if abs(t) < eps:
            first_dropped = abs(-t*x*x*(2*k+1) / ((2*k+2)*(2*k+3)**2))
            return s, first_dropped
        k += 1


def Si_frac(x, prec):
    return Si_frac_ex(x, prec)[0]


def dec(fr, digits):
    """按小数截断打印精确有理数（走整数除法，不经浮点，所以 digits 位都是真的）。
    注意是截断不是四舍五入：残差落在 [0, 10^-digits)，上界是 1 个单位末位而非半个。"""
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    scaled = fr.numerator * 10**digits // fr.denominator
    s = str(scaled).rjust(digits+1, "0")
    return "%s%s.%s" % (sign, s[:-digits], s[-digits:])


# ---------- 2. 精确 Laurent 系数（零误差） ----------

_TAB = {}


def htab(K):
    """H[0..K-1]: 1/(sin x - x) = sum_n H_n x^(2n-3)
    h_j = (-1)^(j+1)/(2j+3)! 是 (sin x - x)/x^3 的系数；
    形式幂级数求逆满足 sum_{i=0..n} h_i H_(n-i) = delta_(n0)。
    手算双路径交叉核对：
        H_0=-6, H_1=-3/10, H_2=-11/1400, H_3=-17/126000。"""
    if K not in _TAB:
        h = [Fr((-1)**(j+1), factorial(2*j+3)) for j in range(K)]
        H = [Fr(0)]*K
        H[0] = 1/h[0]
        for n in range(1, K):
            H[n] = -sum(h[i]*H[n-i] for i in range(1, n+1)) / h[0]
        _TAB[K] = H
    return _TAB[K]


NK = 24
SAFE = Fr(1, 2)     # R 的级数/直算分支阈值


def eps_series(x, K=NK):
    """x <= SAFE 走 R_series 时的级数截断上界；与 prec 无关，x > SAFE 恒为 0。

    R_series 保留 H_n 到 n = K-1，首个被丢项 T_K = H_K x^(2K-3)。
    自 n >= K 起 |T_(n+1)/T_n| <= x^2（x <= 1/2 时 <= 1/4，见 R_series 注释），
    余项 <= |T_K|·(1 + x^2 + x^4 + ...) = |H_K|·x^(2K-3)·1/(1-x^2)。
    x <= 1/2 时 1/(1-x^2) <= 4/3。"""
    if x > SAFE:
        return Fr(0)
    h = abs(htab(K + 1)[K])                # |H_K| = 首丢项系数
    return h * abs(x)**(2*K - 3) * Fr(4, 3)


def R_series(x, K=NK):
    """R(x) = sum_{n>=2} H_n x^(2n-3)
            = sum_{n=1..K-1} H_n x^(2n-3) + (3/10)/x

    这里 n 从 1 起写、末尾补 (3/10)/x，是因为 H_1*x^-1 = -3/(10x) 与之精确相消，
    结果等于「从 n=2 起」的级数。而 -6/x^3 已经被 R 的定义减掉了，绝不能再加回来。
    首项 -11x/1400；首个被丢项 H_K x^(2K-3)。"""
    Hk = htab(K)
    return sum(Hk[n]*x**(2*n-3) for n in range(1, K)) + Fr(3, 10)/x


def R(x, prec, pi, safe=Fr(1, 2)):
    """x <= safe 走精确级数（零舍入），x > safe 走直算。
    safe=1/2 的依据：直算路径在 x 处的相消比 |F|/|R| = 764/x^4，
        x=1e-4 -> 7.6e18（丢 19 位），x=1/2 -> 1.2e4（只丢 4 位）；
    级数路径在 x=1/2 处的截断量级 ≈3.2e-56（NK=24，见 eps_series/selfcheck [2]）。"""
    if x <= safe:
        return R_series(x)
    return 1/(sin_frac(x, prec, pi) - x) + 6/x**3 + Fr(3, 10)/x
def dR(x, prec, pi):
    """R'(x) = (1 - cos x)/(sin x - x)^2 - 18/x^4 - (3/10)/x^2"""
    assert x <= 10, "dR 里 cos_frac 只补偿 prec+5，x 再大归约会吃掉更多位"
    sx = sin_frac(x, prec, pi)
    return (1 - cos_frac(x, prec, pi))/(sx - x)**2 - 18/x**4 - Fr(3, 10)/x**2


def d2R(x, prec, pi):
    """R''(x) = sin x/(sin x - x)^2 + 2(1-cos x)^2/(sin x - x)^3 + 72/x^5 + 3/(5x^3)
    给 (e) 块定标用：末端残差 |R'(x*)| ~ |R''(x*)|·(区间宽)。"""
    assert x <= 10, "d2R 里 cos_frac 只补偿 prec+5，x 再大归约会吃掉更多位"
    sx, cx = sin_frac(x, prec, pi), cos_frac(x, prec, pi)
    return sx/(sx - x)**2 + 2*(1-cx)**2/(sx - x)**3 + 72/x**5 + Fr(3, 5)/x**3


def d2Phi(x, prec, pi):
    """Phi'' = R'' - (sin x / x)' = R'' - (cos x/x - sin x/x^2)"""
    return d2R(x, prec, pi) - cos_frac(x, prec, pi)/x + sin_frac(x, prec, pi)/x**2


def noise_dR(x, prec, pi, bound=False):
    """dR 在 x 处的求值噪声上界：
        |Δ(dR)| <= |∂dR/∂sin|·Δsin + |∂dR/∂cos|·Δcos
        |∂dR/∂sin| = 2(1-cos x)/|sin x - x|^3
        |∂dR/∂cos| = 1/(sin x - x)^2
    bound=False（默认）：Δsin、Δcos 用实际的首个被丢项（真上界，比 1e-prec 紧）；
    bound=True          ：退回 v6 的保守口径 1e-prec / 1e-(prec+5)，只作对照。
    返回 (总噪声, sens_sin, sens_cos, Δsin, Δcos)。"""
    if bound:
        sx, ds = sin_frac(x, prec, pi), Fr(1, 10**prec)
        cx, dc = cos_frac(x, prec, pi), Fr(1, 10**(prec+5))
    else:
        sx, ds = sin_frac_ex(x, prec, pi)
        cx, dc = cos_frac_ex(x, prec, pi)
    sens_sin = 2*(1-cx)/abs(sx-x)**3
    sens_cos = 1/(sx-x)**2
    return sens_sin*ds + sens_cos*dc, sens_sin, sens_cos, ds, dc


def noise_phi(x, prec, pi, bound=False):
    """Phi' = R' - sin x/x，故 ∂Phi'/∂sin = ∂R'/∂sin - 1/x，∂Phi'/∂cos = ∂R'/∂cos。
    两个 sin 灵敏度取绝对值之和（保守）。(d) 里 sin x/x 那部分用的是 prec+5 档的级数，
    这里统一用 prec 档的 Δsin，只会更大，所以仍是上界。"""
    if bound:
        sx, ds = sin_frac(x, prec, pi), Fr(1, 10**prec)
        cx, dc = cos_frac(x, prec, pi), Fr(1, 10**(prec+5))
    else:
        sx, ds = sin_frac_ex(x, prec, pi)
        cx, dc = cos_frac_ex(x, prec, pi)
    sens_sin = 2*(1-cx)/abs(sx-x)**3 + 1/x
    sens_cos = 1/(sx-x)**2
    return sens_sin*ds + sens_cos*dc, sens_sin, sens_cos, ds, dc


# ---------- 2b. 双侧窗口：余项括号 + station 的噪声通道 ----------

def _series_bracket(xr, prec, extra):
    """返回 (S_N, |T_(N+1)|, m, lo, hi)。

    S_N 是交错 Taylor 到判停处的前 N 项和；d = |T_(N+1)| 就是 v7/v8 用的首丢项界。
    再用 extra 个已知项把余项 R_N = sum_{k>N} T_k 夹成 [lo, hi]：
        R_N = (T_(N+1)+...+T_(N+extra)) + tail，
        tail 的首两项为 u、u2，交错且单调降 => tail ∈ [u+u2, u]，
        lo/hi = known + 该区间（按大小排序）。
    extra=0 时 lo/hi = [T_(N+1)+T_(N+2), T_(N+1)]，且 |R_N| <= |T_(N+1)| = d。
    m = d / |T_N| 是停止比，用来解释「为什么实际误差不是 1e-prec」。"""
    eps = Fr(1, 10**prec)
    t, s, n = xr, xr, 1
    while abs(t) > eps:
        t = -t*xr*xr / ((2*n)*(2*n+1))
        s += t
        n += 1
    last = abs(t)
    d = abs(-t*xr*xr / ((2*n)*(2*n+1)))
    m = (d/last) if last else Fr(0)

    cur = -t*xr*xr / ((2*n)*(2*n+1))        # 带符号的 T_(N+1)
    known = Fr(0)
    idx = n
    for _ in range(extra):
        known += cur
        idx += 1
        cur = -cur*xr*xr / ((2*idx)*(2*idx+1))
    u2 = -cur*xr*xr / ((2*(idx+1))*(2*(idx+1)+1))
    a, b = known + cur + u2, known + cur
    lo, hi = (a, b) if a < b else (b, a)
    return s, d, m, lo, hi


def _sin_full(x, prec, pi, extra):
    """(值, 首丢项, 收紧比, 停止比, 余项下界, 余项上界, 归约误差)。
    收紧比 r := max(|lo|,|hi|)/d <= 1，故 |真余项| <= r*d 仍是证明。"""
    xr, k = reduce_x(x, pi)
    s, d, m, lo, hi = _series_bracket(xr, prec, extra)
    red = (abs(k)+1)*2*dpi_bound(prec)
    r = (max(abs(lo), abs(hi))/d) if d else Fr(0)
    return s, d, r, m, lo, hi, red


def _cos_full(x, prec, pi, extra):
    """cos x = sin(x+pi/2)：级数按 prec+5，归约按 prec，另加 pi/2 的偏差。"""
    xr, k = reduce_x(x + pi/2, pi)
    s, d, m, lo, hi = _series_bracket(xr, prec+5, extra)
    red = (abs(k)+1)*2*dpi_bound(prec) + dpi_bound(prec)/2
    r = (max(abs(lo), abs(hi))/d) if d else Fr(0)
    return s, d, r, m, lo, hi, red


def sin_tight(x, prec, pi, extra=4):
    """按你片段里的四元组：(收紧比, 停止比, 余项下界, 余项上界)。"""
    _s, _d, r, m, lo, hi, _red = _sin_full(x, prec, pi, extra)
    return r, m, lo, hi


def cos_tight(x, prec, pi, extra=4):
    _s, _d, r, m, lo, hi, _red = _cos_full(x, prec, pi, extra)
    return r, m, lo, hi


def station_terms(x, prec, pi, extra=4):
    """dR 在 x 处的噪声通道（对应你贴的那段 dict），外加认证带所需的键。

    n_cons  1e-prec 口径（v6，保守证明上界）
    n_est   首丢项口径（v7/v8）
    n_tight 余项括号口径 = ss*dsin*rsin + sc*dcos*rcos（仍为证明，且更紧）
    n_red   pi 归约误差传播（上面的 n_* 只覆盖级数截断，不含它）
    认证总数看 n_tight + n_red。"""
    s_sin, dsin, rsin, msin, lo_s, hi_s, red_s = _sin_full(x, prec, pi, extra)
    s_cos, dcos, rcos, mcos, lo_c, hi_c, red_c = _cos_full(x, prec, pi, extra)
    ss = 2*(1 - s_cos)/abs(s_sin - x)**3        # |∂dR/∂sin|
    sc = 1/(s_sin - x)**2                       # |∂dR/∂cos|
    dsin_c, dcos_c = Fr(1, 10**prec), Fr(1, 10**(prec+5))
    return {"ss": ss, "sc": sc, "dsin": dsin, "dcos": dcos,
            "dsin_c": dsin_c, "dcos_c": dcos_c,
            "rsin": rsin, "rcos": rcos, "msin": msin, "mcos": mcos,
            "nsinlo": lo_s, "nsinhi": hi_s, "ncoslo": lo_c, "ncoshi": hi_c,
            "n_est": ss*dsin + sc*dcos,
            "n_cons": ss*dsin_c + sc*dcos_c,
            "n_tight": ss*dsin*rsin + sc*dcos*rcos,
            "n_red": ss*red_s + sc*red_c,
            "sinv": s_sin, "cosv": s_cos,
            "sinlo": s_sin + lo_s - red_s, "sinhi": s_sin + hi_s + red_s,
            "coslo": s_cos + lo_c - red_c, "coshi": s_cos + hi_c + red_c,
            "red_s": red_s, "red_c": red_c}


def _R_from_sin(x, sv):
    return 1/(sv - x) + 6/x**3 + Fr(3, 10)/x


def _dR_from_sc(x, sv, cv):
    """把 dR 写成 (sin x, cos x) 的函数：dR = (1-cos)/(sin-x)^2 - 18/x^4 - 3/(10x^2)。
    于是给定 sin/cos 的认证区间，R' 的认证带就是矩形四角上的最小/最大值。"""
    return (1 - cv)/(sv - x)**2 - 18/x**4 - Fr(3, 10)/x**2


# ---------- 3. 变号二分（没变号就明确报失败，不静默） ----------

def bisect(a, b, f, it=150):
    """变号二分，返回中点。
    it 决定末端区间宽 (b-a)/2^it，因而决定末端残差 |f(x*)| ≈ |f''|·宽 ——
    那衡量的是「二分停在哪」，不是精度。真正的定位精度由 prec 决定（见 scan 的 (e) 块）。
    默认 150 只为兼容旧调用；scan 里按 it = 4*prec 传。"""
    fa, fb = f(a), f(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if (fa > 0) == (fb > 0):
        return None
    for _ in range(it):
        m = (a + b)/2
        fm = f(m)
        if fm == 0:
            return m
        if (fa > 0) == (fm > 0):
            a, fa = m, fm
        else:
            b, fb = m, fm
    return (a + b)/2


# ---------- 4. 模式 ----------

def do_exact(K=NK):
    print("(sin x - x)/x^3 = sum_j h_j x^(2j),   h_j = (-1)^(j+1)/(2j+3)!")
    for j in range(K):
        print("    h_%-2d  x^%-2d : %s" % (j, 2*j, Fr((-1)**(j+1), factorial(2*j+3))))
    print()
    print("F(x) = 1/(sin x - x) = sum_n H_n x^(2n-3),   sum_{i=0..n} h_i H_(n-i) = delta_(n0)")
    for n, c in enumerate(htab(K)):
        print("    H_%-2d  x^%+3d : %s" % (n, 2*n-3, c))
    print()
    print("R(x) = F + 6/x^3 + 3/(10x) = sum_{n>=2} H_n x^(2n-3)")
    print("    上面 H_0 x^-3 与 H_1 x^-1 两行就是被剥掉的发散项（H_1 = -3/10 与 +3/(10x) 精确相消）")
    print("    首项 = H_2 x = -11x/1400；首个被丢项 = H_%d x^%d" % (K, 2*K-3))


def do_probe(prec=30):
    pi = pi_frac(prec)
    # 采样表里 7.3242 贴近 R 的极小点，3.182 贴近 Phi 的驻点
    xs = [Fr(1,1000), Fr(1,100), Fr(1,10), Fr(1,2), Fr(1), Fr(2),
          Fr(3182,1000), Fr(73242,10000), Fr(20)]
    print("prec=%d   （x<=1/2 走级数，x>1/2 走直算）" % prec)
    for x in xs:
        v = R(x, prec, pi)
        print("    x=%-8s  R = %-22.15g  Phi = R - Si = %.15g"
              % (float(x), float(v), float(v - Si_frac(x, prec))))


def do_selfcheck(prec=30):
    pi = pi_frac(prec)
    K = NK
    tail = abs(htab(K+1)[K])          # H_K = 级数首个被丢项的系数

    print("[1] 两条独立路径对撞：级数（精确有理）vs 直算（高精度 sin）")
    for sf in [Fr(1,2), Fr(1), Fr(3,2), Fr(2)]:
        direct = 1/(sin_frac(sf, prec, pi) - sf) + 6/sf**3 + Fr(3,10)/sf
        ser = R_series(sf, K)
        print("    x=%-5s  R_series(K=%d) - R_direct = %-13.4g   首个被丢项 |H_%d| x^%d = %.3g"
              % (float(sf), K, float(ser - direct), K, 2*K-3, float(tail*sf**(2*K-3))))

    print("[2] 级数分支的截断量级 |H_%d| x^%d —— safe 能取多大由它决定" % (K, 2*K-3))
    for sf in [Fr(1,10), Fr(1,4), Fr(1,2), Fr(1), Fr(3,2), Fr(2)]:
        print("    x=%-5s  |H_%d| x^%d = %.3g" % (float(sf), K, 2*K-3, float(tail*sf**(2*K-3))))
    print("    代码取 safe = 1/2：该处截断已低于 prec=%d 的位面，再往上交给直算" % prec)

    print("[3] Si(pi) 对账单（参考值抄到 30 位；算出来的是精确有理数截断打印，没经浮点）")
    ref = Fr("1.851937051982466170361053370157")
    v, sitail = Si_frac_ex(pi, prec)
    print("    算得  %s" % dec(v, 30))
    print("    参考  1.851937051982466170361053370157")
    print("    |差| = %.3g" % float(abs(v - ref)))
    print("      ^ 这里含参考值自身截断到 30 位小数的误差，其上界是 <1e-30：")
    print("        截断（不是四舍五入）的残差 ∈ [0, 1e-30)，所以 |差| 落在 1e-30 以内都正常，")
    print("        不该拿 5e-31 去卡它（5e-31 是四舍五入的界，本文件是截断）。")
    print("    本级数首个被丢项 = %.3g  —— 真误差 <= 它，而不是 1e-%d" % (float(sitail), prec))
    print("    想更紧的话把参考值补到 35 位再对一次；本文件不硬写 30 位以外的数字，免得抄错。")
def do_bound():
    """R(x) < 0 中 x >= 4 的那半边：精确代数，不是抽样。
        |sin x| <= 1 且 x >= 4  =>  t := sin x / x in [-1/4, 1/4]
        F = 1/(sin x - x) = -1/x * 1/(1 - t),   而 1/(1-t) >= 1/(1 + 1/4) = 4/5
        =>  F <= -4/(5x)
        =>  R = F + 3/(10x) + 6/x^3 <= -1/(2x) + 6/x^3
        而 -1/(2x) + 6/x^3 < 0  <=>  x^2 > 12  <=>  x > 2*sqrt(3) = 3.4641..."""
    print("R(x) < 0 对 x >= 4 严格成立：")
    print("    |sin x| <= 1, x >= 4  =>  t = sin x / x in [-1/4, 1/4]")
    print("    F = -1/x * 1/(1-t),   1/(1-t) >= 1/(1+1/4) = 4/5")
    print("    R <= -4/(5x) + 3/(10x) + 6/x^3 = -1/(2x) + 6/x^3")
    for x in [Fr(4), Fr(5), Fr(10), Fr(100), Fr(10**6)]:
        lhs = Fr(-1, 2)/x + 6/x**3
        print("      x=%-10s  -1/(2x) + 6/x^3 = %-16.10g   < 0 ? %s"
              % (float(x), float(lhs), lhs < 0))
    print("    -1/(2x) + 6/x^3 < 0  <=>  x^2 > 12  <=>  x > 2*sqrt(3) = 3.46410...")
    print("    结论：R 在 [4, +inf) 上严格为负；x < 4 的部分仍属抽样。")


def do_station(x, prec=30, extra=4):
    """把 station_terms 的通道全部落地，并给出 R'(x) 的认证带。

    注意：ss、sc 是 ∂R'/∂sin、∂R'/∂cos，所以 n_* 界的是 R'(x) 的求值噪声，
    不是 R(x) 的。R'(x) 写成 (sin x, cos x) 的函数后对每个变量单调（或至少
    极值在矩形角上），故 sin/cos 的认证区间 -> R' 的认证区间就是四角最值。
    """
    pi = pi_frac(prec)
    T = station_terms(x, prec, pi, extra)
    sinlo, sinhi = T["sinlo"], T["sinhi"]
    coslo, coshi = T["coslo"], T["coshi"]
    v = T["sinv"]
    c = T["cosv"]

    print("station: x = %s   prec = %d   extra = %d" % (x, prec, extra))
    print("[1] 灵敏度（dR 对 sin/cos 的偏导绝对值）")
    print("    ss = 2(1-cos)/|sin-x|^3 = %.6g" % float(T["ss"]))
    print("    sc = 1/(sin-x)^2       = %.6g" % float(T["sc"]))
    print("[2] sin/cos 的余项窗口（首丢项 -> extra 项收紧）")
    print("    sin: 首丢项 %.4g  收紧比 r=%.6f  停止比 m=%.4g  认证 |余项| <= %.4g  括号宽 %.3g"
          % (float(T["dsin"]), float(T["rsin"]), float(T["msin"]),
             float(max(abs(T["nsinlo"]), abs(T["nsinhi"]))),
             float(T["nsinhi"] - T["nsinlo"])))
    print("    cos: 首丢项 %.4g  收紧比 r=%.6f  停止比 m=%.4g  认证 |余项| <= %.4g  括号宽 %.3g"
          % (float(T["dcos"]), float(T["rcos"]), float(T["mcos"]),
             float(max(abs(T["ncoslo"]), abs(T["ncoshi"]))),
             float(T["ncoshi"] - T["ncoslo"])))
    print("    pi 归约传播 red_s = %.4g,  red_c = %.4g" % (float(T["red_s"]), float(T["red_c"])))
    print("[3] 三个口径的噪声上界（都是 |ΔR'|）")
    print("    n_cons  = ss*1e-prec + sc*1e-(prec+5) = %.4g" % float(T["n_cons"]))
    print("    n_est   = ss*Δsin     + sc*Δcos       = %.4g" % float(T["n_est"]))
    print("    n_tight = ss*Δsin*r_s + sc*Δcos*r_c   = %.4g" % float(T["n_tight"]))
    print("    n_red   = ss*red_s    + sc*red_c      = %.4g  （一阶，方向性见下）"
          % float(T["n_red"]))
    print("    单调性 n_tight <= n_est <= n_cons：%s"
          % (T["n_tight"] <= T["n_est"] <= T["n_cons"]))
    print("    级数截断认证总数 n_tight = %.4g；把 pi 归约也算进带里（下面的带已含它）"
          % float(T["n_tight"]))
    print("[4] 认证带（把 sin/cos 区间直接推进 R'，不是线性化）")
    print("    sin x = %s...   窗口 ±%.3g" % (dec(v, prec), float(max(abs(T["nsinlo"]), abs(T["nsinhi"])) + T["red_s"])))
    print("    cos x = %s...   窗口 ±%.3g" % (dec(c, prec), float(max(abs(T["ncoslo"]), abs(T["ncoshi"])) + T["red_c"])))
    corners = [_dR_from_sc(x, sv, cv)
               for sv in (sinlo, sinhi) for cv in (coslo, coshi)]
    lo, hi = min(corners), max(corners)
    mid = _dR_from_sc(x, v, c)
    r_lo = hi - mid
    r_hi = mid - lo
    half = r_lo if r_lo > r_hi else r_hi
    print("    R'(x) = %.15g  ±  %.4g   （带 [%.15g - %.3g, %.15g + %.3g]）"
          % (float(mid), float(half), float(mid), float(half), float(mid), float(half)))
    rv = _R_from_sin(x, v)
    rhalf = max(_R_from_sin(x, sinhi) - rv, rv - _R_from_sin(x, sinlo))
    print("    R(x)  = %.15g  ±  %.4g   （R 关于 sin x 单调降，x>0 时 sin x < x）"
          % (float(rv), float(rhalf)))

    p2 = pi_frac(prec + 15)
    ref = _dR_from_sc(x, sin_frac(x, prec + 15, p2), cos_frac(x, prec + 15, p2))
    print("[5] 独立复核（prec+15 档，只当参考，不进证书）")
    print("    R'(x) 参考 = %.15g   |参考-带中心| = %.4g" % (float(ref), float(abs(ref - _dR_from_sc(x, v, c)))))
    print("    落在认证带内？ %s" % (lo <= ref <= hi))
    s2 = v/(v - x)**2 + 2*(1 - c)**2/(v - x)**3 + 72/x**5 + Fr(3, 5)/x**3
    print("    |R''(x)| = %.6g  ->  驻点定位 |Δx*| <= 半宽/|R''| ~ %.3g"
          % (float(abs(s2)), float(half/abs(s2))))
    print("    解读：band 是证明（sin/cos 区间已含 pi 归约），n_tight 只是截断那一段；")
    print("          两者之差就是去掉归约项后多出来的余量。")


def scan(N, prec):
    """抽样（不是全域证明）。
    (a) 1e-6 ~ 1e6 对数抽样：R 的整体轮廓与量级
    (b) [7.0, 7.8] 加密：只负责覆盖极小点所在区间
    (c) 二分 R' 的零点            -> R 的最小点（位置由二分给出，不由网格给出）
    (d) 二分 Phi' = R' - sin x/x  -> Phi 的驻点
    (e) 把「二分残差」和「求值噪声」分开：前者由 it 定，后者由 prec 定
    (f) 同一驻点、逐档 prec 复跑：噪声与 Δsin 为什么不按 1e-5 缩放

    长跑示例：scan 2000000 40。开销大头是 (b)：N/2 次精确有理 sin。
    精度由 (c)(d) 的二分给出，不靠 (b) 的密度。"""
    M = N//2
    if M < 1:
        sys.exit("N 太小（至少 2）")
    pi = pi_frac(prec)
    it = 4*prec                       # 二分次数与工作精度挂钩（见 (e)）
    t0 = time.time()
    vmin = vmax = xmin = xmax = None

    for i in range(M+1):                                  # (a)
        x = Fr(10**(-6 + 12.0*i/M))
        v = R(x, prec, pi)
        if vmin is None or v < vmin: vmin, xmin = v, x
        if vmax is None or v > vmax: vmax, xmax = v, x
    print("(a) %d 点对数抽样 [1e-6, 1e6]  minR=%.12g @ %.9g   maxR=%.12g @ %.9g  [%.0fs]"
          % (M+1, float(vmin), float(xmin), float(vmax), float(xmax), time.time()-t0))

    for i in range(M+1):                                  # (b)
        x = Fr(7.0 + 0.8*i/M)
        v = R(x, prec, pi)
        if v < vmin: vmin, xmin = v, x
        if v > vmax: vmax, xmax = v, x
    print("(b) 加密 [7.0,7.8] 后           minR=%.15g @ %.12g   maxR=%.15g @ %.12g  [%.0fs]"
          % (float(vmin), float(xmin), float(vmax), float(xmax), time.time()-t0))

    s1 = bisect(Fr(7), Fr(8), lambda t: dR(t, prec, pi), it=it)
    if s1 is None:
        print("(c) R 驻点：区间内无变号 —— 检查区间与 prec，别当成没有驻点")
    else:
        print("(c) R   驻点 x* = %.15g   R(x*)  = %.15g   |R'(x*)| = %.3g  [随 it 变，见 (e)]"
              % (float(s1), float(R(s1, prec, pi)), float(abs(dR(s1, prec, pi)))))
        print("    变号依据（区间两端对照）：R'(7) = %+.3g, R'(8) = %+.3g"
              % (float(dR(Fr(7), prec, pi)), float(dR(Fr(8), prec, pi))))

    s2 = bisect(Fr(3), Fr(4),
                lambda t: dR(t, prec, pi) - sin_frac(t, prec+5, pi)/t, it=it)
    if s2 is None:
        print("(d) Phi 驻点：区间内无变号 —— 检查区间与 prec")
    else:
        print("(d) Phi 驻点 x* = %.15g   Phi(x*) = %.15g   |Phi'(x*)| = %.3g  [随 it 变，见 (e)]"
              % (float(s2), float(R(s2, prec, pi) - Si_frac(s2, prec)),
                 float(abs(dR(s2, prec, pi) - sin_frac(s2, prec+5, pi)/s2))))
        print("    变号依据（区间两端对照）：Phi'(3) = %+.3g, Phi'(4) = %+.3g"
              % (float(dR(Fr(3), prec, pi) - sin_frac(Fr(3), prec+5, pi)/3),
                 float(dR(Fr(4), prec, pi) - sin_frac(Fr(4), prec+5, pi)/4)))

    print("(e) 零点定位精度：把「二分残差」和「求值噪声」分开")
    print("    二分跑了 it = 4*prec = %d 次，末端区间宽 = 2^-%d = %.3g；"
          % (it, it, float(Fr(1, 2**it))))
    print("    (c)(d) 打印的 |f'(x*)| ~ |f''(x*)|·(区间宽)，只说明「二分停在哪」：")
    print("    it 变大它就变小，与 prec 无关 —— 那不是精度证书。")
    print("    真正的定位精度 = 求值噪声 / |f''|。噪声按通道拆开：")
    for tag, s, d2f, nf in (("R", s1, d2R, noise_dR), ("Phi", s2, d2Phi, noise_phi)):
        if s is None:
            continue
        nz, ss, sc, ds, dc = nf(s, prec, pi)
        nzb = nf(s, prec, pi, bound=True)[0]
        f2 = abs(d2f(s, prec, pi))
        print("      %s 驻点：|%s''| = %.4g" % (tag, tag, float(f2)))
        print("        sin 通道  %.4g × %.3g(首个被丢项) = %.3g" % (float(ss), float(ds), float(ss*ds)))
        print("        cos 通道  %.4g × %.3g(首个被丢项) = %.3g" % (float(sc), float(dc), float(sc*dc)))
        print("        合计 %.3g  ->  |Δx| ~ %.2e   （退回保守界 1e-prec：%.3g -> |Δx| ~ %.2e）"
              % (float(nz), float(nz/f2), float(nzb), float(nzb/f2)))

    if s1 is not None:
        print("(f) 同一驻点、逐档 prec 复跑（pi 也按同一 prec 重算；x0 = (c) 定位到的 x*）")
        print("    prec   Δsin(首丢项+归约)  噪声(实际)        |Δx|~噪声/|R''|   噪声(保守1e-prec)   Δsin 低于1e-prec  比上一行")
        base = abs(d2R(s1, prec, pi))          # 固定归一化因子，隔离 prec 的影响
        prev = None
        lo = hi = None
        win_lo = win_hi = None
        for p in (10, 15, 20, 25, 30, 40):
            pp = pi_frac(p)
            ds_p = sin_frac_ex(s1, p, pp)[1]        # 含 pi 归约项（约 8e-(p+9)，远小于级数尾项）
            xr_p, _k = reduce_x(s1, pp)
            _s, dropped, last = _sin_series(xr_p, p)   # 纯级数口径：m = 首丢项/末项
            m = dropped/last
            nz_p = noise_dR(s1, p, pp)[0]
            nzb_p = noise_dR(s1, p, pp, bound=True)[0]
            orders = math.log10(float(Fr(1, 10**p)/ds_p))
            lo = orders if lo is None else min(lo, orders)
            hi = orders if hi is None else max(hi, orders)
            lgm = math.log10(float(1/m))
            win_lo = lgm if win_lo is None else min(win_lo, lgm)
            win_hi = 2*lgm if win_hi is None else max(win_hi, 2*lgm)
            ratio = "" if prev is None else "%.2g" % float(nz_p/prev)
            print("    %-5d  %-18.4g  %-15.4g  %-15.4g  %-15.4g  %6.2f 个数量级  %s"
                  % (p, float(ds_p), float(nz_p), float(nz_p/base), float(nzb_p), orders, ratio))
            prev = nz_p
        print("    实测跨度：Δsin 比 1e-prec 低 %.2f ~ %.2f 个数量级（不是固定倍数）。" % (lo, hi))
        print("    对照：保守列（1e-prec 的定义）每加 5 位严格缩 1e-5 —— 那是定义使然，不代表真实误差。")
        print("    机制：Δsin = T_(N+1) = m·T_N，m = xr^2/((2N)(2N+1))（xr 为归约后的自变量），")
        print("          而 m·eps < |T_N| <= eps（末项若连 m·eps 都不到，前一轮就该退出了），故")
        print("          log10(eps/Δsin) ∈ [log10(1/m), 2·log10(1/m))，本机该窗口 %.2f ~ %.2f 个数量级。"
              % (win_lo, win_hi))
        print("          窗口下端 =「末项刚过 eps 就停」，上端 =「末项贴着 eps 停」，")
        print("          逐档落点由 N 的相位决定 —— 所以相邻档的实测比可大可小，上表最后一列就是它。")

    print("抽样结论（不是全域证明）：minR=%.12g   maxR=%.12g" % (float(vmin), float(vmax)))
    print("参考（本机 prec=%d 实测；另经 mpmath 60 dps 独立复核）：" % prec)
    print("    R   驻点 x* ≈ 7.32415600462945,   R(x*)  ≈ -0.0985368172486791")
    print("        (b) 的加密网格只负责覆盖；驻点位置由 (c) 的二分给出。")
    print("        x = 7.35 不是驻点（R'(7.35) ≈ +6.1e-4），它只是网格上的一个采样点。")
    print("    Phi 驻点 x* ≈ 3.18194875344635,   Phi(x*) ≈ -1.8814967004858")
    print("    maxR ≈ -7.86e-9 @ x ≈ 1e-6    <- 对数抽样的左端点。R(0+)=0，")
    print("                                     所以 maxR 完全由窗口左端决定；")
    print("                                     x=1e6 那端 ≈ -0.7/x ≈ -7e-7，更负")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "exact"
    if mode == "exact":
        do_exact()
    elif mode == "probe":
        do_probe(int(sys.argv[2]) if len(sys.argv) > 2 else 30)
    elif mode == "selfcheck":
        do_selfcheck(int(sys.argv[2]) if len(sys.argv) > 2 else 30)
    elif mode == "bound":
        do_bound()
    elif mode == "station":
        if len(sys.argv) < 3:
            sys.exit("用法: station X prec [extra]   例如 station 7.3242 30 4")
        do_station(Fr(sys.argv[2]),
                   int(sys.argv[3]) if len(sys.argv) > 3 else 30,
                   int(sys.argv[4]) if len(sys.argv) > 4 else 4)
    elif mode == "scan":
        if len(sys.argv) < 4:
            sys.exit("用法: scan N prec   例如 scan 200000 30")
        scan(int(sys.argv[2]), int(sys.argv[3]))
    else:
        sys.exit("未知模式 %r；可用: exact | probe [prec] | selfcheck [prec] | bound | station X prec [extra] | scan N prec" % mode)