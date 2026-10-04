# -*- coding: utf-8 -*-
"""
CMP Particle Pressure Simulator  (pad asperity contact + single-abrasive contact)
=================================================================================
Pad 경도(Shore A / Shore D / Asker C), Pad Modulus, Head 압력, 슬러리 Solid contents를
입력받아 다음을 계산한다.

  [1] Pad 돌기(asperity) 접촉이 탄성인지 소성인지
        - 소성지수 ψ (Greenwood–Williamson)
        - 압력에 따른 소성 전이 (Chang–Etsion–Bogy 탄소성 모델) → 임계압력 P_crit
  [2] 입자 하나가 받는 하중 F_p 와 입자간 간격
        - Head 압력 분담 하중  vs  Pad 가 입자 하나를 밀 수 있는 최대 하중(경도 한계)
        - F_p 가 경도 한계에 걸리면 "Head 압력이 아닌 재료 경도에 의존" 하는 영역
  [3] 입자–Wafer 접촉압력 (Hertz) 과 막질 항복 여부
  [4] 스크래치/균열 임계하중 (Lawn–Evans 날카로운 압입, Auerbach/Hertz 둥근 압입)

모든 계산은 파이썬 표준 라이브러리만 사용한다 (tkinter 포함). 외부 패키지 불필요.
    python cmp_particle_pressure_sim_R01.py

모델/상수의 출처와 가정은 README_particle_sim.md 참조.
"""

import json
import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

APP_TITLE = "CMP Particle Pressure Simulator  R01"

# ------------------------------------------------------------------ 단위
PSI = 6894.757                        # Pa
UNIT_TO_PA = {"psi": PSI, "hPa": 100.0, "kPa": 1000.0}

# ------------------------------------------------------------------ 재료 기본값
# 입자: (밀도 g/cm3, E GPa, H GPa, ν, K_IC MPa·m^0.5)  — 대표값. 측정값으로 교체 권장.
PARTICLES = {
    "실리카 (SiO2)":    (2.20,  70.0,  8.0, 0.17, 0.75),
    "세리아 (CeO2)":    (7.20, 180.0,  7.0, 0.30, 1.00),
    "알루미나 (Al2O3)": (3.95, 380.0, 20.0, 0.23, 3.50),
    "직접 입력":        None,
}
# 막질: (E GPa, H GPa, ν, K_IC MPa·m^0.5 또는 None=연성)
FILMS = {
    "SiO2 (산화막)": ( 70.0,  8.0, 0.17, 0.75),
    "Cu":            (120.0,  1.5, 0.34, None),
    "W":             (410.0,  8.0, 0.28, 10.0),
    "poly-Si":       (160.0, 12.0, 0.22, 0.90),
    "Si3N4":         (260.0, 19.0, 0.25, 1.50),
    "직접 입력":     None,
}

# Asker C ↔ Shore A 근사 대응표 (듀로미터 비교표 기반, 근사치)
ASKER_C_TABLE = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
SHORE_A_TABLE = [ 0,  8, 18, 28, 40, 50, 62, 75, 88, 100]


# ================================================================== 경도 환산
def asker_c_to_shore_a(c):
    c = min(max(c, ASKER_C_TABLE[0]), ASKER_C_TABLE[-1])
    for i in range(len(ASKER_C_TABLE) - 1):
        c0, c1 = ASKER_C_TABLE[i], ASKER_C_TABLE[i + 1]
        if c0 <= c <= c1:
            a0, a1 = SHORE_A_TABLE[i], SHORE_A_TABLE[i + 1]
            return a0 + (a1 - a0) * (c - c0) / (c1 - c0)
    return SHORE_A_TABLE[-1]


def shore_a_to_E_gent(s):
    """Gent (1958): Shore A → Young's modulus [MPa]. 유효 범위 약 20~80."""
    s = min(max(s, 1.0), 99.0)
    return 0.0981 * (56 + 7.62336 * s) / (0.137505 * (254 - 2.54 * s))


def shore_to_E_ruess(s, scale):
    """Ruess 상관식: log10 E[MPa] = 0.0235·S − 0.6403 (Shore A), S→S+50 (Shore D)."""
    if scale == "Shore D":
        s = s + 50.0
    return 10 ** (0.0235 * s - 0.6403)


def hardness_to_modulus(value, scale):
    """경도 입력 → 환산 탄성계수 [MPa] 와 환산 경로 설명."""
    if scale == "Shore A":
        return shore_a_to_E_gent(value), f"Shore A {value:.0f} → Gent 식"
    if scale == "Shore D":
        return shore_to_E_ruess(value, "Shore D"), f"Shore D {value:.0f} → Ruess 식"
    if scale == "Asker C":
        sa = asker_c_to_shore_a(value)
        return shore_a_to_E_gent(sa), f"Asker C {value:.0f} ≈ Shore A {sa:.0f} → Gent 식"
    raise ValueError(scale)


# ================================================================== 수치 도구
def phi(s):
    return math.exp(-0.5 * s * s) / math.sqrt(2.0 * math.pi)


def simpson(f, a, b, n=400):
    if b <= a:
        return 0.0
    if n % 2:
        n += 1
    h = (b - a) / n
    tot = f(a) + f(b)
    for i in range(1, n):
        tot += (4 if i % 2 else 2) * f(a + i * h)
    return tot * h / 3.0


def plane_strain_E(E, nu):
    return E / (1.0 - nu * nu)


def contact_modulus(E1, nu1, E2, nu2):
    return 1.0 / ((1 - nu1 ** 2) / E1 + (1 - nu2 ** 2) / E2)


# ================================================================== Pad 돌기 접촉 (GW + CEB)
class PadContact:
    """Greenwood–Williamson 통계 접촉 + Chang–Etsion–Bogy 탄소성 돌기 모델.

    sigma: 돌기 높이 표준편차 [m], beta: 돌기 끝단 반경 [m], eta: 돌기 밀도 [1/m²]
    E_star: Pad–Wafer 접촉 계수 [Pa],  H: Pad 경도 [Pa],  nu: Pad 포아송비
    """

    def __init__(self, E_star, H, sigma, beta, eta, nu):
        self.E_star, self.H = E_star, H
        self.sigma, self.beta, self.eta = sigma, beta, eta
        K = 0.454 + 0.41 * nu
        self.omega_c = (math.pi * K * H / (2.0 * E_star)) ** 2 * beta   # 임계 간섭 [m]
        self.wc = self.omega_c / sigma
        self.psi = (E_star / H) * math.sqrt(sigma / beta)
        self._build_table()

    def _state(self, h):
        """무차원 분리 h 에서 단위 공칭면적당 하중·면적·접촉 돌기 수."""
        wc = self.wc
        up = max(h + 8.0, 8.0)
        el_P = simpson(lambda s: (s - h) ** 1.5 * phi(s), h, h + wc)
        el_A = simpson(lambda s: (s - h) * phi(s), h, h + wc)
        el_n = simpson(phi, h, h + wc)
        pl_A = simpson(lambda s: (2.0 * (s - h) - wc) * phi(s), h + wc, up)
        pl_n = simpson(phi, h + wc, up)
        sb = self.sigma * self.beta
        P_el = (4.0 / 3.0) * self.E_star * self.eta * sb * math.sqrt(self.sigma / self.beta) * el_P
        A_el = math.pi * self.eta * sb * el_A
        A_pl = math.pi * self.eta * sb * pl_A
        P_pl = self.H * A_pl
        return dict(h=h, P=P_el + P_pl, P_el=P_el, P_pl=P_pl,
                    A=A_el + A_pl, A_el=A_el, A_pl=A_pl,
                    n=self.eta * (el_n + pl_n), n_pl=self.eta * pl_n)

    def _build_table(self):
        self.table = []
        h = 10.0
        while h >= -4.0 - 1e-9:
            st = self._state(h)
            if st["P"] > 0:
                self.table.append(st)
            h -= 0.05
        # table: h 감소 순 = P 증가 순

    def _interp(self, key, target, field="P"):
        """table 에서 field==target 이 되는 지점의 모든 값을 보간."""
        tb = self.table
        vals = [t[field] for t in tb]
        if target <= vals[0]:
            return dict(tb[0]), "below"
        if target >= vals[-1]:
            return dict(tb[-1]), "above"
        for i in range(len(tb) - 1):
            v0, v1 = vals[i], vals[i + 1]
            if v0 <= target <= v1:
                f = (target - v0) / (v1 - v0) if v1 > v0 else 0.0
                out = {k: tb[i][k] + f * (tb[i + 1][k] - tb[i][k]) for k in tb[i]}
                return out, "ok"
        return dict(tb[-1]), "above"

    def at_pressure(self, P):
        st, flag = self._interp("P", P)
        st["flag"] = flag
        st["f_pl"] = st["P_pl"] / st["P"] if st["P"] > 0 else 0.0
        st["p_a"] = st["P"] / st["A"] if st["A"] > 0 else float("nan")
        return st

    def critical_pressure(self, frac=0.5):
        """소성 돌기가 분담하는 하중 비율이 frac 에 도달하는 공칭압력 [Pa]."""
        tb = self.table
        fr = [t["P_pl"] / t["P"] for t in tb]
        if fr[0] >= frac:
            return 0.0, "always"           # 극저압에서도 이미 소성 (ψ 지배)
        if fr[-1] < frac:
            return float("inf"), "never"   # 표 범위 내 압력에서는 탄성 유지
        for i in range(len(tb) - 1):
            if fr[i] <= frac <= fr[i + 1]:
                f = (frac - fr[i]) / (fr[i + 1] - fr[i])
                return tb[i]["P"] + f * (tb[i + 1]["P"] - tb[i]["P"]), "ok"
        return float("inf"), "never"


# ================================================================== 입자 수준 계산
def compute(inp):
    """inp: dict (SI 변환 전 GUI 값). 반환: dict(결과) — 모두 SI."""
    r = {}

    # ---- Pad 물성
    E_sh, path = hardness_to_modulus(inp["hard_val"], inp["hard_scale"])
    r["E_shore_MPa"], r["hard_path"] = E_sh, path
    nu_p = inp["nu_pad"]
    E_pad = inp["E_pad"] * 1e6 if inp["E_pad"] > 0 else E_sh * 1e6
    r["E_pad_from_input"] = inp["E_pad"] > 0
    src = inp["H_source"]
    if src == "직접 입력 (MPa)":
        H_pad = inp["H_pad_direct"] * 1e6
    elif src == "입력 Modulus × (H/E)":
        H_pad = inp["kH"] * E_pad
    else:  # 경도 환산 E × (H/E)
        H_pad = inp["kH"] * E_sh * 1e6
    r["E_pad"], r["H_pad"] = E_pad, H_pad

    # ---- 막질 / 입자
    E_w, H_w, nu_w, K_w = inp["E_w"] * 1e9, inp["H_w"] * 1e9, inp["nu_w"], inp["K_w"]
    rho_p, E_p, H_p, nu_p_ = inp["rho_p"] * 1e3, inp["E_p"] * 1e9, inp["H_p"] * 1e9, inp["nu_p"]
    d = inp["d_nm"] * 1e-9
    R = d / 2.0

    E_star_pw = plane_strain_E(E_pad, nu_p)                 # Pad–Wafer (wafer 강체 가정)
    E_star_pp = contact_modulus(E_pad, nu_p, E_p, nu_p_)    # Pad–입자
    E_star_wp = contact_modulus(E_w, nu_w, E_p, nu_p_)      # Wafer–입자
    r.update(E_star_pw=E_star_pw, E_star_pp=E_star_pp, E_star_wp=E_star_wp)

    # ---- Pad 돌기 접촉
    P = inp["P_val"] * UNIT_TO_PA[inp["P_unit"]]
    r["P"] = P
    pad = PadContact(E_star_pw, H_pad, inp["sigma_um"] * 1e-6, inp["beta_um"] * 1e-6,
                     inp["eta_mm2"] * 1e6, nu_p)
    st = pad.at_pressure(P)
    r["pad"] = pad
    r["psi"] = pad.psi
    r["omega_c"] = pad.omega_c
    r["st"] = st
    r["P_crit"], r["P_crit_flag"] = pad.critical_pressure(0.5)
    r["P_10"], _ = pad.critical_pressure(0.1)
    r["P_90"], _ = pad.critical_pressure(0.9)
    plastic = st["f_pl"] >= 0.5
    r["pad_plastic"] = plastic

    # ---- 슬러리 → 입자 수밀도
    w = inp["wt"] / 100.0
    rho_f = inp["rho_f"] * 1e3
    phi_v = (w / rho_p) / (w / rho_p + (1 - w) / rho_f) if 0 < w < 1 else 0.0
    V_p = math.pi * d ** 3 / 6.0
    n_v = phi_v / V_p if V_p > 0 else 0.0
    n_a = n_v ** (2.0 / 3.0) * inp["active_frac"]
    r.update(phi_v=phi_v, n_v=n_v, n_a=n_a,
             L_bulk=n_v ** (-1.0 / 3.0) if n_v > 0 else float("inf"),
             L_areal=n_a ** (-0.5) if n_a > 0 else float("inf"))

    # ---- 입자당 하중
    p_a = st["p_a"]
    F_share = p_a / n_a if n_a > 0 else float("inf")
    if plastic:
        F_emb = H_pad * math.pi * R * R
        emb_label = "Pad 소성 한계  F = H_pad·πR²"
    else:
        F_emb = (4.0 / 3.0) * E_star_pp * R * R
        emb_label = "Pad 탄성 완전매립 한계  F = (4/3)E*·R²"
    F_p = min(F_share, F_emb)
    r.update(p_a=p_a, F_share=F_share, F_emb=F_emb, F_p=F_p, emb_label=emb_label,
             hardness_controlled=F_emb <= F_share)
    r["N_active_per_cm2"] = n_a * st["A"] * 1e-4

    # Pad 매립 깊이
    if plastic:
        x = F_p / (math.pi * H_pad)
        delta_p = R - math.sqrt(max(R * R - x, 0.0)) if x < R * R else R
    else:
        delta_p = min((9.0 * F_p ** 2 / (16.0 * E_star_pp ** 2 * R)) ** (1.0 / 3.0), R)
    r["delta_p"] = delta_p

    # ---- 입자–Wafer 접촉 (Hertz → 소성 캡)
    a = (3.0 * F_p * R / (4.0 * E_star_wp)) ** (1.0 / 3.0)
    p_m = F_p / (math.pi * a * a)
    Y_w = H_w / 2.8
    p_yield = 1.1 * Y_w
    if p_m >= H_w:
        wafer_mode = "완전 소성 (p_m ≥ H_w)"
        a = math.sqrt(F_p / (math.pi * H_w))
        p_m = H_w
        delta_w = R - math.sqrt(max(R * R - a * a, 0.0))
    elif p_m >= p_yield:
        wafer_mode = "항복 시작 (1.1·Y_w ≤ p_m < H_w)"
        delta_w = a * a / R
    else:
        wafer_mode = "탄성 (p_m < 1.1·Y_w)"
        delta_w = a * a / R
    r.update(a_w=a, p_m=p_m, p_yield=p_yield, wafer_mode=wafer_mode, delta_w=delta_w)

    # 경도 지배 영역에서의 wafer 접촉압력 (입자 크기 무관)
    p_w_hc = (16.0 * H_pad * E_star_wp ** 2 / (9.0 * math.pi ** 2)) ** (1.0 / 3.0)
    H_pad_crit_y = 9.0 * math.pi ** 2 * p_yield ** 3 / (16.0 * E_star_wp ** 2)
    H_pad_crit_H = 9.0 * math.pi ** 2 * H_w ** 3 / (16.0 * E_star_wp ** 2)
    r.update(p_w_hc=p_w_hc, H_pad_crit_y=H_pad_crit_y, H_pad_crit_H=H_pad_crit_H)

    # ---- 균열/스크래치 임계하중 (Lawn)
    if K_w is not None and K_w > 0:
        Kc = K_w * 1e6
        P_sharp = 2.2e4 * Kc ** 4 / H_w ** 3                       # Lawn–Evans 1977
        P_blunt = inp["A_auer"] * Kc ** 2 * R / E_star_wp          # Auerbach / Frank–Lawn
        if plastic:
            d_c_sharp = 2.0 * math.sqrt(P_sharp / (math.pi * H_pad))
            d_c_blunt = 2.0 * inp["A_auer"] * Kc ** 2 / (math.pi * H_pad * E_star_wp)
        else:
            d_c_sharp = 2.0 * math.sqrt(3.0 * P_sharp / (4.0 * E_star_pp))
            d_c_blunt = 2.0 * 3.0 * inp["A_auer"] * Kc ** 2 / (4.0 * E_star_pp * E_star_wp)
        r.update(P_sharp=P_sharp, P_blunt=P_blunt, d_c_sharp=d_c_sharp, d_c_blunt=d_c_blunt)
    else:
        r.update(P_sharp=None, P_blunt=None, d_c_sharp=None, d_c_blunt=None)

    # ---- 압력 sweep (그래프용): F_p, f_pl vs P
    sweep = []
    for i in range(0, 46):
        Ps = 10 ** (-3.0 + 4.5 * i / 45.0) * PSI       # 0.001 ~ 30 psi
        s2 = pad.at_pressure(Ps)
        pl2 = s2["f_pl"] >= 0.5
        Fe = H_pad * math.pi * R * R if pl2 else (4.0 / 3.0) * E_star_pp * R * R
        Fs = s2["p_a"] / n_a if n_a > 0 else float("inf")
        sweep.append((Ps, s2["f_pl"], min(Fs, Fe)))
    r["sweep"] = sweep
    return r


# ================================================================== 결과 서식
def fmt_P(Pa, unit):
    if Pa == float("inf"):
        return "∞"
    return f"{Pa / UNIT_TO_PA[unit]:.4g} {unit}  ({Pa / PSI:.4g} psi)"


def fmt_len(m):
    if m == float("inf"):
        return "∞"
    if m >= 1e-3:
        return f"{m * 1e3:.3g} mm"
    if m >= 1e-6:
        return f"{m * 1e6:.3g} μm"
    return f"{m * 1e9:.3g} nm"


def fmt_F(N):
    if N == float("inf"):
        return "∞"
    if N >= 1e-3:
        return f"{N * 1e3:.3g} mN"
    if N >= 1e-6:
        return f"{N * 1e6:.3g} μN"
    return f"{N * 1e9:.3g} nN"


def fmt_p(Pa):
    if Pa >= 1e9:
        return f"{Pa / 1e9:.3g} GPa"
    if Pa >= 1e6:
        return f"{Pa / 1e6:.3g} MPa"
    return f"{Pa / 1e3:.3g} kPa"


def regime_label(psi):
    if psi < 0.6:
        return "탄성 (ψ < 0.6)"
    if psi <= 1.0:
        return "천이 (0.6 ≤ ψ ≤ 1.0)"
    return "소성 (ψ > 1.0)"


def build_report(inp, r):
    u = inp["P_unit"]
    st = r["st"]
    L = []
    A = L.append
    A("=" * 66)
    A(" [1] PAD 돌기 접촉 판정")
    A("=" * 66)
    A(f" 경도 환산 탄성계수     : {r['E_shore_MPa']:.3g} MPa   ({r['hard_path']})")
    A(f" 계산에 쓴 Pad E        : {r['E_pad'] / 1e6:.3g} MPa   "
      f"({'입력 Modulus' if r['E_pad_from_input'] else '경도 환산값 사용'})")
    ratio = r["E_pad"] / (r["E_shore_MPa"] * 1e6)
    if r["E_pad_from_input"] and (ratio > 3 or ratio < 1 / 3):
        A(f"   ⚠ 경도 환산 E 와 입력 Modulus 가 {max(ratio, 1 / ratio):.1f} 배 차이 — 경도 척도·측정 조건을 확인하세요.")
    A(f" Pad 경도 H_pad         : {fmt_p(r['H_pad'])}   (산출 기준: {inp['H_source']}, H/E = {inp['kH']})")
    A(f" Pad–Wafer 접촉계수 E*  : {fmt_p(r['E_star_pw'])}")
    A(f" 돌기 σ / β / η         : {inp['sigma_um']} μm / {inp['beta_um']} μm / {inp['eta_mm2']} /mm²")
    A("")
    A(f" ▶ 소성지수 ψ = (E*/H)·√(σ/β) = {r['psi']:.3g}   →  {regime_label(r['psi'])}")
    A(f"   임계 간섭 ω_c         : {fmt_len(r['omega_c'])}  (ω_c/σ = {r['pad'].wc:.3g})")
    A("")
    A(f" 공칭압력 P              : {fmt_P(r['P'], u)}")
    if st["flag"] == "above":
        A("   ⚠ 압력이 모델 표 범위를 넘어 최대값으로 고정됨")
    A(f"   실접촉면적비 A_r/A_n  : {st['A']:.3g}  ({st['A'] * 100:.3g} %)")
    A(f"   돌기 평균 접촉압력 p_a: {fmt_p(st['p_a'])}")
    A(f"   접촉 돌기 수          : {st['n'] * 1e-6:.3g} /mm²  (그중 소성 {st['n_pl'] * 1e-6:.3g} /mm²)")
    A(f"   소성 돌기 하중 분담률 : {st['f_pl'] * 100:.1f} %")
    A("")
    pc, flag = r["P_crit"], r["P_crit_flag"]
    if flag == "always":
        A(" ▶ 임계압력 P_crit(Pad)  : ≈ 0  — 극저압에서도 돌기 접촉이 소성 (ψ 지배). 모든 공정 압력에서 소성.")
    elif flag == "never":
        A(" ▶ 임계압력 P_crit(Pad)  : 모델 범위(수백 psi) 안에서 소성 전이 없음 — 탄성 접촉 유지.")
    else:
        A(f" ▶ 임계압력 P_crit(Pad)  : {fmt_P(pc, u)}   (소성 돌기 하중분담 50 % 기준)")
        A(f"   전이 구간 10 %~90 %    : {fmt_P(r['P_10'], u)}  ~  {fmt_P(r['P_90'], u)}")
    above = r["P"] >= pc
    A(f" ▶ 판정: 현재 압력은 임계압력 {'이상' if above else '이하'} → Pad 돌기 접촉 "
      f"{'소성' if above else '탄성'}  ⇒ 입자당 하중은 "
      f"{'Pad 경도(H_pad)에 의해 결정 (Head 압력 무관)' if above else 'Pad 탄성(E*)과 압력에 의해 결정'}")
    A("")
    A("=" * 66)
    A(" [2] 슬러리 · 입자당 하중")
    A("=" * 66)
    A(f" 입자 {inp['particle']}  d = {inp['d_nm']} nm, ρ = {inp['rho_p']} g/cm³")
    A(f" Solid contents {inp['wt']} wt%  →  부피분율 φ = {r['phi_v'] * 100:.3g} vol%")
    A(f" 입자 수밀도 n_v         : {r['n_v'] * 1e-18:.3g} ×10¹⁸ /m³   평균 간격(벌크) {fmt_len(r['L_bulk'])}")
    A(f" 계면 단층 면밀도 n_a    : {r['n_a'] * 1e-12:.3g} /μm²        평균 간격(면)   {fmt_len(r['L_areal'])}")
    A(f" 활성 입자 수 (실접촉 내): {r['N_active_per_cm2']:.3g} /cm²")
    A("")
    A(f" Head 압력 분담 하중 F_share = p_a / n_a      : {fmt_F(r['F_share'])}")
    A(f" Pad 매립 한계 하중 F_emb                     : {fmt_F(r['F_emb'])}   ({r['emb_label']})")
    A(f" ▶ 입자당 하중 F_p = min(F_share, F_emb)       : {fmt_F(r['F_p'])}")
    if r["hardness_controlled"]:
        A("   → F_emb 가 지배: 입자가 Pad 에 완전 매립. 하중이 Head 압력이 아닌 Pad 경도로 결정됨.")
    else:
        A("   → F_share 가 지배: 입자가 부분 매립. 하중이 Head 압력·입자 농도에 따라 변함.")
    A(f"   Pad 매립 깊이 δ_p      : {fmt_len(r['delta_p'])}  (R = {inp['d_nm'] / 2:.3g} nm)")
    A("")
    A("=" * 66)
    A(" [3] 입자–WAFER 접촉")
    A("=" * 66)
    A(f" 막질 {inp['film']}: E = {inp['E_w']} GPa, H = {inp['H_w']} GPa, ν = {inp['nu_w']}, "
      f"K_IC = {inp['K_w'] if inp['K_w'] else '—(연성)'} MPa·√m")
    A(f" Wafer–입자 접촉계수 E*  : {fmt_p(r['E_star_wp'])}")
    A(f" 접촉 반경 a            : {fmt_len(r['a_w'])}")
    A(f" ▶ 입자당 접촉압력 p_m   : {fmt_p(r['p_m'])}   vs  항복 시작 1.1·Y_w = {fmt_p(r['p_yield'])},  H_w = {fmt_p(inp['H_w'] * 1e9)}")
    A(f" ▶ 막질 변형 모드        : {r['wafer_mode']}")
    A(f"   Wafer 압입 깊이 δ_w    : {fmt_len(r['delta_w'])}")
    A("")
    A(f" 경도 지배 영역의 p_m (입자 크기 무관) = (16·H_pad·E*²/9π²)^(1/3) = {fmt_p(r['p_w_hc'])}")
    A(f"   막질 항복 시작에 필요한 Pad 경도 : {fmt_p(r['H_pad_crit_y'])}  /  완전 소성: {fmt_p(r['H_pad_crit_H'])}")
    A("")
    A("=" * 66)
    A(" [4] 스크래치 · 균열 임계하중 (Lawn)")
    A("=" * 66)
    if r["P_sharp"] is None:
        A(" 막질이 연성(K_IC 미정의) — 균열 기준 해당 없음. 소성 압입(스크래치 홈) 기준은 [3] 참조.")
    else:
        A(f" 날카로운 압입 균열 시작  P* = 2.2×10⁴·K_IC⁴/H³ : {fmt_F(r['P_sharp'])}")
        A(f" 둥근 압입 Hertz 콘균열   P_c = A·K_IC²·R/E*    : {fmt_F(r['P_blunt'])}   (A = {inp['A_auer']:.3g})")
        A(f" ▶ 현재 F_p / P*  = {r['F_p'] / r['P_sharp']:.3g}      F_p / P_c = {r['F_p'] / r['P_blunt']:.3g}")
        ok = r["F_p"] < min(r["P_sharp"], r["P_blunt"])
        A(f" ▶ 판정: {'균열 임계 이하 — 단일 입자로는 균열/스크래치 불가' if ok else '균열 임계 이상 — 스크래치 위험'}")
        A(f"   균열을 일으키는 임계 입자(응집체) 지름: 날카로운 기준 {fmt_len(r['d_c_sharp'])}, "
          f"둥근 기준 {fmt_len(r['d_c_blunt'])}")
    A("")
    A(" 가정: GW/CEB 돌기 통계접촉, 계면 단층 입자 n_a = n_v^(2/3)·활성비, Hertz 구-평면 접촉,")
    A("       Pad 완전매립 한계 = H_pad·πR² (소성) 또는 (4/3)E*R² (탄성). 상세는 README 참조.")
    return "\n".join(L)


# ================================================================== GUI
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1280x820")
        self.minsize(1100, 700)
        self.vars = {}
        self._build()
        self.on_particle_select()
        self.on_film_select()
        self.on_compute()
        self.after(150, self._nudge_redraw)

    def _nudge_redraw(self):
        self.update_idletasks()
        w, h = self.winfo_width(), self.winfo_height()
        x, y = self.winfo_x(), self.winfo_y()
        self.geometry(f"{w + 1}x{h + 1}+{x}+{y}")
        self.after(50, lambda: self.geometry(f"{w}x{h}+{x}+{y}"))

    # ---------- 위젯 도우미
    def _entry(self, parent, row, label, key, default, width=10, unit=""):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=1)
        v = tk.StringVar(value=str(default))
        self.vars[key] = v
        ttk.Entry(parent, textvariable=v, width=width).grid(row=row, column=1, sticky="w", pady=1)
        if unit:
            ttk.Label(parent, text=unit, foreground="#666").grid(row=row, column=2, sticky="w")
        return v

    def _combo(self, parent, row, label, key, values, default, cmd=None, width=16):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=1)
        v = tk.StringVar(value=default)
        self.vars[key] = v
        cb = ttk.Combobox(parent, textvariable=v, values=values, state="readonly", width=width)
        cb.grid(row=row, column=1, columnspan=2, sticky="w", pady=1)
        if cmd:
            cb.bind("<<ComboboxSelected>>", lambda e: cmd())
        return v

    def _build(self):
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ---- Pad
        f = ttk.LabelFrame(left, text="Pad", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._combo(f, 0, "경도 척도", "hard_scale", ["Shore A", "Shore D", "Asker C"], "Shore D")
        self._entry(f, 1, "경도 값", "hard_val", 57)
        self._entry(f, 2, "Modulus (RT)", "E_pad", 300, unit="MPa  (0 = 경도 환산값 사용)")
        self._entry(f, 3, "포아송비 ν", "nu_pad", 0.35)
        self._combo(f, 4, "H_pad 산출", "H_source",
                    ["경도 환산 E × (H/E)", "입력 Modulus × (H/E)", "직접 입력 (MPa)"], "경도 환산 E × (H/E)",
                    width=22)
        self._entry(f, 5, "H/E 비", "kH", 0.1)
        self._entry(f, 6, "H_pad 직접 입력", "H_pad_direct", 30, unit="MPa")

        # ---- 거칠기
        f = ttk.LabelFrame(left, text="Pad 돌기 거칠기 (GW 파라미터)", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._entry(f, 0, "σ 돌기 높이 표준편차", "sigma_um", 5.0, unit="μm")
        self._entry(f, 1, "β 돌기 끝단 반경", "beta_um", 30.0, unit="μm")
        self._entry(f, 2, "η 돌기 밀도", "eta_mm2", 100.0, unit="/mm²")
        ttk.Separator(f, orient="horizontal").grid(row=3, column=0, columnspan=3, sticky="ew", pady=4)
        ttk.Label(f, text="AFM / ISO 25178 측정값으로 환산", foreground="#444").grid(
            row=4, column=0, columnspan=3, sticky="w")
        self._entry(f, 5, "Sq (= Rq)", "Sq_um", 5.0, unit="μm   → σ")
        self._entry(f, 6, "Spc 평균 봉우리 곡률", "Spc", 0.0333, unit="1/μm → β = 1/Spc")
        self._entry(f, 7, "Spd 봉우리 밀도", "Spd", 100.0, unit="/mm² → η")
        ttk.Button(f, text="환산값 적용 (σ, β, η 채우기)", command=self.on_apply_afm).grid(
            row=8, column=0, columnspan=3, sticky="ew", pady=(2, 0))

        # ---- 압력
        f = ttk.LabelFrame(left, text="공정 압력", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._entry(f, 0, "Head 압력", "P_val", 3.0)
        self._combo(f, 1, "단위", "P_unit", ["psi", "hPa", "kPa"], "psi", width=8)

        # ---- 슬러리
        f = ttk.LabelFrame(left, text="슬러리 · 입자", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._entry(f, 0, "Solid contents", "wt", 5.0, unit="wt%")
        self._entry(f, 1, "액상 밀도", "rho_f", 1.0, unit="g/cm³")
        self._combo(f, 2, "입자 재질", "particle", list(PARTICLES.keys()), "실리카 (SiO2)",
                    cmd=self.on_particle_select)
        self._entry(f, 3, "입자 지름", "d_nm", 100.0, unit="nm")
        self._entry(f, 4, "입자 밀도", "rho_p", 2.2, unit="g/cm³")
        self._entry(f, 5, "입자 E", "E_p", 70.0, unit="GPa")
        self._entry(f, 6, "입자 H", "H_p", 8.0, unit="GPa")
        self._entry(f, 7, "입자 ν", "nu_p", 0.17)
        self._entry(f, 8, "활성 입자 비율", "active_frac", 1.0, unit="(계면 단층 중 하중 분담 비율)")

        # ---- 막질
        f = ttk.LabelFrame(left, text="Wafer 막질", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._combo(f, 0, "막질", "film", list(FILMS.keys()), "SiO2 (산화막)", cmd=self.on_film_select)
        self._entry(f, 1, "막질 E", "E_w", 70.0, unit="GPa")
        self._entry(f, 2, "막질 H", "H_w", 8.0, unit="GPa")
        self._entry(f, 3, "막질 ν", "nu_w", 0.17)
        self._entry(f, 4, "K_IC", "K_w", 0.75, unit="MPa·√m  (0 = 연성, 균열기준 생략)")
        self._entry(f, 5, "Auerbach 상수 A", "A_auer", 8600, unit="(Hertz 콘균열)")

        # ---- 버튼
        b = ttk.Frame(left)
        b.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(b, text="계산", command=self.on_compute).pack(fill=tk.X, pady=2)
        ttk.Button(b, text="결과 저장 (TXT)", command=self.on_save_txt).pack(fill=tk.X, pady=2)
        ttk.Button(b, text="조건 저장 (JSON)", command=self.on_save_json).pack(fill=tk.X, pady=2)
        ttk.Button(b, text="조건 불러오기 (JSON)", command=self.on_load_json).pack(fill=tk.X, pady=2)

        # ---- 오른쪽: 결과 텍스트 + 그래프
        self.summary = tk.StringVar(value="")
        ttk.Label(right, textvariable=self.summary, font=("Consolas", 11, "bold"),
                  justify="left").pack(anchor="w", pady=(0, 4))

        self.txt = tk.Text(right, font=("Consolas", 10), wrap="none", height=28)
        self.txt.pack(fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(self.txt, orient="vertical", command=self.txt.yview)
        self.txt.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        self.canvas = tk.Canvas(right, height=230, bg="white", highlightthickness=1,
                                highlightbackground="#bbb")
        self.canvas.pack(fill=tk.X, pady=(6, 0))
        self.canvas.bind("<Configure>", lambda e: self._draw_plot())
        self.result = None

    # ---------- 이벤트
    def on_particle_select(self):
        p = PARTICLES[self.vars["particle"].get()]
        if p:
            for k, v in zip(("rho_p", "E_p", "H_p", "nu_p"), p[:4]):
                self.vars[k].set(str(v))

    def on_film_select(self):
        p = FILMS[self.vars["film"].get()]
        if p:
            for k, v in zip(("E_w", "H_w", "nu_w"), p[:3]):
                self.vars[k].set(str(v))
            self.vars["K_w"].set(str(p[3] if p[3] is not None else 0))

    def on_apply_afm(self):
        try:
            sq, spc, spd = (float(self.vars[k].get()) for k in ("Sq_um", "Spc", "Spd"))
            if spc <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("입력 오류", "Sq, Spc(>0), Spd 는 숫자여야 합니다.")
            return
        self.vars["sigma_um"].set(f"{sq:g}")
        self.vars["beta_um"].set(f"{1.0 / spc:g}")
        self.vars["eta_mm2"].set(f"{spd:g}")
        self.on_compute()

    def _read_inputs(self):
        num_keys = ["hard_val", "E_pad", "nu_pad", "kH", "H_pad_direct", "sigma_um", "beta_um",
                    "eta_mm2", "P_val", "wt", "rho_f", "d_nm", "rho_p", "E_p", "H_p", "nu_p",
                    "active_frac", "E_w", "H_w", "nu_w", "K_w", "A_auer"]
        inp = {}
        for k in num_keys:
            try:
                inp[k] = float(self.vars[k].get())
            except ValueError:
                raise ValueError(f"'{k}' 입력이 숫자가 아닙니다.")
        for k in ["hard_scale", "H_source", "P_unit", "particle", "film"]:
            inp[k] = self.vars[k].get()
        if inp["K_w"] <= 0:
            inp["K_w"] = None
        for k in ["sigma_um", "beta_um", "eta_mm2", "d_nm", "rho_p", "E_w", "H_w", "E_p"]:
            if inp[k] <= 0:
                raise ValueError(f"'{k}' 는 0보다 커야 합니다.")
        if not (0 <= inp["wt"] < 100):
            raise ValueError("Solid contents 는 0 이상 100 미만이어야 합니다.")
        return inp

    def on_compute(self):
        try:
            inp = self._read_inputs()
            r = compute(inp)
        except ValueError as e:
            messagebox.showerror("입력 오류", str(e))
            return
        except Exception as e:  # 계산 중 예외는 메시지로
            messagebox.showerror("계산 오류", repr(e))
            return
        self.inp, self.result = inp, r
        u = inp["P_unit"]
        pc = r["P_crit"]
        pc_txt = ("≈ 0 (항상 소성)" if r["P_crit_flag"] == "always"
                  else "없음 (탄성 유지)" if r["P_crit_flag"] == "never"
                  else f"{pc / UNIT_TO_PA[u]:.3g} {u}")
        above = r["P"] >= pc
        self.summary.set(
            f"소성지수 ψ = {r['psi']:.3g}  [{regime_label(r['psi'])}]      "
            f"임계압력 P_crit = {pc_txt}      "
            f"현재 P = {inp['P_val']:g} {u} → {'임계 이상 (경도 지배)' if above else '임계 이하 (탄성·압력 지배)'}\n"
            f"입자당 하중 F_p = {fmt_F(r['F_p'])}   Wafer 접촉압력 p_m = {fmt_p(r['p_m'])}  [{r['wafer_mode']}]"
        )
        self.txt.configure(state="normal")
        self.txt.delete("1.0", tk.END)
        self.txt.insert(tk.END, build_report(inp, r))
        self.txt.configure(state="disabled")
        self._draw_plot()

    # ---------- 그래프 (tk Canvas, log-x)
    def _draw_plot(self):
        c = self.canvas
        c.delete("all")
        if not self.result:
            return
        W, H = c.winfo_width(), c.winfo_height()
        if W < 50 or H < 50:
            return
        r, inp = self.result, self.inp
        sw = r["sweep"]
        ml, mr, mt, mb = 70, 70, 18, 34
        x0, x1 = ml, W - mr
        y0, y1 = mt, H - mb
        lx_min, lx_max = -3.0, math.log10(30.0)

        def X(P):
            return x0 + (math.log10(max(P / PSI, 1e-9)) - lx_min) / (lx_max - lx_min) * (x1 - x0)

        Fs = [f for _, _, f in sw if f != float("inf")]
        Fmax = max(Fs) * 1.15 if Fs else 1.0

        def YF(F):
            return y1 - (F / Fmax) * (y1 - y0)

        def Yf(fr):
            return y1 - fr * (y1 - y0)

        # 축
        c.create_rectangle(x0, y0, x1, y1, outline="#888")
        for e in range(-3, 2):
            xx = X(10 ** e * PSI)
            c.create_line(xx, y1, xx, y1 + 4, fill="#444")
            c.create_text(xx, y1 + 12, text=f"{10 ** e:g}", font=("Consolas", 8))
        c.create_text((x0 + x1) / 2, H - 6, text="Head 압력 P [psi, log]", font=("Consolas", 9))
        for k in range(0, 5):
            yy = Yf(k / 4)
            c.create_text(x0 - 6, yy, text=f"{k * 25}%", anchor="e", font=("Consolas", 8), fill="#1f77b4")
            c.create_text(x1 + 6, yy, text=fmt_F(Fmax * k / 4), anchor="w", font=("Consolas", 8), fill="#d62728")
        c.create_text(x0 - 6, y0 - 8, text="소성 분담률", anchor="e", font=("Consolas", 8), fill="#1f77b4")
        c.create_text(x1 + 6, y0 - 8, text="입자당 하중 F_p", anchor="w", font=("Consolas", 8), fill="#d62728")

        # 곡선
        pts_f = [(X(P), Yf(fr)) for P, fr, _ in sw]
        pts_F = [(X(P), YF(F)) for P, _, F in sw if F != float("inf")]
        if len(pts_f) > 1:
            c.create_line(*sum(pts_f, ()), fill="#1f77b4", width=2)
        if len(pts_F) > 1:
            c.create_line(*sum(pts_F, ()), fill="#d62728", width=2)

        # 임계압력 · 현재 압력
        pc = r["P_crit"]
        if 0 < pc < float("inf"):
            xx = X(pc)
            c.create_line(xx, y0, xx, y1, fill="#2ca02c", dash=(4, 3), width=2)
            c.create_text(xx + 3, y0 + 8, text="P_crit", anchor="w", fill="#2ca02c", font=("Consolas", 8))
        xx = X(r["P"])
        c.create_line(xx, y0, xx, y1, fill="#555", dash=(2, 2))
        c.create_text(xx + 3, y1 - 8, text="현재 P", anchor="w", fill="#555", font=("Consolas", 8))

    # ---------- 저장/불러오기
    def on_save_txt(self):
        if not self.result:
            return
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text", "*.txt")],
                                            initialfile="cmp_particle_pressure_result.txt")
        if path:
            with open(path, "w", encoding="utf-8-sig") as f:
                f.write(APP_TITLE + "\n\n" + self.summary.get() + "\n\n" + build_report(self.inp, self.result))

    def on_save_json(self):
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json")],
                                            initialfile="cmp_particle_pressure_inputs.json")
        if path:
            with open(path, "w", encoding="utf-8") as f:
                json.dump({k: v.get() for k, v in self.vars.items()}, f, ensure_ascii=False, indent=2)

    def on_load_json(self):
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json")])
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
        except Exception as e:
            messagebox.showerror("불러오기 오류", repr(e))
            return
        for k, v in data.items():
            if k in self.vars:
                self.vars[k].set(str(v))
        self.on_compute()


if __name__ == "__main__":
    App().mainloop()
