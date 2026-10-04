# -*- coding: utf-8 -*-
"""
CMP Particle Pressure Simulator  R02  (pad asperity contact + single-abrasive contact)
======================================================================================
Pad 경도(Shore A / Shore D / Asker C), Pad Bulk Modulus, Head 압력, 슬러리 Solid contents를
입력받아 다음을 계산한다.

  [1] Pad 돌기(asperity) 접촉이 탄성인지 소성인지
        - Bulk Modulus → 돌기 끝단 Modulus 추정 (기공률 보정 + 수화 보정)
        - 유효 Pad 경도 H_eff = (H/E)·E_tip × f_rate (속도·글레이징 경화) 또는 직접 입력
        - 소성지수 ψ (Greenwood–Williamson), 압력에 따른 소성 전이 (CEB) → 임계압력 P_crit
  [2] 입자 하나가 받는 하중 F_p, 입자간 간격, 크기 분포 상위 분위수, LPC 응집체 접촉 수
        - Head 압력 분담 하중  vs  Pad–입자–Wafer 하중 평형(매립 깊이)으로 구한 Pad 한계 하중
  [3] 입자–Wafer 접촉압력 (Hertz), 매립 깊이 δ_p / δ_w
  [4] 스크래치 판정
        - 1차: 수화층(t_h, H_h) 관통 → 나노스크래치 (폭 출력)
        - 2차: 벌크 막질(H_w) 잔류 압흔 깊이 ≥ δ_th → 깊은 스크래치
        - 임계 응집체 지름, 상위 분위수 입자·LPC 응집체 판정

모든 계산은 파이썬 표준 라이브러리만 사용한다 (tkinter 포함).
    python cmp_particle_pressure_sim_R02.py

모델/상수의 출처와 가정은 README.md 참조.
"""

import json
import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

APP_TITLE = "CMP Particle Pressure Simulator  R02"

# ------------------------------------------------------------------ 단위
PSI = 6894.757                        # Pa
UNIT_TO_PA = {"psi": PSI, "hPa": 100.0, "kPa": 1000.0}
WAFER_R = 0.150                       # 300 mm 웨이퍼 반경 [m] (LPC 접촉 수 추정용)

# ------------------------------------------------------------------ 재료 기본값
# 입자: (밀도 g/cm3, E GPa, H GPa, ν)  — 대표값. 측정값으로 교체 권장.
PARTICLES = {
    "실리카 (SiO2)":    (2.20,  70.0,  8.0, 0.17),
    "세리아 (CeO2)":    (7.20, 180.0,  7.0, 0.30),
    "알루미나 (Al2O3)": (3.95, 380.0, 20.0, 0.23),
    "직접 입력":        None,
}
# 막질: (E GPa, H GPa, ν, 수화층 두께 nm, 수화층 경도 GPa)
FILMS = {
    "SiO2 (산화막)": ( 70.0,  8.0, 0.17, 2.0, 1.0),
    "Cu":            (120.0,  1.5, 0.34, 2.0, 0.5),
    "W":             (410.0,  8.0, 0.28, 2.0, 1.0),
    "poly-Si":       (160.0, 12.0, 0.22, 2.0, 1.0),
    "Si3N4":         (260.0, 19.0, 0.25, 2.0, 1.0),
    "직접 입력":     None,
}

# Asker C ↔ Shore A 근사 대응표 (듀로미터 비교표 기반, 근사치)
ASKER_C_TABLE = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
SHORE_A_TABLE = [ 0,  8, 18, 28, 40, 50, 62, 75, 88, 100]

H_SRC_AUTO = "E_tip × (H/E) × f_rate"
H_SRC_DIRECT = "직접 입력 (글레이징·나노압입 값, MPa)"


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
    """Gent (1958): Shore A → Young's modulus [MPa]."""
    s = min(max(s, 1.0), 99.0)
    return 0.0981 * (56 + 7.62336 * s) / (0.137505 * (254 - 2.54 * s))


def shore_to_E_ruess(s, scale):
    """Ruess 상관식: log10 E[MPa] = 0.0235·S − 0.6403 (Shore A), S→S+50 (Shore D)."""
    if scale == "Shore D":
        s = s + 50.0
    return 10 ** (0.0235 * s - 0.6403)


def hardness_to_modulus(value, scale):
    if scale == "Shore A":
        return shore_a_to_E_gent(value), f"Shore A {value:.0f} → Gent 식"
    if scale == "Shore D":
        return shore_to_E_ruess(value, "Shore D"), f"Shore D {value:.0f} → Ruess 식"
    if scale == "Asker C":
        sa = asker_c_to_shore_a(value)
        return shore_a_to_E_gent(sa), f"Asker C {value:.0f} ≈ Shore A {sa:.0f} → Gent 식"
    raise ValueError(scale)


# ================================================================== 거칠기 환산 (GW ↔ ISO 25178)
def gw_to_afm(sigma_um, beta_um, eta_mm2):
    return sigma_um * 1e3, 1e3 / beta_um, eta_mm2


def afm_to_gw(Sq_nm, Spc_per_mm, Spd_mm2):
    return Sq_nm / 1e3, 1e3 / Spc_per_mm, Spd_mm2


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


def bisect(f, lo, hi, it=80):
    """f(lo) 와 f(hi) 부호가 다를 때 근을 찾는다."""
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def plane_strain_E(E, nu):
    return E / (1.0 - nu * nu)


def contact_modulus(E1, nu1, E2, nu2):
    return 1.0 / ((1 - nu1 ** 2) / E1 + (1 - nu2 ** 2) / E2)


def cap_area(R, d):
    """반경 R 구가 깊이 d 만큼 박힐 때 투영 접촉면적 (d > R 이면 πR² 로 포화)."""
    d = max(d, 0.0)
    if d >= R:
        return math.pi * R * R
    return math.pi * (2.0 * R * d - d * d)


def cap_depth(R, A):
    """투영면적 A 에 해당하는 캡 깊이 (A ≥ πR² 이면 R)."""
    a2 = A / math.pi
    if a2 >= R * R:
        return R
    return R - math.sqrt(R * R - a2)


# ================================================================== Pad 돌기 접촉 (GW + CEB)
class PadContact:
    """Greenwood–Williamson 통계 접촉 + Chang–Etsion–Bogy 탄소성 돌기 모델."""

    def __init__(self, E_star, H, sigma, beta, eta, nu):
        self.E_star, self.H = E_star, H
        self.sigma, self.beta, self.eta = sigma, beta, eta
        K = 0.454 + 0.41 * nu
        self.omega_c = (math.pi * K * H / (2.0 * E_star)) ** 2 * beta
        self.wc = self.omega_c / sigma
        self.psi = (E_star / H) * math.sqrt(sigma / beta)
        self._build_table()

    def _state(self, h):
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

    def _interp(self, target, field="P"):
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
                return {k: tb[i][k] + f * (tb[i + 1][k] - tb[i][k]) for k in tb[i]}, "ok"
        return dict(tb[-1]), "above"

    def at_pressure(self, P):
        st, flag = self._interp(P)
        st["flag"] = flag
        st["f_pl"] = st["P_pl"] / st["P"] if st["P"] > 0 else 0.0
        st["p_a"] = st["P"] / st["A"] if st["A"] > 0 else float("nan")
        return st

    def critical_pressure(self, frac=0.5):
        tb = self.table
        fr = [t["P_pl"] / t["P"] for t in tb]
        if fr[0] >= frac:
            return 0.0, "always"
        if fr[-1] < frac:
            return float("inf"), "never"
        for i in range(len(tb) - 1):
            if fr[i] <= frac <= fr[i + 1]:
                f = (frac - fr[i]) / (fr[i + 1] - fr[i])
                return tb[i]["P"] + f * (tb[i + 1]["P"] - tb[i]["P"]), "ok"
        return float("inf"), "never"


# ================================================================== 입자 하중 평형 (Pad ↔ Wafer)
def particle_balance(R, pad_plastic, H_eff, E_pp, E_wp, gap=0.0):
    """Pad–입자–Wafer 하중 평형.

    Pad 반력  F_pad(δ_p) = H_eff·A_cap(δ_p)            (소성)  또는 (4/3)E*_pp √R δ_p^1.5 (탄성, δ_p ≤ R)
    Wafer 반력 F_waf(δ_w) = (4/3)E*_wp √R δ_w^1.5        (Hertz)
    기하 닫힘  δ_p + δ_w = 2R − gap   (gap: 입자 위치의 Pad 표면–Wafer 간극. 소성 돌기 접촉 내부는 0)
    반환: F, δ_p, δ_w
    """
    span = 2.0 * R - gap
    if span <= 0:
        return 0.0, 0.0, 0.0

    def F_pad(dp):
        if pad_plastic:
            return H_eff * cap_area(R, dp)
        d = min(dp, R)
        return (4.0 / 3.0) * E_pp * math.sqrt(R) * d ** 1.5

    def F_waf(dw):
        return (4.0 / 3.0) * E_wp * math.sqrt(R) * max(dw, 0.0) ** 1.5

    g = lambda dp: F_pad(dp) - F_waf(span - dp)
    dp = bisect(g, 1e-15, span - 1e-15)
    F = F_pad(dp)
    return F, dp, span - dp


# ================================================================== 메인 계산
def compute(inp):
    r = {}

    # ---- Pad: Bulk E → 돌기 끝단 E_tip → 유효 경도 H_eff
    E_sh, path = hardness_to_modulus(inp["hard_val"], inp["hard_scale"])
    r["E_shore_MPa"], r["hard_path"] = E_sh, path
    nu_p = inp["nu_pad"]
    E_bulk = inp["E_pad"] * 1e6 if inp["E_pad"] > 0 else E_sh * 1e6
    r["E_bulk"], r["E_bulk_from_input"] = E_bulk, inp["E_pad"] > 0

    p = inp["porosity"] / 100.0
    C_tip = inp["f_wet"] / (1.0 - p) ** inp["n_ga"]
    if inp["E_tip_direct"] > 0:
        E_tip, tip_path = inp["E_tip_direct"] * 1e6, "직접 입력 (나노압입 등)"
    else:
        E_tip = E_bulk * C_tip
        tip_path = (f"E_bulk / (1−p)^n × f_wet = ×{C_tip:.3g}  "
                    f"(p = {inp['porosity']:g} %, n = {inp['n_ga']:g}, f_wet = {inp['f_wet']:g})")
    r.update(E_tip=E_tip, tip_path=tip_path, C_tip=C_tip)

    # 속도·경화 보정: 고분자는 변형률 속도가 오르면 E 와 H 가 함께 오르므로 둘 다 f_rate 로 스케일한다.
    # 직접 입력(글레이징) 모드에서는 경도와 모순되지 않도록 E ≥ H_eff/(H/E) 로 끌어올린다.
    H_static = inp["kH"] * E_tip
    if inp["H_source"] == H_SRC_DIRECT:
        H_eff = inp["H_pad_direct"] * 1e6
        E_mech = max(E_tip, H_eff / inp["kH"])
        H_note = "직접 입력 (글레이징·나노압입 값)"
        E_note = ("= E_tip" if E_mech == E_tip
                  else f"= H_eff/(H/E) — 입력 H_eff 와 일관되도록 E_tip({fmt_p(E_tip)}) 에서 상향")
    else:
        E_mech = E_tip * inp["f_rate"]
        H_eff = H_static * inp["f_rate"]
        H_note = f"(H/E = {inp['kH']:g}) × E_tip × f_rate {inp['f_rate']:g}"
        E_note = f"= E_tip × f_rate {inp['f_rate']:g}"
    r.update(H_static=H_static, H_eff=H_eff, H_note=H_note, E_mech=E_mech, E_note=E_note)

    # ---- 막질 / 입자
    E_w, H_w, nu_w = inp["E_w"] * 1e9, inp["H_w"] * 1e9, inp["nu_w"]
    t_h, H_h = inp["t_h_nm"] * 1e-9, inp["H_h"] * 1e9
    rho_p, E_p, nu_pp = inp["rho_p"] * 1e3, inp["E_p"] * 1e9, inp["nu_p"]
    d = inp["d_nm"] * 1e-9
    R = d / 2.0

    E_star_pw = plane_strain_E(E_mech, nu_p)
    E_star_pp = contact_modulus(E_mech, nu_p, E_p, nu_pp)
    E_star_wp = contact_modulus(E_w, nu_w, E_p, nu_pp)
    r.update(E_star_pw=E_star_pw, E_star_pp=E_star_pp, E_star_wp=E_star_wp)

    # ---- 거칠기
    if inp["rough_mode"] == "AFM":
        sigma_um, beta_um, eta_mm2 = afm_to_gw(inp["Sq_nm"], inp["Spc"], inp["Spd"])
    else:
        sigma_um, beta_um, eta_mm2 = inp["sigma_um"], inp["beta_um"], inp["eta_mm2"]
    r.update(sigma_um=sigma_um, beta_um=beta_um, eta_mm2=eta_mm2)

    # ---- Pad 돌기 접촉
    P = inp["P_val"] * UNIT_TO_PA[inp["P_unit"]]
    r["P"] = P
    pad = PadContact(E_star_pw, H_eff, sigma_um * 1e-6, beta_um * 1e-6, eta_mm2 * 1e6, nu_p)
    st = pad.at_pressure(P)
    r.update(pad=pad, psi=pad.psi, omega_c=pad.omega_c, st=st)
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

    # ---- 입자당 하중: 분담 하중 vs Pad–Wafer 하중 평형
    p_a = st["p_a"]
    F_share = p_a / n_a if n_a > 0 else float("inf")
    F_bal, dp_bal, dw_bal = particle_balance(R, plastic, H_eff, E_star_pp, E_star_wp)
    k_emb = F_bal / (R * R)                       # 경도 지배 하중 계수 F = k·R²
    if F_bal <= F_share:
        F_p, hardness_controlled = F_bal, True
        delta_p, delta_w_el = dp_bal, dw_bal
    else:
        F_p, hardness_controlled = F_share, False
        delta_p = cap_depth(R, F_p / H_eff) if plastic else \
            min((9.0 * F_p ** 2 / (16.0 * E_star_pp ** 2 * R)) ** (1.0 / 3.0), R)
        delta_w_el = (9.0 * F_p ** 2 / (16.0 * E_star_wp ** 2 * R)) ** (1.0 / 3.0)
    r.update(p_a=p_a, F_share=F_share, F_bal=F_bal, F_p=F_p, k_emb=k_emb,
             hardness_controlled=hardness_controlled, delta_p=delta_p,
             emb_label=("Pad 소성 반력 H_eff·A_cap(δ_p)" if plastic else "Pad 탄성 Hertz 반력"))
    r["N_active_per_cm2"] = n_a * st["A"] * 1e-4

    # ---- 입자–Wafer 접촉 (Hertz, 벌크 기준)
    a = (3.0 * F_p * R / (4.0 * E_star_wp)) ** (1.0 / 3.0)
    p_m = F_p / (math.pi * a * a)
    Y_w = H_w / 2.8
    p_yield = 1.1 * Y_w
    if p_m >= H_w:
        wafer_mode = "완전 소성 (p_m ≥ H_w)"
    elif p_m >= p_yield:
        wafer_mode = "항복 시작 (1.1·Y_w ≤ p_m < H_w)"
    else:
        wafer_mode = "탄성 (p_m < 1.1·Y_w)"
    r.update(a_w=a, p_m=p_m, p_yield=p_yield, wafer_mode=wafer_mode, delta_w=delta_w_el)

    p_w_hc = (16.0 * k_emb * E_star_wp ** 2 / (9.0 * math.pi ** 3)) ** (1.0 / 3.0)
    H_pad_crit_y = 9.0 * math.pi ** 2 * p_yield ** 3 / (16.0 * E_star_wp ** 2)
    r.update(p_w_hc=p_w_hc, H_pad_crit_y=H_pad_crit_y)

    # ---- 스크래치 1차: 수화층 관통 (나노스크래치)
    def layer_pen(F, Rx):
        """수화층 경도 H_h 기준 압흔 반경·깊이."""
        a2 = F / (math.pi * H_h)
        if a2 >= Rx * Rx:
            return Rx, Rx
        return math.sqrt(a2), Rx - math.sqrt(Rx * Rx - a2)

    p_y_h = 1.1 * H_h / 2.8
    if p_m >= p_y_h:
        a_h, dep_h = layer_pen(F_p, R)
    else:
        a_h, dep_h = 0.0, 0.0
    nano = dep_h >= t_h
    r.update(a_h=a_h, dep_h=dep_h, nano=nano, p_y_h=p_y_h)
    # 수화층 관통에 필요한 하중 (R 기준)
    A_need = cap_area(R, t_h)
    F_need = H_h * A_need
    r.update(F_need=F_need, p_need=F_need / (math.pi * R * R), w_need=2.0 * math.sqrt(A_need / math.pi))
    # 임계 응집체 지름 (경도 지배 F = k·R²): δ_h/R = 1 − √(1 − k/πH_h)
    ratio_h = k_emb / (math.pi * H_h)
    if p_w_hc < p_y_h:
        d_c_nano, d_c_nano_note = None, "경도 지배 접촉압력이 수화층 항복 미만"
    elif ratio_h >= 1.0:
        d_c_nano, d_c_nano_note = 2.0 * t_h, "Pad 하중 한계가 수화층 경도 초과 (δ = R)"
    else:
        g_h = 1.0 - math.sqrt(1.0 - ratio_h)
        d_c_nano, d_c_nano_note = 2.0 * t_h / g_h, f"δ_h/R = {g_h:.3g}"
    r.update(d_c_nano=d_c_nano, d_c_nano_note=d_c_nano_note)

    # ---- 스크래치 2차: 벌크 잔류 압흔 (깊은 스크래치)
    d_th = inp["scratch_nm"] * 1e-9
    if p_m >= p_yield:
        a2 = F_p / (math.pi * H_w)
        a_res = min(math.sqrt(a2), R)
        dep_res = R - math.sqrt(max(R * R - a2, 0.0)) if a2 < R * R else R
    else:
        a_res, dep_res = 0.0, 0.0
    deep = dep_res >= d_th
    ratio_w = k_emb / (math.pi * H_w)
    if p_w_hc < p_yield:
        d_c_deep, d_c_deep_note = None, "경도 지배 접촉압력이 벌크 항복 미만"
    elif ratio_w >= 1.0:
        d_c_deep, d_c_deep_note = 2.0 * d_th, "Pad 하중 한계가 벌크 경도 초과 (δ = R)"
    else:
        g_w = 1.0 - math.sqrt(1.0 - ratio_w)
        d_c_deep, d_c_deep_note = 2.0 * d_th / g_w, f"δ_res/R = {g_w:.3g}"
    r.update(a_res=a_res, dep_res=dep_res, d_th=d_th, deep=deep,
             d_c_deep=d_c_deep, d_c_deep_note=d_c_deep_note)

    # ---- 입자 크기 분포 상위 분위수 (로그정규, CV)
    cv = inp["cv_pct"] / 100.0
    s_ln = math.sqrt(math.log(1.0 + cv * cv)) if cv > 0 else 0.0
    d_999 = d * math.exp(3.0902 * s_ln)           # 99.9 % 분위 (z = 3.0902)
    R_999 = d_999 / 2.0
    F_999 = min(F_share, k_emb * R_999 * R_999)
    a_999 = (3.0 * F_999 * R_999 / (4.0 * E_star_wp)) ** (1.0 / 3.0)
    p_999 = F_999 / (math.pi * a_999 * a_999)
    if p_999 >= p_y_h:
        a_h999, dep_h999 = layer_pen(F_999, R_999)
    else:
        a_h999, dep_h999 = 0.0, 0.0
    r.update(d_999=d_999, F_999=F_999, dep_h999=dep_h999, a_h999=a_h999, nano_999=dep_h999 >= t_h)

    # ---- LPC 응집체
    d_L = inp["d_L_um"] * 1e-6
    lpc = inp["lpc"] * 1e6                        # counts/mL → /m³
    n_aL = lpc ** (2.0 / 3.0) if lpc > 0 else 0.0
    N_L = n_aL * st["A"] * math.pi * WAFER_R ** 2
    r.update(d_L=d_L, N_L=N_L,
             L_nano=(d_c_nano is not None and d_L >= d_c_nano),
             L_deep=(d_c_deep is not None and d_L >= d_c_deep))

    # ---- 압력 sweep (그래프용)
    sweep = []
    for i in range(0, 46):
        Ps = 10 ** (-3.0 + 4.5 * i / 45.0) * PSI
        s2 = pad.at_pressure(Ps)
        pl2 = s2["f_pl"] >= 0.5
        Fb, _, _ = particle_balance(R, pl2, H_eff, E_star_pp, E_star_wp)
        Fs = s2["p_a"] / n_a if n_a > 0 else float("inf")
        sweep.append((Ps, s2["f_pl"], min(Fs, Fb)))
    r["sweep"] = sweep
    return r


# ================================================================== 결과 서식
def fmt_P(Pa, unit):
    if Pa == float("inf"):
        return "∞"
    return f"{Pa / UNIT_TO_PA[unit]:.4g} {unit}  ({Pa / PSI:.4g} psi)"


def fmt_len(m):
    if m is None:
        return "—"
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
    A("=" * 72)
    A(" [1] PAD 돌기 접촉 판정")
    A("=" * 72)
    A(f" 경도 환산 탄성계수 (bulk) : {r['E_shore_MPa']:.3g} MPa   ({r['hard_path']})")
    A(f" Bulk Modulus E_bulk        : {r['E_bulk'] / 1e6:.3g} MPa   "
      f"({'입력 Modulus' if r['E_bulk_from_input'] else '경도 환산값 사용'})")
    ratio = r["E_bulk"] / (r["E_shore_MPa"] * 1e6)
    if r["E_bulk_from_input"] and (ratio > 3 or ratio < 1 / 3):
        A(f"   ⚠ 경도 환산 E 와 입력 Modulus 가 {max(ratio, 1 / ratio):.1f} 배 차이 — 경도 척도·측정 조건을 확인하세요.")
    A(f" ▶ 돌기 끝단 Modulus E_tip  : {fmt_p(r['E_tip'])}   [{r['tip_path']}]")
    A(f"   정적 Pad 경도 H_static   : {fmt_p(r['H_static'])}   (= E_tip × H/E {inp['kH']:g})")
    A(f" ▶ 유효 Pad 경도 H_eff      : {fmt_p(r['H_eff'])}   [{r['H_note']}]")
    A(f" ▶ 접촉역학용 Pad E_dyn     : {fmt_p(r['E_mech'])}   [{r['E_note']}]")
    A(f" Pad–Wafer 접촉계수 E*      : {fmt_p(r['E_star_pw'])}  (E_dyn 기준)")
    A(f" 돌기 σ / β / η             : {r['sigma_um']:.3g} μm / {r['beta_um']:.3g} μm / {r['eta_mm2']:.3g} /mm²"
      f"   [{'AFM 환산' if inp['rough_mode'] == 'AFM' else 'GW 직접'}]")
    A("")
    A(f" ▶ 소성지수 ψ = (E*/H_eff)·√(σ/β) = {r['psi']:.3g}   →  {regime_label(r['psi'])}")
    A(f"   임계 간섭 ω_c             : {fmt_len(r['omega_c'])}  (ω_c/σ = {r['pad'].wc:.3g})")
    A("")
    A(f" 공칭압력 P                  : {fmt_P(r['P'], u)}")
    if st["flag"] == "above":
        A("   ⚠ 압력이 모델 표 범위를 넘어 최대값으로 고정됨")
    A(f"   실접촉면적비 A_r/A_n      : {st['A']:.3g}  ({st['A'] * 100:.3g} %)")
    A(f"   돌기 평균 접촉압력 p_a    : {fmt_p(st['p_a'])}   (소성 접촉에서는 정의상 p_a = H_eff)")
    A(f"   접촉 돌기 수              : {st['n'] * 1e-6:.3g} /mm²  (그중 소성 {st['n_pl'] * 1e-6:.3g} /mm²)")
    A(f"   소성 돌기 하중 분담률     : {st['f_pl'] * 100:.1f} %")
    A("")
    pc, flag = r["P_crit"], r["P_crit_flag"]
    if flag == "always":
        A(" ▶ 임계압력 P_crit(Pad)      : ≈ 0  — 극저압에서도 돌기 접촉이 소성 (ψ 지배). 모든 공정 압력에서 소성.")
    elif flag == "never":
        A(" ▶ 임계압력 P_crit(Pad)      : 모델 범위(수백 psi) 안에서 소성 전이 없음 — 탄성 접촉 유지.")
    else:
        A(f" ▶ 임계압력 P_crit(Pad)      : {fmt_P(pc, u)}   (소성 돌기 하중분담 50 % 기준)")
        A(f"   전이 구간 10 %~90 %        : {fmt_P(r['P_10'], u)}  ~  {fmt_P(r['P_90'], u)}")
    above = r["P"] >= pc
    A(f" ▶ 판정: 현재 압력은 임계압력 {'이상' if above else '이하'} → Pad 돌기 접촉 "
      f"{'소성' if above else '탄성'}  ⇒ 입자당 하중은 "
      f"{'Pad 경도(H_eff)에 의해 결정 (Head 압력 무관)' if above else 'Pad 탄성(E*)과 압력에 의해 결정'}")
    A("")
    A("=" * 72)
    A(" [2] 슬러리 · 입자당 하중")
    A("=" * 72)
    A(f" 입자 {inp['particle']}  d = {inp['d_nm']:g} nm (중앙값), ρ = {inp['rho_p']:g} g/cm³, 크기 분포 CV {inp['cv_pct']:g} %")
    A(f" Solid contents {inp['wt']:g} wt%  →  부피분율 φ = {r['phi_v'] * 100:.3g} vol%")
    A(f" 입자 수밀도 n_v             : {r['n_v'] * 1e-18:.3g} ×10¹⁸ /m³   평균 간격(벌크) {fmt_len(r['L_bulk'])}")
    A(f" 계면 활성 면밀도 n_a        : {r['n_a'] * 1e-12:.3g} /μm²   (활성비 {inp['active_frac']:g})   평균 간격(면) {fmt_len(r['L_areal'])}")
    A(f" 활성 입자 수 (실접촉 내)    : {r['N_active_per_cm2']:.3g} /cm²")
    A("")
    A(f" Head 압력 분담 하중 F_share = p_a / n_a            : {fmt_F(r['F_share'])}")
    A(f" Pad–Wafer 하중 평형 하중 F_bal (δ_p + δ_w = 2R)     : {fmt_F(r['F_bal'])}   ({r['emb_label']}, k = F/R² = {r['k_emb']:.3g} Pa)")
    A(f" ▶ 입자당 하중 F_p = min(F_share, F_bal)              : {fmt_F(r['F_p'])}")
    if r["hardness_controlled"]:
        A("   → F_bal 이 지배: 입자가 Pad 에 매립되어 하중이 Head 압력이 아닌 Pad 경도로 결정됨.")
    else:
        A("   → F_share 가 지배: 입자가 부분 매립. 하중이 Head 압력·입자 농도에 따라 변함.")
    A(f"   Pad 매립 깊이 δ_p          : {fmt_len(r['delta_p'])}  (R = {inp['d_nm'] / 2:.3g} nm{'; δ_p ≥ R 이면 투영면적 πR² 로 포화' if r['delta_p'] >= 0.999 * inp['d_nm'] / 2 * 1e-9 else ''})")
    A(f" 상위 99.9 % 입자 d_99.9     : {fmt_len(r['d_999'])}  →  F_p = {fmt_F(r['F_999'])}")
    A("")
    A("=" * 72)
    A(" [3] 입자–WAFER 접촉 (벌크 기준)")
    A("=" * 72)
    A(f" 막질 {inp['film']}: E = {inp['E_w']:g} GPa, H = {inp['H_w']:g} GPa, ν = {inp['nu_w']:g}  |  수화층 t_h = {inp['t_h_nm']:g} nm, H_h = {inp['H_h']:g} GPa")
    A(f" Wafer–입자 접촉계수 E*      : {fmt_p(r['E_star_wp'])}")
    A(f" Hertz 접촉 반경 a          : {fmt_len(r['a_w'])}    탄성 압입 깊이 δ_w : {fmt_len(r['delta_w'])}")
    A(f" ▶ 입자당 접촉압력 p_m       : {fmt_p(r['p_m'])}   vs  벌크 항복 시작 1.1·Y_w = {fmt_p(r['p_yield'])},  H_w = {fmt_p(inp['H_w'] * 1e9)}")
    A(f" ▶ 벌크 변형 모드            : {r['wafer_mode']}")
    A(f" 경도 지배 영역의 p_m (크기 무관) = (16·k·E*²/9π³)^(1/3) = {fmt_p(r['p_w_hc'])}")
    A(f"   벌크 항복 시작에 필요한 Pad 경도 H_eff : {fmt_p(r['H_pad_crit_y'])}")
    A("")
    A("=" * 72)
    A(" [4] 스크래치 판정")
    A("=" * 72)
    A(f" ① 나노스크래치 — 수화층 관통 기준 (t_h = {inp['t_h_nm']:g} nm, H_h = {inp['H_h']:g} GPa)")
    A(f"    관통에 필요한 하중 (R = {inp['d_nm'] / 2:g} nm) : F_need = {fmt_F(r['F_need'])},  Pad 국소 압력 F/πR² = {fmt_p(r['p_need'])},  홈 폭 = {fmt_len(r['w_need'])}")
    if r["dep_h"] <= 0:
        A("    단일 입자 수화층 압흔             : 없음 (접촉압력이 수화층 항복 미만)")
    else:
        A(f"    단일 입자 수화층 압흔             : 깊이 {fmt_len(r['dep_h'])},  폭 {fmt_len(2 * r['a_h'])}")
    A(f"    ▶ 판정 (중앙값 입자)              : {'나노스크래치 발생 (δ_h ≥ t_h)' if r['nano'] else '불가 (δ_h < t_h)'}")
    tail999 = "" if r["dep_h999"] <= 0 else f" — 깊이 {fmt_len(r['dep_h999'])}, 폭 {fmt_len(2 * r['a_h999'])}"
    A(f"    ▶ 판정 (상위 99.9 % 입자 {fmt_len(r['d_999'])}) : "
      f"{'나노스크래치 발생' if r['nano_999'] else '불가'}{tail999}")
    if r["d_c_nano"] is None:
        A(f"    ▶ 임계 응집체 지름 d_c(나노)        : 없음 — {r['d_c_nano_note']}")
    else:
        A(f"    ▶ 임계 응집체 지름 d_c(나노)        : {fmt_len(r['d_c_nano'])}   ({r['d_c_nano_note']})")
    A("")
    A(f" ② 깊은 스크래치 — 벌크 잔류 압흔 기준 (δ_th = {fmt_len(r['d_th'])}, H_w = {inp['H_w']:g} GPa)")
    if r["dep_res"] <= 0:
        A("    단일 입자 벌크 잔류 압흔          : 없음 (접촉압력이 벌크 항복 미만 → 탄성 회복)")
    else:
        A(f"    단일 입자 벌크 잔류 압흔          : 깊이 {fmt_len(r['dep_res'])},  폭 {fmt_len(2 * r['a_res'])}")
    A(f"    ▶ 판정                            : {'깊은 스크래치 발생 (δ_res ≥ δ_th)' if r['deep'] else '불가 (δ_res < δ_th)'}")
    if r["d_c_deep"] is None:
        A(f"    ▶ 임계 응집체 지름 d_c(깊은)        : 없음 — {r['d_c_deep_note']}")
    else:
        A(f"    ▶ 임계 응집체 지름 d_c(깊은)        : {fmt_len(r['d_c_deep'])}   ({r['d_c_deep_note']})")
    A("")
    A(f" ③ LPC 응집체 — d_L = {fmt_len(r['d_L'])}, LPC = {inp['lpc']:g} counts/mL")
    A(f"    실접촉 영역 내 순간 접촉 수 (300 mm 웨이퍼, 추정) : {r['N_L']:.3g} 개")
    A(f"    ▶ d_L vs d_c(나노) : {'나노스크래치 유발 크기' if r['L_nano'] else '미만'}      "
      f"d_L vs d_c(깊은) : {'깊은 스크래치 유발 크기' if r['L_deep'] else '미만'}")
    A("")
    A(" 가정: GW/CEB 돌기 통계접촉(E_dyn, H_eff 기준), 계면 활성 입자 n_a = n_v^(2/3)·활성비, Hertz 구-평면 접촉,")
    A("       입자 하중 평형 δ_p + δ_w = 2R (소성 돌기 접촉 내부 간극 0), 압흔 = 경도 정의 H = F/πa². 상세는 README 참조.")
    return "\n".join(L)


# ================================================================== GUI
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1380x940")
        self.minsize(1200, 760)
        self.vars = {}
        self.result = None
        self._build()
        self.on_particle_select()
        self.on_film_select()
        self.on_rough_mode()
        self.on_compute()
        self.after(150, self._nudge_redraw)

    def _nudge_redraw(self):
        self.update_idletasks()
        w, h = self.winfo_width(), self.winfo_height()
        x, y = self.winfo_x(), self.winfo_y()
        self.geometry(f"{w + 1}x{h + 1}+{x}+{y}")
        self.after(50, lambda: self.geometry(f"{w}x{h}+{x}+{y}"))

    def _entry(self, parent, row, label, key, default, width=10, unit="", col=0):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", pady=1)
        v = tk.StringVar(value=str(default))
        self.vars[key] = v
        e = ttk.Entry(parent, textvariable=v, width=width)
        e.grid(row=row, column=col + 1, sticky="w", pady=1)
        if unit:
            ttk.Label(parent, text=unit, foreground="#666").grid(row=row, column=col + 2, sticky="w")
        return e

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
        # 왼쪽 입력 패널을 두 열로 나눈다
        left = ttk.Frame(self, padding=8)
        left.pack(side=tk.LEFT, fill=tk.Y)
        colA = ttk.Frame(left)
        colA.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 6))
        colB = ttk.Frame(left)
        colB.pack(side=tk.LEFT, fill=tk.Y)
        right = ttk.Frame(self, padding=(0, 8, 8, 8))
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ---- Pad
        f = ttk.LabelFrame(colA, text="Pad", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._combo(f, 0, "경도 척도", "hard_scale", ["Shore A", "Shore D", "Asker C"], "Shore D")
        self._entry(f, 1, "경도 값", "hard_val", 57)
        self._entry(f, 2, "Bulk Modulus (RT)", "E_pad", 300, unit="MPa  (0 = 경도 환산값)")
        self._entry(f, 3, "포아송비 ν", "nu_pad", 0.35)
        ttk.Separator(f, orient="horizontal").grid(row=4, column=0, columnspan=3, sticky="ew", pady=4)
        ttk.Label(f, text="돌기 끝단 Modulus  E_tip = E_bulk / (1−p)^n × f_wet", foreground="#444").grid(
            row=5, column=0, columnspan=3, sticky="w")
        self._entry(f, 6, "기공률 p", "porosity", 33, unit="vol%")
        self._entry(f, 7, "Gibson–Ashby 지수 n", "n_ga", 2.0, unit="(개방 2, 폐쇄 1~2)")
        self._entry(f, 8, "수화·온도 보정 f_wet", "f_wet", 1.0, unit="(건조 1.0)")
        self._entry(f, 9, "E_tip 직접 입력", "E_tip_direct", 0, unit="MPa  (0 = 자동)")
        ttk.Separator(f, orient="horizontal").grid(row=10, column=0, columnspan=3, sticky="ew", pady=4)
        ttk.Label(f, text="유효 Pad 경도  H_eff", foreground="#444").grid(row=11, column=0, columnspan=3, sticky="w")
        self._combo(f, 12, "H_eff 산출", "H_source", [H_SRC_AUTO, H_SRC_DIRECT], H_SRC_AUTO, width=30)
        self._entry(f, 13, "H/E 비", "kH", 0.1)
        self._entry(f, 14, "동적 보정 f_rate", "f_rate", 5.0, unit="(정적 1, 슬라이딩 1~20; E·H 동시 적용)")
        self._entry(f, 15, "H_eff 직접 입력", "H_pad_direct", 100, unit="MPa  (글레이징 표면 등)")

        # ---- 거칠기
        f = ttk.LabelFrame(colA, text="Pad 돌기 거칠기", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self.vars["rough_mode"] = tk.StringVar(value="GW")
        mf = ttk.Frame(f)
        mf.grid(row=0, column=0, columnspan=5, sticky="w", pady=(0, 4))
        ttk.Label(mf, text="입력 방식:").pack(side=tk.LEFT)
        ttk.Radiobutton(mf, text="GW 직접", variable=self.vars["rough_mode"], value="GW",
                        command=self.on_rough_mode).pack(side=tk.LEFT, padx=(6, 10))
        ttk.Radiobutton(mf, text="AFM / ISO 25178", variable=self.vars["rough_mode"], value="AFM",
                        command=self.on_rough_mode).pack(side=tk.LEFT)
        ttk.Label(f, text="GW", foreground="#1f77b4").grid(row=1, column=1, sticky="w")
        ttk.Label(f, text="AFM", foreground="#d62728").grid(row=1, column=3, sticky="w")
        rows = [
            ("σ 돌기 높이 표준편차 [μm]", "sigma_um", 5.0, "Sq_nm", 5000.0, "Sq [nm]"),
            ("β 돌기 끝단 반경 [μm]",     "beta_um", 30.0, "Spc", 33.333, "Spc [1/mm] (β = 1/Spc)"),
            ("η 돌기 밀도 [/mm²]",        "eta_mm2", 100.0, "Spd", 100.0, "Spd [/mm²]"),
        ]
        self.gw_entries, self.afm_entries = [], []
        for i, (gl, gk, gd, ak, ad, al) in enumerate(rows, start=2):
            e1 = self._entry(f, i, gl, gk, gd, width=9)
            ttk.Label(f, text="⇄").grid(row=i, column=2, padx=4)
            v = tk.StringVar(value=str(ad))
            self.vars[ak] = v
            e2 = ttk.Entry(f, textvariable=v, width=9)
            e2.grid(row=i, column=3, sticky="w", pady=1)
            ttk.Label(f, text=al, foreground="#666").grid(row=i, column=4, sticky="w")
            self.gw_entries.append(e1)
            self.afm_entries.append(e2)
        ttk.Label(f, text="비활성 쪽에는 활성 쪽 입력의 환산값이 자동 표시됩니다.", foreground="#888").grid(
            row=5, column=0, columnspan=5, sticky="w", pady=(3, 0))

        # ---- 압력
        f = ttk.LabelFrame(colA, text="공정 압력", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._entry(f, 0, "Head 압력", "P_val", 3.0)
        self._combo(f, 1, "단위", "P_unit", ["psi", "hPa", "kPa"], "psi", width=8)

        # ---- 슬러리
        f = ttk.LabelFrame(colB, text="슬러리 · 입자", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._entry(f, 0, "Solid contents", "wt", 5.0, unit="wt%")
        self._entry(f, 1, "액상 밀도", "rho_f", 1.0, unit="g/cm³")
        self._combo(f, 2, "입자 재질", "particle", list(PARTICLES.keys()), "실리카 (SiO2)",
                    cmd=self.on_particle_select)
        self._entry(f, 3, "입자 지름 (중앙값)", "d_nm", 100.0, unit="nm")
        self._entry(f, 4, "크기 분포 CV", "cv_pct", 20.0, unit="%  (로그정규, 99.9 % 분위 계산)")
        self._entry(f, 5, "입자 밀도", "rho_p", 2.2, unit="g/cm³")
        self._entry(f, 6, "입자 E", "E_p", 70.0, unit="GPa")
        self._entry(f, 7, "입자 H", "H_p", 8.0, unit="GPa")
        self._entry(f, 8, "입자 ν", "nu_p", 0.17)
        self._entry(f, 9, "활성 입자 비율", "active_frac", 0.05, unit="(계면 단층 중 하중 분담 비율)")
        ttk.Separator(f, orient="horizontal").grid(row=10, column=0, columnspan=3, sticky="ew", pady=4)
        ttk.Label(f, text="LPC (대입자/응집체)", foreground="#444").grid(row=11, column=0, columnspan=3, sticky="w")
        self._entry(f, 12, "응집체 지름 d_L", "d_L_um", 0.5, unit="μm")
        self._entry(f, 13, "LPC", "lpc", 1000.0, unit="counts/mL")

        # ---- 막질
        f = ttk.LabelFrame(colB, text="Wafer 막질 · 스크래치 기준", padding=6)
        f.pack(fill=tk.X, pady=(0, 6))
        self._combo(f, 0, "막질", "film", list(FILMS.keys()), "SiO2 (산화막)", cmd=self.on_film_select)
        self._entry(f, 1, "막질 E", "E_w", 70.0, unit="GPa")
        self._entry(f, 2, "막질 H (벌크)", "H_w", 8.0, unit="GPa")
        self._entry(f, 3, "막질 ν", "nu_w", 0.17)
        ttk.Separator(f, orient="horizontal").grid(row=4, column=0, columnspan=3, sticky="ew", pady=4)
        self._entry(f, 5, "수화층 두께 t_h", "t_h_nm", 2.0, unit="nm  (① 나노스크래치 기준)")
        self._entry(f, 6, "수화층 경도 H_h", "H_h", 1.0, unit="GPa")
        self._entry(f, 7, "깊은 스크래치 깊이 δ_th", "scratch_nm", 10.0, unit="nm  (② 벌크 잔류 압흔 기준)")

        # ---- 버튼
        b = ttk.Frame(colB)
        b.pack(fill=tk.X, pady=(4, 0))
        ttk.Button(b, text="계산", command=self.on_compute).pack(fill=tk.X, pady=2)
        ttk.Button(b, text="결과 저장 (TXT)", command=self.on_save_txt).pack(fill=tk.X, pady=2)
        ttk.Button(b, text="조건 저장 (JSON)", command=self.on_save_json).pack(fill=tk.X, pady=2)
        ttk.Button(b, text="조건 불러오기 (JSON)", command=self.on_load_json).pack(fill=tk.X, pady=2)

        # ---- 오른쪽
        self.summary = tk.StringVar(value="")
        ttk.Label(right, textvariable=self.summary, font=("Consolas", 11, "bold"),
                  justify="left").pack(anchor="w", pady=(0, 4))
        self.txt = tk.Text(right, font=("Consolas", 10), wrap="none", height=30)
        self.txt.pack(fill=tk.BOTH, expand=True)
        sb = ttk.Scrollbar(self.txt, orient="vertical", command=self.txt.yview)
        self.txt.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas = tk.Canvas(right, height=220, bg="white", highlightthickness=1,
                                highlightbackground="#bbb")
        self.canvas.pack(fill=tk.X, pady=(6, 0))
        self.canvas.bind("<Configure>", lambda e: self._draw_plot())

    # ---------- 이벤트
    def on_particle_select(self):
        p = PARTICLES[self.vars["particle"].get()]
        if p:
            for k, v in zip(("rho_p", "E_p", "H_p", "nu_p"), p):
                self.vars[k].set(str(v))

    def on_film_select(self):
        p = FILMS[self.vars["film"].get()]
        if p:
            for k, v in zip(("E_w", "H_w", "nu_w", "t_h_nm", "H_h"), p):
                self.vars[k].set(str(v))

    def on_rough_mode(self):
        afm = self.vars["rough_mode"].get() == "AFM"
        for e in self.gw_entries:
            e.configure(state="disabled" if afm else "normal")
        for e in self.afm_entries:
            e.configure(state="normal" if afm else "disabled")
        self._sync_roughness()

    def _sync_roughness(self):
        try:
            if self.vars["rough_mode"].get() == "AFM":
                sq, spc, spd = (float(self.vars[k].get()) for k in ("Sq_nm", "Spc", "Spd"))
                s, b, n = afm_to_gw(sq, spc, spd)
                self.vars["sigma_um"].set(f"{s:.4g}")
                self.vars["beta_um"].set(f"{b:.4g}")
                self.vars["eta_mm2"].set(f"{n:.4g}")
            else:
                s, b, n = (float(self.vars[k].get()) for k in ("sigma_um", "beta_um", "eta_mm2"))
                sq, spc, spd = gw_to_afm(s, b, n)
                self.vars["Sq_nm"].set(f"{sq:.4g}")
                self.vars["Spc"].set(f"{spc:.4g}")
                self.vars["Spd"].set(f"{spd:.4g}")
        except (ValueError, ZeroDivisionError):
            pass

    def _read_inputs(self):
        num_keys = ["hard_val", "E_pad", "nu_pad", "porosity", "n_ga", "f_wet", "E_tip_direct", "kH",
                    "f_rate", "H_pad_direct", "sigma_um", "beta_um", "eta_mm2", "Sq_nm", "Spc", "Spd",
                    "P_val", "wt", "rho_f", "d_nm", "cv_pct", "rho_p", "E_p", "H_p", "nu_p", "active_frac",
                    "d_L_um", "lpc", "E_w", "H_w", "nu_w", "t_h_nm", "H_h", "scratch_nm"]
        inp = {}
        for k in num_keys:
            try:
                inp[k] = float(self.vars[k].get())
            except ValueError:
                raise ValueError(f"'{k}' 입력이 숫자가 아닙니다.")
        for k in ["hard_scale", "H_source", "P_unit", "particle", "film", "rough_mode"]:
            inp[k] = self.vars[k].get()
        pos = ["d_nm", "rho_p", "E_w", "H_w", "E_p", "scratch_nm", "f_wet", "f_rate", "t_h_nm", "H_h",
               "active_frac", "d_L_um"]
        pos += ["Sq_nm", "Spc", "Spd"] if inp["rough_mode"] == "AFM" else ["sigma_um", "beta_um", "eta_mm2"]
        for k in pos:
            if inp[k] <= 0:
                raise ValueError(f"'{k}' 는 0보다 커야 합니다.")
        if inp["H_source"] == H_SRC_DIRECT and inp["H_pad_direct"] <= 0:
            raise ValueError("H_eff 직접 입력값은 0보다 커야 합니다.")
        if not (0 <= inp["porosity"] < 100):
            raise ValueError("기공률은 0 이상 100 미만이어야 합니다.")
        if not (0 <= inp["wt"] < 100):
            raise ValueError("Solid contents 는 0 이상 100 미만이어야 합니다.")
        return inp

    def on_compute(self):
        self._sync_roughness()
        try:
            inp = self._read_inputs()
            r = compute(inp)
        except ValueError as e:
            messagebox.showerror("입력 오류", str(e))
            return
        except Exception as e:
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
            f"ψ = {r['psi']:.3g} [{regime_label(r['psi'])}]   P_crit = {pc_txt}   "
            f"현재 P = {inp['P_val']:g} {u} → {'임계 이상 (경도 지배)' if above else '임계 이하 (탄성·압력 지배)'}   "
            f"E_tip = {fmt_p(r['E_tip'])}   H_eff = {fmt_p(r['H_eff'])}\n"
            f"F_p = {fmt_F(r['F_p'])}   p_m = {fmt_p(r['p_m'])}   "
            f"① 나노스크래치: {'발생' if r['nano'] else '불가'} (99.9 % 입자: {'발생' if r['nano_999'] else '불가'}; d_c {fmt_len(r['d_c_nano'])})   "
            f"② 깊은 스크래치: {'발생' if r['deep'] else '불가'} (d_c {fmt_len(r['d_c_deep'])})"
        )
        self.txt.configure(state="normal")
        self.txt.delete("1.0", tk.END)
        self.txt.insert(tk.END, build_report(inp, r))
        self.txt.configure(state="disabled")
        self._draw_plot()

    # ---------- 그래프
    def _draw_plot(self):
        c = self.canvas
        c.delete("all")
        if not self.result:
            return
        W, H = c.winfo_width(), c.winfo_height()
        if W < 50 or H < 50:
            return
        r = self.result
        sw = r["sweep"]
        ml, mr, mt, mb = 70, 70, 18, 34
        x0, x1, y0, y1 = ml, W - mr, mt, H - mb
        lx_min, lx_max = -3.0, math.log10(30.0)

        def X(P):
            return x0 + (math.log10(max(P / PSI, 1e-9)) - lx_min) / (lx_max - lx_min) * (x1 - x0)

        Fs = [f for _, _, f in sw if f != float("inf")]
        Fmax = max(Fs) * 1.15 if Fs else 1.0

        def YF(F):
            return y1 - (F / Fmax) * (y1 - y0)

        def Yf(fr):
            return y1 - fr * (y1 - y0)

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

        pts_f = [(X(P), Yf(fr)) for P, fr, _ in sw]
        pts_F = [(X(P), YF(F)) for P, _, F in sw if F != float("inf")]
        if len(pts_f) > 1:
            c.create_line(*sum(pts_f, ()), fill="#1f77b4", width=2)
        if len(pts_F) > 1:
            c.create_line(*sum(pts_F, ()), fill="#d62728", width=2)
        # 수화층 관통에 필요한 하중 수준
        if 0 < r["F_need"] < Fmax:
            yy = YF(r["F_need"])
            c.create_line(x0, yy, x1, yy, fill="#ff7f0e", dash=(3, 3))
            c.create_text(x1 - 4, yy - 7, text="F_need (수화층 관통)", anchor="e", fill="#ff7f0e", font=("Consolas", 8))

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
        self.on_rough_mode()
        self.on_compute()


if __name__ == "__main__":
    App().mainloop()
