# 1. 연소가스 쪽 모델

코드: `regencool/gas.py`, `regencool/cantera_gas.py`, `regencool/contour.py`

연소가스 모델은 윤곽의 각 위치 x에서 다음 값을 1D 냉각 해석에 넘깁니다.

| 값 | 의미 |
|---|---|
| `M` | 국소 Mach 수 |
| `p`, `T` | 정압, 정온 |
| `Taw` | 단열벽 온도 (회복 온도) |
| `h0` | 경계층 물성 인자를 뺀 Bartz 계수 |
| `g1` | T0 / T (Bartz σ에 쓰임) |

## 1.1 윤곽

- 좌표: 분사기면이 x = 0, 노즐 출구 쪽으로 x가 커집니다. 길이 단위는 mm입니다.
- r(x)는 고온가스가 닿는 내벽 반경입니다. 점 사이는 PCHIP(형상 보존 3차)으로 보간합니다.
- **매개변수 윤곽** (`type = "parametric"`): 원통 → 필렛 R2 → 수렴 원뿔 → 노즐목 상류 원호 R1 →
  노즐목 하류 원호 Rn → 벨(2차 Bézier) 또는 원뿔. R1, R2, Rn은 노즐목 반경 대비 비율로 줍니다.
  벨 길이는 같은 면적비 15° 원뿔 길이의 `bell_fraction`배입니다.
- **CSV 윤곽** (`type = "csv"`): 헤더가 있는 CSV에서 `x_mm`, `r_mm` 열을 읽습니다. 측정 스캔이나
  CAD에서 뽑은 윤곽을 쓸 때는 이 방식을 씁니다. 점이 4개 이상이어야 합니다.
- **노즐목 곡률반경 rc** (Bartz에 쓰임): `rc_throat_mm`로 주거나, 매개변수 윤곽이면 ½(R1 + Rn),
  CSV면 노즐목 ±r_t/2 구간의 2차 맞춤으로 추정합니다.

## 1.2 이상기체 모델 (`model = "ideal"`)

비열비 γ와 분자량 MW가 일정한 기체입니다. 손 계산 확인이나 빠른 검토에 씁니다.

- cp = γR/(γ − 1), R = R_u/MW
- 점도 기본값 (Bartz 추정식): μ = 1.184×10⁻⁷ · MW^0.5 · T0^0.6 [Pa·s]
- Prandtl 수 기본값 (Eucken): Pr = 4γ / (9γ − 5)
- 면적비 A/A* = ε(x)에서 등엔트로피 관계를 이분법으로 풀어 M을 구합니다. 노즐목 상류는 아음속,
  하류는 초음속 해를 씁니다.
- T = T0 / g1, p = pc · g1^(−γ/(γ−1)), g1 = 1 + (γ − 1)/2 · M²
- c*는 `cstar_m_s`로 직접 주거나 이상 c*에 `eta_cstar`를 곱합니다.

## 1.3 Cantera 모델 (`model = "cantera"`)

추진제 조합이 정해졌으면 이 모델을 씁니다.

**연소실.** 압력 pc에서 단열 HP 평형을 풉니다. 공급 엔탈피는 추진제마다 다음 중 하나로 줍니다.

- 메커니즘에 있는 화학종 (`species`): 공급 온도 `T_K`에서 기체상 엔탈피 − 기화열 `dh_vap_kJ_kg`
  (기체로 공급하면 0)
- 화학식 (`formula`, C/H/O/N만): 298.15 K 생성 엔탈피 `hf_kJ_mol` + 액체 비열 `cp_J_kgK` × (T − 298.15)

극저온 추진제는 `T_K = 298.15`로 두고 `dh_vap_kJ_kg = h_gas(298.15 K) − h_liquid(T_feed)`를 넣습니다.

**노즐 팽창.** 엔트로피 s = s0를 유지하며 pc에서 압력을 낮춥니다.

- `mode = "eq"`: 매 압력에서 SP 평형 (이동 평형)
- `mode = "fr"`: 연소실 조성으로 동결

노즐목은 질량유속 G = ρu (u = √(2(h0 − h)))가 최대인 압력이고, 각 상태의 면적비는 ε = G_t / G입니다.
비열비를 따로 가정하지 않고 실제 엔탈피 차이로 속도를 구하므로, 온도에 따라 바뀌는 비열비가 자동으로
반영됩니다. 회복 온도는 Taw = T(h + r u²/2, p), r = Pr^(1/3)입니다.

**Bartz에 쓰는 물성.** 연소실의 동결 cp, 점도, Prandtl 수(Cantera 혼합 평균 수송 물성)를 씁니다.
c*는 `eta_cstar` × Cantera c*입니다. T0는 효율 100 %의 평형 온도이므로 벽 열유속에는 보수적입니다.

**재생 결합** (`regen_coupling = true`). 냉각제가 연료이면 재킷에서 받은 열만큼 연료 공급 엔탈피가
올라갑니다.

    Δh_fuel = Δh_coolant × ṁ_coolant / ṁ_fuel

1D 해석 → Δh_fuel 갱신 → 연소실 평형 재계산을 Δh_fuel 변화가 5 kJ/kg 이하가 될 때까지 (최대 3회)
반복합니다. 냉각제 일부만 연소실로 가면 `fuel_mdot_kg_s`를 실제 연료 유량으로 둡니다.

Cantera 표는 `결과폴더/gas_cache/`에 저장되어 같은 조건이면 다시 계산하지 않습니다.

## 1.4 Bartz 열전달계수

Bartz (1957):

    h0 = 0.026 / Dt^0.2 · (μ^0.2 cp / Pr^0.6) · (pc / c*)^0.8 · (Dt / rc)^0.1 · (At / A)^0.9
    σ  = 1 / [ (½ · Twg/T0 · g1 + ½)^(0.8 − 0.2ω) · g1^(0.2ω) ],   ω = 0.6
    hg = C_hg · h0 · σ(Twg)

- σ는 벽온도 Twg에 의존하므로 1D 해석의 벽온도 반복 안에서 매번 다시 계산합니다.
- `C_hg`는 보정 배율입니다. 시험 벽온도나 CFD 열유속과 맞도록 조정하세요. 보정 전에는 1.0을 쓰고
  결과의 절대값에 ±20–30 % 불확도가 있다고 보세요.
- 열유속은 q = hg (Taw − Twg)입니다.

## 1.5 노즐 박리 판정

`gas`와 `run`은 출구 압력 pe와 주변 압력 pa(`p_ambient_bar`)를 비교합니다.

- Summerfield: pe/pa < 0.4이면 박리 위험으로 표시
- Schmucker: p_sep = pa (1.88 Me − 1)^(−0.64), pe < p_sep이면 박리

판정만 하고 박리 후의 유동은 계산하지 않습니다.

## 1.6 한계

- 연소실 안의 온도·조성 분포는 추정하지 않습니다. 분사기면부터 균일한 완전 연소 가스이고, 위치에 따른
  변화는 면적 변화에 의한 등엔트로피 가속뿐입니다. 연소 효율은 c*에만 들어가고, 벽 근처 혼합비 변화,
  막냉각, 미연소 구간은 반영하지 않습니다.
- 노즐은 1D 등엔트로피입니다. 충격파, 박리, 2D 경계층 성장, 벽면 촉매 반응은 없습니다.
- 과팽창 노즐을 해면에서 시험하면 박리 이후 구간의 hg는 실제와 다릅니다.

## 참고문헌

- Bartz, D. R., "A Simple Equation for Rapid Estimation of Rocket Nozzle Convective Heat Transfer
  Coefficients," *Jet Propulsion* 27(1), 1957.
- Goodwin, D. G. et al., *Cantera: An object-oriented software toolkit for chemical kinetics,
  thermodynamics, and transport processes*, https://cantera.org
- Summerfield, M., Foster, C. R., Swan, W. C., "Flow Separation in Overexpanded Supersonic Exhaust Nozzles,"
  *Jet Propulsion* 24, 1954.
- Schmucker, R. H., "Flow Processes in Overexpanding Chemical Rocket Nozzles," NASA TM-77396, 1984.
- Huzel, D. K., Huang, D. H., *Modern Engineering for Design of Liquid-Propellant Rocket Engines*, AIAA, 1992.
