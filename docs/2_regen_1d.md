# 2. 1D 재생냉각 해석

코드: `regencool/solver1d.py`, `regencool/channel.py`, `regencool/correlations.py`, `regencool/fluids.py`

## 2.1 흐름 방향과 압력 기준

- **대향류.** 냉각제는 노즐 출구 쪽 끝(x 최대)으로 들어가 분사기 쪽 끝(x = 0)으로 나옵니다.
  해석은 출구 쪽에서 분사기 쪽으로 x 간격 `dx_mm`씩 진행합니다.
- `P_in_bar`는 입구 매니폴드의 정압입니다. 채널 입구에서 K_in · ½ρv²를 잃고 채널로 들어갑니다.
- 채널 출구에서 K_out · ½ρv²를 더 잃습니다.
- 출력하는 압력강하는 두 가지입니다.
  - `dP_ch`: 채널 안 (마찰 + 가속)
  - `dP_tot`: 입구 매니폴드 정압 − 채널 출구 정압 + 출구 손실

## 2.2 채널 형상 정의

모든 치수는 mm이고, t와 s는 **윤곽에 수직한 방향**으로 잽니다.

| 기호 | 키 | 의미 |
|---|---|---|
| t(x) | `t_mm` | 고온벽 두께: 가스 쪽 내벽에서 채널 바닥까지 |
| s(x) | `s_mm` | 채널 크기 (형상마다 아래 표의 의미) |
| n | `n` | 채널 개수 |
| β(x) | `helix_deg` | 나선각 (0 = 축 방향) |
| rib | `rib_mm` | rect 채널 사이 리브 폭 |

t, s, β는 상수 또는 `[[x_mm, 값], ...]` 표로 줍니다. 표는 PCHIP으로 보간하고, 범위 밖에서는 끝 값을 유지합니다.

| 형상 | s의 의미 | 폭 | 단면적 A | 비고 |
|---|---|---|---|---|
| `rect` | 높이 | w = 피치 − rib | w · s | 폭이 크면 분할 리브로 나눔 (아래) |
| `circ` | 지름 d | d | πd²/4 | web = 피치 − d |
| `tear` | 바닥에서 꼭짓점까지 전체 높이 | 2R | (3π/4 + 1)R² | R = s/(1 + √2), 원 위에 45° 지붕. LPBF에서 지지대 없이 출력 가능 |

- **채널 중심 반경**: r_m = r + (t + s/2) cos α, 여기서 α는 벽 기울기각 (cos α = 1/√(1 + r′²)).
- **피치** (채널 경로에 수직): 피치 = 2π r_m / N_eff · cos β
- **N_eff** = n · 360 / (360 − n_land · land_deg). 볼트나 계측 포트용 솔리드 랜드가 있으면 채널이 차지하는
  각도가 줄어 실효 피치 개수가 커집니다.
- **web 최소 두께** (`web_min_mm`, circ/tear만): 피치가 좁아지는 노즐목 근처에서 s를 자동으로 줄여
  채널 사이 벽이 web_min 이상이 되게 합니다.
  - circ: s ≤ 피치 − web_min
  - tear: s ≤ (피치 − web_min)/2 · (1 + √2)
  - 어느 경우든 s는 0.3 mm 아래로 줄이지 않습니다.
- **rect 분할**: `split_x_mm` 하류(x ≥ split_x)에서 채널을 둘로 나눕니다 (분할 리브 폭 `split_rib_mm`).
  `w_split_mm`을 주면 하위 채널 폭이 그 값 이하가 될 때까지 자동으로 더 나눕니다.
- **경로 길이와 곡률**: 나선 경로의 x당 길이 ds/dx, 법선 곡률 κ_n, 측지 곡률 κ_g를 계산합니다.
  κ_n은 열전달의 Dean 인자에, √(κ_n² + κ_g²)는 마찰 보정에 씁니다.

## 2.3 열저항 직렬 모델

각 x 위치에서 피치 하나, 축 길이 1 mm당 열유량 q′는 세 열저항의 직렬 연결로 구합니다.

    q′ = (Taw − Tb) / (R_g + R_w + R_c)
    R_g = 1 / (hg · A_hot)                       A_hot = 2πr / N_eff · √(1 + r′²) · dx
    R_w = t_eff / (k_w · A_mid)                   A_mid: 벽 두께 중앙에서 같은 식
    R_c = 1 / (hc · n_s · (P_d + η_fin · P_f) · ds)

- **P_d**는 고온벽을 직접 마주하는 젖은 둘레, **P_f**는 핀(리브 옆면) 둘레, **n_s**는 피치당 하위 채널 수입니다.
- **핀 효율**: η_fin = tanh(mH)/(mH), m = √(2hc / (k_w t_f)). H는 핀 높이, t_f는 핀 두께입니다.

  | 형상 | P_d | P_f | H | t_f | t_eff |
  |---|---|---|---|---|---|
  | rect | w | 2s | s | rib | t |
  | circ | πd/2 | πd/2 | d/2 | 피치 − πd/4 | t + 0.107d |
  | tear | πR | 둘레 − πR | √2 R | 피치 − A_윗부분/H | t + 0.107 · 2R |

  t_eff의 0.107d 항은 둥근 채널 바닥 때문에 늘어난 평균 전도 거리를 근사합니다.
- **벽온도 반복**: Twg와 Twc(냉각제 쪽 벽온도)를 고정점 반복으로 풉니다. 이완은 0.6/0.4, 수렴 기준은 0.05 K,
  최대 80회입니다. hg의 σ(Twg)와 냉각제 벽면 물성이 벽온도에 의존하기 때문에 반복이 필요합니다.
- **에너지 보존**: 냉각제는 온도가 아니라 엔탈피로 진행합니다.

      H ← H + q′ · N_eff / n · dx / ṁ_ch

  그래서 액체 → 유사임계 → 초임계 구간에서 비열이 크게 변해도 연속적으로 처리합니다.
  테스트 `test_energy_balance`가 벽 통과 열량과 냉각제 엔탈피 상승이 1 % 안에서 같은지 확인합니다.

## 2.4 냉각제 쪽 열전달과 마찰

Re, Pr은 bulk 온도의 물성으로 계산하고, 수력직경은 D_h = 4A/P입니다.

**층류 (Re < 2300)**

    Nu = max( Nu_완전발달(a),  1.86 · Gz^(1/3) · (μ_b/μ_w)^0.14 ),   Gz = Re · Pr · D_h / L

- 첫 항은 Shah & London의 사각 덕트 완전발달 값(종횡비 a)입니다.
- 둘째 항은 Sieder–Tate 열적 입구 영역 식이고, L은 냉각제 입구부터의 경로 길이입니다.
- 마찰계수도 Shah & London의 f·Re(a)를 씁니다.
- circ/tear는 원관 상수로 바꿉니다: f × 64/56.9, Nu × 4.36/3.61.

**난류 (Re > 4000)**

    Nu = Nu_Gnielinski · E · C_prop
    Nu_Gnielinski = (f_s/8)(Re − 1000)Pr / [1 + 12.7 √(f_s/8) (Pr^(2/3) − 1)],   f_s = (0.79 ln Re − 1.64)⁻²

- **거칠기 증가**: E = √(f_rough / f_s), 최대 2.
  f_rough는 거칠기 `roughness_um`의 Colebrook–White 식(Haaland 초기값)으로 구합니다.
  적층 제조 구리 채널은 10–50 µm가 보통입니다.
- **물성 보정 C_prop**
  - T_b < 0.9 T_pc (과냉 액체): (μ_b/μ_w)^0.11, 최대 1.35. μ_w는 min(T_w, 0.97 T_pc)에서 평가합니다.
  - T_b ≥ 0.9 T_pc (유사임계 근처·초임계): Jackson 형식 (ρ_w/ρ_b)^0.3 · (c̄p/cp_b)^0.4, 0.5–1.5로 제한.
    c̄p = (H(T_w) − H(T_b)) / (T_w − T_b)입니다.
  - T_pc는 물성표에서 cp가 최대인 온도(유사임계온도)입니다.
- **천이 (2300–4000)**: 층류와 난류 값을 선형으로 섞습니다.

**곡률 (Dean 형식, 난류만)**

    Dn = Re · (D_h / 2R_c)²,   Dn > 1이면 hc × Dn^(±0.05)

곡률반경 R_c = 1/|κ_n|입니다. 고온벽이 굽은 경로의 바깥쪽이면 +(증가), 안쪽이면 −(감소)를 씁니다.

- 노즐목(r(x)가 최소, r″ > 0)에서는 고온벽이 바깥쪽이라 hc가 커집니다.
- 원통에서 수렴부로 꺾이는 필렛(r″ < 0)에서는 고온벽이 안쪽이라 hc가 줄어듭니다.
- 나선 채널은 원주 방향 곡률(−cos α / r_m 성분)이 더해져 고온벽 쪽 hc를 줄이는 방향입니다.

**압력강하**

    dP_f = f · F_c · (ds · dx) / D_h · ½ρv²,     F_c = Dn^0.05 (난류, 전체 곡률 √(κ_n² + κ_g²) 사용)
    dP_a = (ṁ/A)² · (1/ρ_다음 − 1/ρ)             (밀도 감소에 의한 가속)

**rect 고온벽 판 굽힘** (참고값): σ = Δp (w/t)² / 2, Δp = 냉각제 압력 − 가스 압력.
정확한 응력은 2D 단면 해석([3장](3_thermal_stress_2d.md))으로 확인합니다.

## 2.5 냉각제 물성

`TabulatedFluid`는 압력 하나에서 만든 온도별 물성표입니다 (H, ρ, μ, k, cp). 기본 압력은 `P_in_bar`이고
`property_pressure_bar`로 바꿀 수 있습니다.

**내장 이소프로판올 (`fluid = "isopropanol"`)**

| 구간 | 모델 |
|---|---|
| T ≤ 450 K (액체) | DIPPR 형식 상관식 (ρ, μ, k, cp) |
| T > 450 K | 체적 보정 Peng–Robinson (엔탈피, 밀도), Chung 희박기체 점도 + Jossi–Stiel–Thodos 고밀도 보정, Stiel–Thodos 열전도도 |
| 450–560 K | 두 모델이 연속되도록 배율을 주고, 560 K까지 배율을 1로 서서히 줄임 |

- 정확도는 액체에서 2–5 %, 유사임계점 근처에서 15–30 % 정도입니다. 열전도도의 임계 증가는 없습니다.
- 포화압 psat은 Lee–Kesler 식이고, 포화 여유 계산에 씁니다.

**다른 냉각제 (`fluid = "csv"`)**

`T_K, H_J_kg, rho_kg_m3, mu_Pa_s, k_W_mK [, cp_J_kgK]` 열을 가진 CSV를 줍니다. 예를 들어 NIST WebBook의
등압 데이터를 그대로 쓸 수 있습니다. `Tc`, `Pc_bar`, `omega`를 주면 임계압 여유와 포화 여유도 계산합니다.

**끓음 방지 확인 (출력 지표)**

- `P_over_Pc`: 채널 최소 압력 − 임계압 [bar]. 양수여야 합니다.
- `psat_margin`: T_b < T_c인 구간에서 min(P − psat) [bar].

끓음이나 막비등 모델은 없습니다. 압력이 임계압 아래로 내려가면 결과를 믿을 수 없습니다.

## 2.6 출력

`solution_1d.csv`의 열은 다음과 같습니다.

`x_mm, r_hot_mm, t_mm, s_mm, beta_deg, pitch_mm, web_mm, n_sub, Tb_K, P_bar, p_gas_bar, Taw_K, Twg_K, Twc_K,
q_MW_m2, q_cool_MW_m2, hg_W_m2K, hc_W_m2K, Re, v_m_s, Dh_mm, eta_fin, Dean_factor, turbulent`

- `q_MW_m2`는 가스 쪽 열유속 hg(Taw − Twg), `q_cool_MW_m2`는 냉각제가 받는 열을 같은 고온벽 면적으로 나눈
  값입니다. 축방향 전도를 끄면 둘이 같고, 켜면 벽이 축방향으로 옮긴 열만큼 다릅니다.
- 축방향 전도를 켜면 `Twg_no_axial_K` (첫 행진, 축방향 전도 없음)가 붙습니다.
- 채널 압력이 임계압 아래로 내려가는 곳이 있으면 끓음 판정 열(`Tsat_K, Tsat_margin_K, q_c_MW_m2, r_onb,
  q_chf_MW_m2, r_chf, r_chf_bound, boil_flag`)이 붙고 `boiling.png`가 생깁니다 (2.9절).

`summary.json`의 `metrics`는 다음을 담습니다.

- 최대 Twg와 그 위치, 노즐목 Twg, 최대 Twc
- 냉각제 출구 온도, Δh, 총 열량 Q, 최대 열유속
- ΔP (전체/채널), 채널 출구 압력, 최대 유속, 최소 Re, 층류 비율
- 최소 web, 최대 폭, 최소 s, 최대 나선각
- 임계압·포화 여유
- 축방향 전도를 켜면: 반복 횟수, 수렴 여부, 축방향 전도 없는 최대 Twg, 가스 쪽 총 열량 `Qg_kW`
- `summary.json`의 `boiling`: 끓음 판정 요약 (2.9절), `axial_hist`: 축방향 반복 기록

## 2.7 한계

- 원주방향 불균일 가열(혼합비 줄무늬, 랜드 옆 채널)은 다루지 않습니다.
- 축방향 벽 열전도는 기본으로 꺼져 있습니다 (2.8절). 끄면 국소 1D 열저항이라, 열저항이 x 방향으로 갑자기 바뀌는
  곳(아래 물성 보정 전환점, 형상 접선점)에서 Twg가 계단처럼 변합니다.
- 물성 보정이 T_b = 0.9 T_pc에서 점도비 → Jackson으로 바뀌어 Twg가 계단처럼 변할 수 있습니다.
  원래 검증된 모델과 같게 두려고 일부러 섞지 않았습니다. 이 근처 결과는 두 식의 차이만큼 불확도가 있다고 보세요.
  축방향 전도를 켜면 계단이 벽 두께의 몇 배 길이로 퍼집니다.
- 끓음은 판정만 하고, 열전달 계산 자체는 단상입니다 (2.9절).
- 매개변수 윤곽의 원호와 직선이 만나는 곳에서는 곡률이 불연속이라 Dean 인자 때문에 hc에 작은 계단이 생깁니다.
- 채널 사이의 유량 불균일은 이 해석에 없고, [4장](4_design_studies.md)의 공차·매니폴드 연구에서 다룹니다.

## 2.8 축방향 벽 열전도 (`axial_conduction = true` 또는 `--axial`)

기본 1D 행진은 각 x 위치를 고립된 가스–벽–냉각제 열저항으로 풉니다. 그래서 국소 열저항이 x 방향으로 갑자기 바뀌면
(물성 보정 전환, 원호–직선 접선점, 채널 크기 꺾임) 그 계단이 Twg에 그대로 나타납니다. 실제 벽은 구리처럼 열전도도가
높으면 축방향으로 열을 퍼뜨려 이런 계단을 벽 두께 몇 배 길이로 완만하게 만듭니다.

**전역 벽 회로.** 행진이 끝나면 모든 위치의 (Twg_j, Twc_j)를 한 번에 푸는 선형 회로를 만듭니다.

- 각 위치의 가스 a = hg·A_hot, 벽 b = k_w·A_mid/t_eff, 냉각제 c = hc·A_wet 전도도(W/K)는 직전 행진의 값으로 고정합니다.
- 이웃한 위치 j, j+1 사이에 축방향 전도도를 둡니다.
  - **고온벽 판(slab)**: K_s = k_w · (w_mid · t) / L_mid. w_mid는 벽 두께 중앙의 피치 폭,
    L_mid는 벽 중앙면(윤곽을 법선 방향으로 t/2만큼 옮긴 곡선)의 자오선 길이입니다.
    두께 방향은 쌍선형 요소 하나로 보고, 축방향 열을 고온 쪽·냉각제 쪽 절점에 1/3, 1/6 가중치로 나눕니다.
  - **채널 사이 고체**(리브, 랜드, web): K_r = k_w · A_sol · cos β / L_path. A_sol = 피치 × 채널 높이 − 채널 단면적.
    나선 채널에서는 고체가 채널 경로를 따라 이어지므로 자오선 방향 전도도에 cos β가 곱해집니다.
    고체의 평균 온도를 Tb + η_fin(Twc − Tb)로 보고(핀 온도 분포의 평균), 그 축방향 열을 냉각제 쪽 절점에 넣습니다.
- 양 끝은 단열입니다. 냉각제 출구 위치(분사기 쪽 끝)도 포함해 끝에 인위적인 계단이 생기지 않게 합니다.

**행진과의 결합.** 전역 해와 행진 해가 같아질 때까지 반복합니다.

1. 첫 행진은 축방향 전도 없이 풉니다.
2. 전역 회로를 풀고, 행진 Twg와의 최대 차이가 `ax_tol` (0.3 K)보다 작으면 끝냅니다.
3. 아니면 이웃 위치의 온도를 전역 해로 고정하고, 각 위치 자신의 온도는 2×2 회로 안에서 음해(implicit)로 두어 다시
   행진합니다. 축방향 전도도가 국소 회로보다 훨씬 큰 곳(얇은 벽, 낮은 hc)에서도 안정적입니다.
4. 이웃 온도는 Aitken 동적 이완(Irons–Tuck, ω ∈ [0.1, 1.5])으로 갱신합니다. 벽면 물성 보정 때문에 hc가 벽온도에
   의존해 단순 대입(Picard)은 주기 2로 진동할 수 있는데, Aitken은 이를 막습니다. 데모 엔진은 6회 안에 수렴합니다.

**에너지.** 축방향 열은 위치 사이에서 주고받기만 하므로 전체 합은 0입니다. 그래서 가스 쪽 총 열량 `Qg_kW`와 냉각제가
받은 총 열량 `Q_kW`가 거의 같습니다 (출구 끝 위치의 작은 차이만 남음).

**비용과 사용처.** 1D 해석이 약 6배 느려집니다 (데모 0.1 → 0.7초). 설계 반복(최적화)에는 끄고, 최종 확인과
Twg 계단이 문제 되는 곳에서 켜는 것을 권합니다. `optimize`, `tolerance`, `manifold`에서도 켜면 그대로 적용됩니다.

**근사.** 핀 안의 열이 핀 높이를 따라 어떻게 냉각제로 빠지는지는 평균 효율 하나로 줄였습니다. 원주 방향 전도와
덮개(jacket) 쪽 전도는 없습니다. 벽 단면의 온도 분포가 중요하면 [3장](3_thermal_stress_2d.md)의 2D 단면 해석을 함께 보세요.

## 2.9 끓음 여유 판정

내장 IPA 물성표는 한 압력에서 만든 단상 물성이라 끓음을 계산하지 못합니다. 채널 압력이 냉각제 임계압 아래로 내려가는
위치가 있으면, `run`이 냉각제 쪽 벽(채널 바닥, Twc)을 그 위치 압력의 포화 상태와 비교해 판정만 합니다.

| 항목 | 식 | 불확도 |
|---|---|---|
| 포화 여유 | Tsat(P) − Twc. 양수면 벽에서 증기가 생기지 않음 | 포화선 ±1–2 K |
| 핵비등 시작 (ONB) | Davis–Anderson: ΔT_ONB = √(8σ·Tsat·v_fg·q_c / (k_l·h_fg)), r_onb = (Twc − Tsat)/ΔT_ONB | ±40 % |
| CHF (최선 추정) | Hall–Mudawar (2000) 과냉 유동비등: q_CHF = G·h_fg·0.0722·We^−0.312·(ρ_f/ρ_g)^−0.644·[1 − 0.900·(ρ_f/ρ_g)^0.724·x_e], We = G²D_h/(ρ_f σ) | ±50 % (물 자료 외삽) |
| CHF (하한) | Zuber 풀비등 × Ivey–Morris 과냉 인자: q_Z = 0.131·ρ_g^½·h_fg·[σg(ρ_f − ρ_g)]^¼·[1 + 0.1·(ρ_f/ρ_g)^0.75·c_p,l·ΔT_sub/h_fg] | ±35 % |

- q_c = hc(Twc − Tb)는 채널 바닥의 냉각제 쪽 열유속, x_e = −c_p,l(Tsat − Tb)/h_fg는 평형 건도(과냉이면 음수)입니다.
- **플래그**: 0 단상, 1 ONB 도달 (r_onb ≥ 1), 2 CHF 접근 = 벽이 Tsat 이상이면서
  (q_c ≥ 0.5 q_CHF,최선 또는 q_c ≥ q_CHF,하한). 벽이 Tsat보다 낮으면 증기가 생기지 않으므로 q_c/q_CHF가 커도
  소손 위험으로 보지 않습니다.
- **포화 물성** (`regencool/boiling.py`):
  - IPA 내장: Lee–Kesler 증기압을 뒤집은 Tsat, Watson 잠열 (정상 끓는점 기준), DIPPR형 액체 밀도·전도도·비열,
    Peng–Robinson 증기 근의 증기 밀도, 대응상태 표면장력 σ ∝ (Tc − T)^(11/9).
  - 다른 냉각제: `[coolant] saturation_table`에 압력별 포화 물성 CSV를 줍니다 ([설정 키](config_reference.md#coolant)).
  - 임계압 이상이거나 h_fg < 1 kJ/kg인 위치는 건너뜁니다.

**해석.** 정상 운전에서 채널 압력이 임계압보다 높으면 판정은 비어 있습니다. 임계압 아래로 내려가는 설계(낮은 공급
압력, 시동·정지 과도)라면 이 판정이 0이어야 단상 결과를 믿을 수 있습니다. 1이나 2가 나오면 이 위치의 Twg는
단상 해석 결과일 뿐이고, 실제로는 핵비등으로 좋아지거나 막비등·근임계 열전달 악화로 크게 나빠질 수 있습니다.
특히 임계압 바로 아래(0.9–1.0 Pc)는 h_fg와 σ가 0으로 가므로 판정 자체의 불확도가 큽니다.

## 참고문헌

- Gnielinski, V., "New equations for heat and mass transfer in turbulent pipe and channel flow,"
  *Int. Chem. Eng.* 16, 1976.
- Petukhov, B. S., "Heat transfer and friction in turbulent pipe flow with variable physical properties,"
  *Advances in Heat Transfer* 6, 1970.
- Colebrook, C. F., "Turbulent flow in pipes...," *J. Inst. Civil Eng.* 11, 1939; Haaland, S. E., *J. Fluids Eng.* 105, 1983.
- Shah, R. K., London, A. L., *Laminar Flow Forced Convection in Ducts*, Academic Press, 1978.
- Sieder, E. N., Tate, G. E., "Heat transfer and pressure drop of liquids in tubes," *Ind. Eng. Chem.* 28, 1936.
- Jackson, J. D., "Fluid flow and convective heat transfer to fluids at supercritical pressure,"
  *Nuclear Engineering and Design* 264, 2013.
- Dean, W. R., "Note on the motion of fluid in a curved pipe," *Phil. Mag.* 4, 1927.
  곡률 보정 형식은 Taylor (NASA TN D-4332, 1968)와 Niino et al. (AIAA 82-1107)에서 쓰인 형태입니다.
- Peng, D.-Y., Robinson, D. B., "A new two-constant equation of state," *Ind. Eng. Chem. Fund.* 15, 1976.
- Chung, T.-H. et al., "Generalized multiparameter correlation for nonpolar and polar fluid transport
  properties," *Ind. Eng. Chem. Res.* 27, 1988.
- Jossi, J. A., Stiel, L. I., Thodos, G., "The viscosity of pure substances in the dense gaseous and liquid
  phases," *AIChE J.* 8, 1962.
- Stiel, L. I., Thodos, G., "The thermal conductivity of nonpolar substances in the dense gaseous and liquid
  regions," *AIChE J.* 10, 1964.
- Lee, B. I., Kesler, M. G., "A generalized thermodynamic correlation based on three-parameter
  corresponding states," *AIChE J.* 21, 1975.
- Poling, B. E., Prausnitz, J. M., O'Connell, J. P., *The Properties of Gases and Liquids*, 5th ed., McGraw-Hill, 2001.
- Irons, B. M., Tuck, R. C., "A version of the Aitken accelerator for computer iteration," *Int. J. Numer. Meth.
  Eng.* 1, 1969; Küttler, U., Wall, W. A., "Fixed-point fluid-structure interaction solvers with dynamic relaxation,"
  *Comput. Mech.* 43, 2008.
- Davis, E. J., Anderson, G. H., "The incipience of nucleate boiling in forced convection flow," *AIChE J.* 12, 1966.
- Hall, D. D., Mudawar, I., "Critical heat flux (CHF) for water flow in tubes — II. Subcooled CHF correlations,"
  *Int. J. Heat Mass Transfer* 43, 2000.
- Zuber, N., "Hydrodynamic aspects of boiling heat transfer," AEC Report AECU-4439, 1959;
  Ivey, H. J., Morris, D. J., "On the relevance of the vapour–liquid exchange mechanism for sub-cooled boiling heat
  transfer at high pressure," UKAEA AEEW-R 137, 1962.
- Watson, K. M., "Thermodynamics of the liquid state," *Ind. Eng. Chem.* 35, 1943.
