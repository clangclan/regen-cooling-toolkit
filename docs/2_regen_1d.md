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
  최대 60회입니다. hg의 σ(Twg)와 냉각제 벽면 물성이 벽온도에 의존하기 때문에 반복이 필요합니다.
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
q_MW_m2, hg_W_m2K, hc_W_m2K, Re, v_m_s, Dh_mm, eta_fin, Dean_factor, turbulent`

`summary.json`의 `metrics`는 다음을 담습니다.

- 최대 Twg와 그 위치, 노즐목 Twg, 최대 Twc
- 냉각제 출구 온도, Δh, 총 열량 Q, 최대 열유속
- ΔP (전체/채널), 채널 출구 압력, 최대 유속, 최소 Re, 층류 비율
- 최소 web, 최대 폭, 최소 s, 최대 나선각
- 임계압·포화 여유

## 2.7 한계

- 축방향 벽 열전도와 원주방향 불균일 가열을 무시합니다 (국소 1D 열저항).
- 물성 보정이 T_b = 0.9 T_pc에서 점도비 → Jackson으로 바뀌어 Twg가 계단처럼 변할 수 있습니다.
  원래 검증된 모델과 같게 두려고 일부러 섞지 않았습니다. 이 근처 결과는 두 식의 차이만큼 불확도가 있다고 보세요.
- 매개변수 윤곽의 원호와 직선이 만나는 곳에서는 곡률이 불연속이라 Dean 인자 때문에 hc에 작은 계단이 생깁니다.
- 채널 사이의 유량 불균일은 이 해석에 없고, [4장](4_design_studies.md)의 공차·매니폴드 연구에서 다룹니다.

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
