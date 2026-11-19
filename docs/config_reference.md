# 설정 파일 키 목록

설정은 TOML 파일 하나입니다. 예시는 `examples/demo_engine.toml` (이상기체)과
`examples/demo_engine_cantera.toml` (Cantera)에 있습니다.

- 단위는 키 이름 끝에 붙어 있습니다 (`_mm`, `_bar`, `_K`, `_C`, `_kg_s` 등).
- "필수"가 아닌 키는 생략하면 기본값을 씁니다.
- 파일 경로(`contour.file`, `coolant.table`)는 TOML 파일이 있는 폴더 기준입니다.
- **축 방향 분포**: 표에 "분포"라고 적힌 키는 상수 또는 `[[x_mm, 값], [x_mm, 값], ...]`로 줄 수 있습니다.
  점 사이는 PCHIP으로 보간하고, 범위 밖에서는 끝 값을 유지합니다.

## `[case]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `name` | `"case"` | 케이스 이름. 기본 결과 폴더 `results/<name>/` |
| `dx_mm` | 0.25 | 1D 해석의 축 방향 간격 |
| `C_hg` | 1.0 | Bartz hg 배율 (시험·CFD로 보정). `--C-hg`가 우선 |
| `p_ambient_bar` | 1.01325 | 노즐 박리 판정용 주변 압력 |
| `axial_conduction` | false | 축방향 벽 열전도 결합 ([2.8절](2_regen_1d.md#28-축방향-벽-열전도-axial_conduction--true-또는---axial)). `--axial`이 우선 |

## `[contour]`

분사기면이 x = 0인 고온가스 쪽 내벽 반경 r(x).

| 키 | 기본값 | 의미 |
|---|---|---|
| `type` | `"parametric"` | `"parametric"` 또는 `"csv"` |
| `rc_throat_mm` | 자동 | Bartz용 노즐목 곡률반경. 매개변수 윤곽은 ½(R1 + Rn), CSV는 노즐목 근처 2차 맞춤 |

`type = "csv"`:

| 키 | 기본값 | 의미 |
|---|---|---|
| `file` | 필수 | 헤더가 있는 CSV |
| `x_col`, `r_col` | `"x_mm"`, `"r_mm"` | 열 이름 |

`type = "parametric"`:

| 키 | 기본값 | 의미 |
|---|---|---|
| `r_chamber_mm` | 필수 | 연소실 반경 |
| `r_throat_mm` | 필수 | 노즐목 반경 |
| `l_cyl_mm` | 필수 | 원통부 길이 (분사기면부터 필렛 시작까지) |
| `conv_half_angle_deg` | 30 | 수렴부 반각 |
| `r1_ratio` | 1.5 | 노즐목 상류 원호 반경 / r_throat |
| `r2_ratio` | 1.0 | 원통–수렴부 필렛 반경 / r_throat |
| `rn_ratio` | 0.382 | 노즐목 하류 원호 반경 / r_throat |
| `eps_exit` | 4.0 | 출구 면적비 |
| `divergent` | `"bell"` | `"bell"` (2차 Bézier) 또는 `"cone"` |
| `theta_n_deg` | 25 | 벨 시작각. `cone`이면 원뿔 반각 |
| `theta_e_deg` | 10 | 벨 출구각 |
| `bell_fraction` | 0.8 | 벨 길이 / 같은 면적비 15° 원뿔 길이 |

## `[gas]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `model` | `"ideal"` | `"ideal"` 또는 `"cantera"` |
| `pc_bar` | 필수 | 연소실 압력 |
| `eta_cstar` | 1.0 | Bartz에 쓰는 c* = eta × 이론 c* |
| `cstar_m_s` | 없음 | c*를 직접 줄 때 (eta보다 우선) |

### `[gas.ideal]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `T0_K` | 필수 | 연소 온도 |
| `gamma` | 필수 | 비열비 (일정) |
| `MW` | 필수 | 분자량 [kg/kmol] |
| `Pr` | 4γ/(9γ − 5) | Prandtl 수 |
| `mu_Pa_s` | 1.184e-7 · MW^0.5 · T0^0.6 | 점도 |

### `[gas.cantera]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `of` | 필수 | 산화제/연료 질량비 |
| `mode` | `"eq"` | `"eq"` 이동 평형, `"fr"` 동결 |
| `mechanism` | `"gri30_highT.yaml"` | Cantera 메커니즘 (C, H, O, N) |
| `npts` | 600 | 팽창 표의 압력 점 수 |
| `regen_coupling` | false | 냉각제 엔탈피 상승을 연료 공급 엔탈피에 더해 반복 |
| `fuel_mdot_kg_s` | 냉각제 유량 | 연소실로 들어가는 연료 유량 (재생 결합용) |

### `[gas.cantera.fuel]`, `[gas.cantera.oxidizer]`

둘 중 한 형식으로 줍니다.

| 형식 | 키 | 의미 |
|---|---|---|
| 메커니즘 화학종 | `species` | 화학종 이름 (예: `"N2O"`, `"O2"`, `"CH4"`) |
| | `T_K` (298.15) | 공급 온도 |
| | `dh_vap_kJ_kg` (0) | 기화열. 기체상 엔탈피에서 뺌. 기체 공급이면 0 |
| 화학식 | `formula` | 예: `"C3H8O"`, `"C2H6O"` (C, H, O, N만) |
| | `hf_kJ_mol` | 298.15 K 액체 생성 엔탈피 |
| | `cp_J_kgK` (0) | 액체 비열 (298.15 K와 T_K 사이 보정) |
| | `T_K` (298.15) | 공급 온도 |

극저온 공급: `T_K = 298.15`, `dh_vap_kJ_kg = h_gas(298.15 K) − h_liquid(T_feed)`.

## `[coolant]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `fluid` | `"isopropanol"` | `"isopropanol"` (`"ipa"`, `"2-propanol"`) 또는 `"csv"` |
| `mdot_kg_s` | 필수 | 전체 냉각제 유량 |
| `T_in_K` | 293.15 | 입구 온도 |
| `P_in_bar` | 필수 | 입구 매니폴드 정압 (노즐 출구 쪽 끝) |
| `roughness_um` | 30 | 채널 벽 거칠기 |
| `K_in`, `K_out` | 0.3, 1.5 | 채널 입구·출구 국부 손실계수 (채널 동압 기준) |
| `property_pressure_bar` | `P_in_bar` | 물성표를 만드는 압력 |
| `saturation_table` | 없음 (IPA는 내장) | 끓음 판정용 포화 물성 CSV. 열: `P_bar, Tsat_K, hfg_J_kg, rho_l, rho_v, k_l, cp_l, sigma_N_m` |

`fluid = "csv"`:

| 키 | 기본값 | 의미 |
|---|---|---|
| `table` | 필수 | `T_K, H_J_kg, rho_kg_m3, mu_Pa_s, k_W_mK [, cp_J_kgK]` 열을 가진 CSV |
| `name` | `"csv"` | 표시 이름 |
| `Tc`, `Pc_bar`, `omega` | 없음 | 임계온도 [K], 임계압, 이심인자. 있으면 임계압·포화 여유를 계산 |

## `[channels]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `shape` | `"tear"` | `"rect"`, `"circ"`, `"tear"` |
| `n` | 필수 | 채널 개수 |
| `t_mm` | 필수, 분포 | 고온벽 두께 (윤곽에 수직) |
| `s_mm` | 필수, 분포 | 채널 크기: rect 높이, circ 지름, tear 전체 높이 |
| `rib_mm` | 1.0 | rect 리브 폭 |
| `web_min_mm` | 0 | circ/tear: 채널 사이 최소 벽 두께. 넘지 않도록 s를 자동으로 줄임 |
| `helix_deg` | 0, 분포 | 나선각 (0 = 축 방향) |
| `land_deg`, `n_land` | 0, 0 | 채널이 없는 솔리드 랜드의 각도와 개수 |
| `split_x_mm` | 없음 | rect: 이 x부터 채널을 둘로 나눔 |
| `split_rib_mm` | `rib_mm` | rect: 분할 리브 폭 |
| `w_split_mm` | 없음 | rect: 하위 채널 폭이 이 값 이하가 될 때까지 자동 분할 |

형상 정의는 [2장](2_regen_1d.md#22-채널-형상-정의)에 있습니다.

## `[material]`

고온벽과 덮개 재료 (온도 무관 상수).

| 키 | 기본값 | 의미 |
|---|---|---|
| `E_GPa` | 120 | 탄성계수 |
| `nu` | 0.33 | 푸아송비 |
| `alpha_per_K` | 17e-6 | 열팽창계수 |
| `k_W_mK` | 320 | 열전도도 (1D·2D 모두) |
| `T_ref_K` | 293.15 | 무응력 기준 온도 |
| `name` | `"wall material"` | 표시 이름 |

## `[limits]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `dP_max_bar` | 없음 | 허용 압력강하 (최적화 제약, 없으면 최적화에서 10 bar) |
| `P_min_bar` | 냉각제 임계압 + 1 bar | 채널 최소 압력 하한 (최적화 제약) |
| `Twg_max_C` | 없음 | 그림에 기준선으로만 그림 |

## `[stress]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `x_mm` | `["throat"]` | 2D 단면 위치 목록. 숫자, `"throat"`, `"max_Twg"` |
| `jacket_mm` | 2.0 | 채널 위 덮개 두께 |
| `h_mm` | 0.015 | 채널 근처 요소 크기 |

## `[optimize]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `knots_x_mm` | 6개 자동 | s(x) 매듭의 x 위치 |
| `dP_max_bar` | `[limits]` 값 | 압력강하 제약 |
| `P_min_bar` | `[limits]` 값 | 출구 압력 제약 |
| `robust_d_mm` | 0 | > 0이면 공차 ±d를 함께 고려 |
| `s_min_mm` | 0.3 | 최소 채널 크기 |
| `maxfev` | 200 | 최대 1D 해석 횟수 (강건이면 한 번 평가에 3회) |

## `[tolerance]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `offsets_mm` | [-0.2, -0.1, -0.05, 0.05, 0.1, 0.2] | 모든 채널 균일 오프셋 |
| `one_channel_mm` | [0.05, 0.1, 0.2] | 채널 하나만 −d (같은 ΔP) |

## `[manifold]`

| 키 | 기본값 | 의미 |
|---|---|---|
| `theta_deg` | 필수 | 헤더 단면 표의 각도 (0 = 입구, 180 = 반대편) |
| `area_mm2` | 필수 | 각도별 헤더 단면적 (한 갈래) |
| `dh_mm` | 필수 | 각도별 헤더 수력직경 |
| `R_mean_mm` | 필수 | 헤더 중심선 반경 |
| `k_r` | 0.8 | 분기 뒤 운동량 회복 계수 |
| `C_b` | 0.5 | 분기 입구 손실 계수 |
| `land0_deg` | 0 | 첫 랜드 위치 (채널 각도 배치용) |
| `d_grid_mm` | [-0.2, 0.2, 0.05] | 응답표의 채널 오프셋 [시작, 끝, 간격] |
| `mf_grid` | [0.4, 1.4, 0.1] | 응답표의 유량 배율 [시작, 끝, 간격] |
| `sigma_mm` | [0.025, 0.05, 0.1] | Monte Carlo 채널 크기 표준편차 목록 |
| `mean_mm` | [0.0] | Monte Carlo 평균 오프셋 목록 |
| `samples` | 2000 | 표본 수 |
| `seed` | 1 | 난수 시드 |

`d_grid_mm`는 Monte Carlo의 ±3σ보다 넓게 잡으세요. 표 밖의 값은 표 끝에서 잘립니다.

## `[radiation]`

`python -m regencool radiation`이 씁니다 ([5장](5_radiation_injector.md)).

| 키 | 기본값 | 의미 |
|---|---|---|
| `X_H2O`, `X_CO2` | Cantera 연소실 조성 | H2O, CO2 몰분율. 이상기체 모델이면 필수 |
| `eps_wall` | 0.8 | 벽 방사율 |
| `T_clip_K` | 2500 | Leckner 상관식을 평가할 최고 온도 (적합 범위 상한) |
| `extrapolate` | true | 자르지 않고 외삽한 값도 함께 보고 |

## `[injector]`

`python -m regencool injector`가 씁니다 ([5장](5_radiation_injector.md)). 하위 표 하나가 추진제 하나입니다
(이름은 자유: `[injector.fuel]`, `[injector.oxidizer]` 등).

| 키 | 기본값 | 의미 |
|---|---|---|
| `P_chamber_bar` | `gas.pc_bar` | 하류(연소실) 압력 |

`[injector.<이름>]`:

| 키 | 기본값 | 의미 |
|---|---|---|
| `model` | `"spi"` | `"spi"` 비압축 단상, `"ipa_real"` IPA 실유체 등엔트로피(HEM), `"n2o_dyer"` 아산화질소 Dyer |
| `mdot_kg_s` | 필수 | 이 추진제의 전체 유량 |
| `T_K`, `P_up_bar` | 필수 (`from_regen`이면 생략) | 오리피스 상류 온도, 압력 |
| `from_regen` | false | 상류 상태를 1D 냉각 해석의 채널 출구(출구 손실 뒤)에서 가져옴 |
| `dP_feed_bar` | 0 | `from_regen`일 때 채널 출구와 오리피스 사이 손실 |
| `rho_kg_m3` | 냉각제 물성표 | `spi`의 밀도 |
| `Cd` | 0.75 | 유량계수 |
| `n_holes` | 1 | 같은 크기 구멍 수 |
