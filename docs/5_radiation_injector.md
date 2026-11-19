# 5. 연소가스 복사 추정과 분사기 오리피스 크기

1D 냉각 해석에 덧붙이는 두 가지 후처리입니다. 둘 다 1차 설계 판단용이며, 결과의 불확도가 큽니다.

```bash
python -m regencool radiation examples/demo_engine.toml
```

```bash
python -m regencool injector examples/demo_engine.toml
```

## 5.1 연소가스 복사 (`radiation`)

### 무엇을 묻는가

Bartz 열유속은 대류만 다룹니다. 소형 엔진은 빔 길이가 짧아 복사가 대류의 몇 % 수준인 경우가 많지만,
H2O·CO2가 많고 연소실 압력이 높으면 무시하기 어려울 수 있습니다. 이 명령은 "복사를 더하면 Twg가 몇 K 오르는가"를
대략 알려 줍니다. 냉각 해석을 다시 풀지는 않습니다.

### 방법 (`regencool/radiation.py`)

1. **가스 방사율.** Leckner (1972) H2O, CO2 전방사율 상관식에 압력 보정과 두 기체 밴드 겹침 보정을 적용합니다
   (Modest, *Radiative Heat Transfer*의 형태).
   - 위치마다 국소 정적 온도 T_g, 정압 p (1D 등엔트로피 가스 해)를 씁니다.
   - 평균 빔 길이 L_m = 0.9 D (무한 원통).
   - 분압 p_a = X·p. 몰분율은 `[radiation]`에 주거나, Cantera 모델이면 연소실 평형 조성을 씁니다.
2. **가스 흡수율.** Hottel 규칙: a_g = ε_g(T_w, p_a L_m T_w/T_g)·(T_g/T_w)^0.5.
3. **벽 열유속.** 회색 벽 (방사율 ε_w):

   q_rad = (1 + ε_w)/2 · σ (ε_g T_g⁴ − a_g T_w⁴)

   (1 + ε_w)/2는 회색 벽과 가스 사이의 유효 방사율 근사입니다.
4. **벽온도 상승.** 1D 해에서 벽+냉각제 열저항 R_c = (Twg − Tb)/q를 꺼내, 추가 열유속이 가스 쪽(Taw − Twg가
   줄어듦)과 냉각제 쪽으로 나뉜다고 봅니다.

   ΔTwg = q_rad · R_c / (1 + h_g R_c)

### 두 가지 값: 자른 값과 외삽한 값

Leckner 상관식은 400–2500 K에서 맞춘 식입니다. 연소실은 보통 이보다 뜨겁기 때문에 두 값을 함께 냅니다.

- **자른 값** (`T_clip_K`, 기본 2500 K): 방사율을 2500 K에서 평가하고, 방출 항의 T_g⁴는 실제 온도를 씁니다.
- **외삽한 값** (`extrapolate = true`): 상관식을 그대로 실제 온도까지 씁니다. 고온에서 ε_g가 빠르게 줄어 보통 더 작습니다.

두 값의 차이가 이 추정의 불확도의 일부입니다. 둘 다 작으면(대류의 몇 %, ΔTwg 수십 K 이하) 복사를 무시해도 되고,
크면 C_hg 보정이나 별도 복사 해석(WSGG, 선형 모델)을 고려하세요.

### 출력

- 화면: 경우별 총 복사열 Q_rad, 대류 총열 Qg, 최대 q_rad, 최대 q_rad/q, 최대 ΔTwg와 그 위치, 복사를 더한 최대 Twg
- `radiation.png`: 왼쪽 대류·복사 열유속 (로그 축), 오른쪽 ΔTwg
- `radiation.json`: 위치별 ε_g, a_g, q_rad, ΔTwg와 요약값

### 한계

- 그을음(soot) 복사가 없습니다. 탄화수소 연료에서 연료 과농이면 과소평가합니다.
- 1D 가스 온도는 코어 온도입니다. 실제로 벽 근처 경계층은 더 차가워 흡수가 생깁니다 (보수적인 쪽).
- 축방향 복사 교환(연소실 → 노즐목)이 없고, 위치마다 무한 원통으로 봅니다.
- 냉각 해석과 결합하지 않습니다 (냉각제 온도 상승분은 무시). 대류의 10 %를 넘으면 결과를 다시 봐야 합니다.

## 5.2 분사기 오리피스 크기 (`injector`)

### 무엇을 묻는가

추진제별 유량, 상류 상태, 연소실 압력, 유량계수 Cd, 구멍 수를 주면 오리피스 지름을 정합니다.

A = ṁ / (C_d · G),  d = √(4A / (n π))

핵심은 질량유속 G입니다. 재생냉각을 거친 연료는 뜨겁고(때로 초임계), 자가가압 산화제는 오리피스 안에서 끓기 때문에
비압축 식(SPI)이 크게 틀릴 수 있습니다.

### 모델 (`regencool/injector.py`)

| `model` | 식 | 쓰는 곳 |
|---|---|---|
| `spi` | G = √(2ρ(P1 − P2)) | 과냉 액체 (차가운 연료, 가압 산화제) |
| `ipa_real` | G = max over P' ∈ [P2, P1] of ρ(P', s1)·√(2(h1 − h(P', s1))) | 재생냉각을 거친 뜨거운 IPA. 막힘(choking)을 잡음 |
| `n2o_dyer` | G = k/(1+k)·G_SPI + 1/(1+k)·G_HEM, k = √((P1 − P2)/(Pv − P2)) | 자가가압 아산화질소 |

**`ipa_real`: 실유체 등엔트로피 팽창.**

- 상태식은 Peng–Robinson입니다. 액체 밀도는 부피 이동(volume translation)으로 보정하는데, 이동량은 293 K
  액체 밀도에 맞춥니다.
- 상류 (T1, P1)에서 엔트로피 s1을 구하고, 하류 압력 P'를 내려가며 같은 엔트로피의 상태를 찾습니다.
- 포화 돔 안으로 들어가면 균질 평형(HEM) 혼합물로 다룹니다.
  - Tsat는 PR 퓨개시티 같음 조건으로 구합니다.
  - 건도는 x = (s1 − s_l)/(s_v − s_l)이고, h와 ρ는 질량·부피 가중합입니다.
  - 돔 밖에서는 액체 쪽이면 액체 근, 증기 쪽이면 증기 근을 씁니다. 임계압 근처(P > 0.999 Pc)에서는 Gibbs 에너지가
    최소인 안정 근을 씁니다.
- G(P')가 P' = P2보다 높은 압력에서 최대가 되면 막힌 유동입니다. 이때는 그 최대값을 쓰고 "(choked)"로 표시합니다.
- 차가운 액체에서는 SPI와 거의 같아야 합니다 (테스트로 확인).

**`n2o_dyer`: Dyer 비균질 비평형 모델.**

- 포화 물성(Pv, ρ_l, ρ_v, h_l, h_v)은 ESDU 91022 상관식입니다.
- G_HEM은 상류 포화 액체에서 포화선을 따라 등엔트로피로 끓어 내려가는 균질 평형 유동의 최대 질량유속입니다.
- k는 기포 성장 시간과 체류 시간의 비를 나타냅니다. 과냉도가 클수록 SPI에 가깝습니다.
- 상류 압력이 증기압보다 낮으면(P1 < 0.99 Pv) 상류가 이미 끓고 있으므로 오류를 냅니다.
  상류 온도나 압력을 확인하세요.

**재생 출구와 결합 (`from_regen = true`).** 연료 오리피스 상류 상태를 1D 냉각 해석의 채널 출구에서 가져옵니다.

- 온도는 냉각제 출구 온도입니다.
- 압력은 채널 출구 압력에서 출구 손실과 `dP_feed_bar`를 뺀 값입니다.

그래서 채널 설계를 바꾸면(압력강하, 출구 온도) 연료 오리피스 크기도 함께 바뀝니다. 공급 압력 예산을 정할 때
"탱크 → 채널 → 오리피스 → 연소실" 압력 사슬 전체를 이 명령으로 확인할 수 있습니다.

### 출력

- 화면: 추진제별 모델, 상류 T1·P1, 하류 P2, G, Cd, 구멍 수 × 지름, 막힘 여부
- `injector.json`: 같은 값과 모델별 부가 정보 (막힌 압력, 포화 압력, k 등)

### 불확도

- **Cd가 가장 큰 미지수입니다.** 날카로운 모서리 짧은 오리피스는 0.6–0.85이고, L/D, 입구 라운드, 캐비테이션에 따라
  달라집니다. 반드시 수류 시험으로 보정하세요. 지름은 Cd^−½에 비례합니다.
- Dyer 모델은 공학적 맞춤식으로, G에 ±20–30 % 오차가 있습니다.
- `ipa_real`의 PR 상태식은 액체 엔탈피·밀도에 몇 % 오차가 있고, 임계점 근처에서 더 큽니다. HEM은 오리피스가 짧으면
  G를 과소평가하는 경향이 있습니다 (비평형 효과).
- 분무 품질, 충돌 각, 압력강하 비(연소 안정성 기준, 보통 ΔP/Pc ≥ 15–20 %)는 다루지 않습니다.

## 참고문헌

- Leckner, B., "Spectral and total emissivity of water vapor and carbon dioxide," *Combustion and Flame* 19, 1972.
- Modest, M. F., *Radiative Heat Transfer*, 3rd ed., Academic Press, 2013 (Ch. 11, total emissivity correlations).
- Hottel, H. C., Sarofim, A. F., *Radiative Transfer*, McGraw-Hill, 1967.
- Dyer, J., Zilliac, G., Sadhwani, A., Karabeyoglu, A., Cantwell, B., "Modeling feed system flow physics for
  self-pressurizing propellants," AIAA 2007-5702.
- ESDU 91022, "Thermophysical properties of nitrous oxide," 1991.
- Peng, D.-Y., Robinson, D. B., "A new two-constant equation of state," *Ind. Eng. Chem. Fundam.* 15, 1976;
  Péneloux, A., Rauzy, E., Fréze, R., "A consistent correction for Redlich–Kwong–Soave volumes,"
  *Fluid Phase Equilib.* 8, 1982.
