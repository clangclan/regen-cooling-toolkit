# regen-cooling-toolkit

소형 액체로켓엔진의 재생냉각 채널을 설계할 때 쓰는 파이썬 툴킷입니다.
TOML 설정 파일 하나에 엔진 윤곽, 연소가스, 냉각제, 채널을 적고, 명령 한 줄로 다음을 계산합니다.

- 연소가스 쪽 열전달 (Bartz, 이상기체 또는 Cantera 평형)
- 1D 대향류 냉각 해석 (벽온도, 냉각제 온도, 압력강하), 선택적으로 축방향 벽 열전도 결합
- 냉각제 끓음 여유 판정 (포화 여유, 핵비등 시작, CHF), 채널 압력이 임계압 아래로 내려갈 때
- 2D 채널 단면 열-구조 해석 (온도장, von Mises 응력)
- 채널 크기 분포 s(x) 최적화
- 제작 공차 민감도
- 입구 매니폴드 유량 분배와 채널별 공차 Monte Carlo
- 연소가스 복사 추정 (H2O/CO2)
- 분사기 오리피스 크기 (비압축, 실유체 IPA 등엔트로피, 아산화질소 Dyer), 재생 출구 상태와 결합

> **이 저장소에는 방법론과 코드만 있습니다.** 실제 엔진의 CAD 모델, 윤곽 스캔, 채널 치수, 최적화 결과,
> 매니폴드 형상, CFD 자료, 해석 결과는 들어 있지 않습니다. `examples/`의 데모 엔진은 사용법을 보여 주기
> 위한 가상의 수치입니다.

## 설치

Python 3.11 이상이 필요합니다 (설정 파일을 표준 라이브러리 `tomllib`로 읽습니다).

```bash
git clone https://github.com/clangclan/regen-cooling-toolkit.git
```

```bash
cd regen-cooling-toolkit
```

```bash
python -m venv .venv
```

가상환경을 켭니다. Windows는 `.venv\Scripts\activate`, macOS/Linux는 `source .venv/bin/activate`입니다.

```bash
pip install -e ".[cantera,test]"
```

Cantera는 `gas.model = "cantera"`를 쓸 때만 필요합니다. 빼고 설치하려면 `pip install -e .`를 쓰세요.

## 빠른 시작

```bash
python -m regencool run examples/demo_engine.toml
```

결과는 기본적으로 `results/<case name>/`에 저장됩니다. `-o 폴더`로 바꿀 수 있습니다.
설치가 제대로 됐는지 확인하려면 데모 결과를 아래 참고값과 비교하세요 (가상 엔진 기준).

| 명령 | 데모 참고값 |
|---|---|
| `run examples/demo_engine.toml` | 최대 Twg 570 °C (x = 29.2 mm), 노즐목 517 °C, ΔP 4.95 bar, 채널 출구 55.27 bar, 냉각제 출구 222 °C |
| `run` 2D 노즐목 단면 | von Mises 최대 351 MPa |
| `run examples/demo_engine.toml --axial` | 6회 반복 수렴, 최대 Twg 554 °C (x = 9.5 mm, 축방향 전도 없으면 570 °C) |
| `optimize examples/demo_engine.toml` | 최대 Twg 512 °C, ΔP 7.96 bar (한계 8 bar) |
| `tolerance examples/demo_engine.toml` | 채널 하나만 −0.1 mm → 유량 ×0.72, Twg 614 °C |
| `manifold examples/demo_engine.toml` | 채널별 σ 0.05 mm → 최대 Twg 99백분위 744 °C |
| `run examples/demo_engine_cantera.toml` | 재생 결합 후 T0 3000 K, 최대 Twg 661 °C, ΔP 9.88 bar |
| `radiation examples/demo_engine.toml` | 2500 K에서 자른 값: Q_rad 8.47 kW (대류 115.7 kW), 최대 q_rad/q 0.137, 최대 ΔTwg 44 K |
| `injector examples/demo_engine.toml` | 재생 출구 495.2 K / 55.05 bar → 연료 8 × 0.891 mm (막힘), 산화제 8 × 1.924 mm |

## 명령

```
python -m regencool <명령> case.toml [-o OUT] [--C-hg X] [--no-stress] [--axial]
```

| 명령 | 하는 일 | 주요 출력 |
|---|---|---|
| `gas` | 윤곽을 따라 연소가스 상태(Mach, p, Taw, Bartz h)를 계산하고 노즐 박리를 판정 | `gas.png`, `gas_summary.json` |
| `run` | 1D 냉각 해석 + `[stress]`에 적은 위치의 2D 단면 열응력 + 끓음 판정 (임계압 아래일 때) | `solution_1d.csv`, `solution_1d.png`, `geometry.png`, `station_x*.png`, `stress_summary.png`, `summary.json`, `boiling.png` |
| `optimize` | ΔP·최소압력 제약 아래 최대 벽온도가 가장 낮은 s(x)를 찾음 | `optimized_s.toml`, `optimization.png`, `optimization.json` |
| `tolerance` | 전체 채널 균일 오프셋 ±d, 채널 하나만 −d (같은 ΔP) | `tolerance.csv`, `tolerance.png`, `tolerance.json` |
| `manifold` | 링 헤더 정압 분포 + 채널 응답표 + 채널별 크기 산포 Monte Carlo | `manifold.png`, `manifold.json` |
| `radiation` | 1D 해 위에 연소가스 복사 열유속과 벽온도 상승 추정 | `radiation.png`, `radiation.json` |
| `injector` | `[injector]`의 추진제별 오리피스 지름 (재생 출구 상태 사용 가능) | `injector.json` |

- `--C-hg X`: Bartz 열전달계수 배율을 덮어씁니다 (`[case] C_hg`보다 우선).
- `--no-stress`: `run`에서 2D 단면 해석을 건너뜁니다 (빠른 반복용).
- `--axial`: 축방향 벽 열전도를 결합합니다 (`[case] axial_conduction`보다 우선, 약 6배 느림).

## 새 엔진에 적용하는 순서

1. **설정 파일 복사.** `examples/demo_engine.toml`을 복사해 이름을 바꿉니다. 모든 키는
   [docs/config_reference.md](docs/config_reference.md)에 있습니다.
2. **윤곽.** 노즐 내벽(고온가스 쪽) 반경 r(x)를 넣습니다. 분사기면이 x = 0입니다.
   매개변수 윤곽(원통 + 수렴 원뿔 + 노즐목 원호 + 벨)을 쓰거나 `x_mm, r_mm` CSV를 줍니다.
3. **연소가스 확인.** `python -m regencool gas case.toml`로 T0, c*, 출구 압력, 박리 판정을 봅니다.
   추진제가 정해졌으면 Cantera 모델을 쓰고, 냉각제가 연료라면 `regen_coupling = true`로 둡니다.
4. **냉각제와 채널 초안.** 유량, 입구 압력·온도, 채널 형상(rect/circ/tear), 개수, 벽 두께 t(x),
   채널 크기 s(x)를 넣고 `run --no-stress`로 돌립니다.
5. **1D 결과 판단.** 최대 Twg가 재료 한계보다 충분히 낮은지, ΔP가 공급 압력 예산 안인지,
   채널 출구 압력이 냉각제 임계압보다 높은지(끓음 방지), 층류 구간이 없는지 봅니다.
   임계압 아래로 내려가면 `run`이 끓음 판정을 출력합니다. 플래그가 0이 아니면 단상 결과를 믿을 수 없습니다.
6. **s(x) 최적화.** `optimize`가 출력한 `s_mm = [...]`를 `[channels]`에 붙여 넣습니다.
   제작 공차가 크면 `robust_d_mm`을 공차 크기로 둡니다.
7. **공차와 매니폴드.** `tolerance`로 채널 하나가 작게 나왔을 때의 고온점을, `manifold`로 헤더 압력
   분포와 채널별 산포의 영향을 확인합니다.
8. **2D 열응력.** `[stress] x_mm`에 노즐목, 최대 Twg 위치 등을 넣고 `run`으로 단면 응력을 봅니다.
   최종 확인에는 `--axial`로 축방향 전도를 켜서 Twg 계단이 완만해진 결과도 봅니다.
9. **복사와 분사기.** `radiation`으로 복사가 무시할 만한지 확인하고, `injector`로 재생 출구 상태에서 연료·산화제
   오리피스 지름과 압력 예산을 확인합니다.
10. **보정.** 연소시험 벽온도나 CFD 열유속이 생기면 `C_hg`로 Bartz를 보정하고 4–9를 다시 돌립니다.

## 문서

| 문서 | 내용 |
|---|---|
| [docs/1_hot_gas.md](docs/1_hot_gas.md) | 연소가스 모델: 이상기체, Cantera 평형/동결, Bartz, 회복온도, 노즐 박리 판정 |
| [docs/2_regen_1d.md](docs/2_regen_1d.md) | 1D 대향류 해석: 채널 형상 정의, 열저항, 냉각제 상관식, 압력강하, 냉각제 물성, 축방향 벽 열전도, 끓음 여유 판정 |
| [docs/3_thermal_stress_2d.md](docs/3_thermal_stress_2d.md) | 2D 채널 단면: 격자, 경계조건, 일반화 평면변형률, 결과 해석 |
| [docs/4_design_studies.md](docs/4_design_studies.md) | 최적화, 공차 민감도, 매니폴드 네트워크와 Monte Carlo |
| [docs/5_radiation_injector.md](docs/5_radiation_injector.md) | 연소가스 복사 추정, 분사기 오리피스 모델 (SPI, 실유체 HEM, Dyer) |
| [docs/6_transient_roadmap.md](docs/6_transient_roadmap.md) | 시동·정지 과도 해석의 모델 구성, 검증 절차, 빠진 물리, 막냉각 사전 검토 (방법론만, 코드 없음) |
| [docs/config_reference.md](docs/config_reference.md) | 설정 파일의 모든 키 |

## 코드 구성

| 모듈 | 역할 |
|---|---|
| `regencool/case.py` | TOML을 읽어 윤곽·가스·냉각제·채널 객체를 만들고 해석을 묶음 |
| `regencool/contour.py` | 내벽 윤곽 (매개변수 / CSV), 면적비, 노즐목 곡률반경 |
| `regencool/gas.py`, `cantera_gas.py` | 연소가스 상태와 Bartz 계수 |
| `regencool/fluids.py` | 냉각제 물성표 (내장 IPA, 또는 CSV) |
| `regencool/channel.py` | 채널 형상, 나선각, 피치, 곡률, 단면 성질 |
| `regencool/correlations.py` | 마찰계수, Nusselt 수, Dean 곡률 인자 |
| `regencool/solver1d.py` | 1D 대향류 행진 해석, 축방향 벽 열전도 결합 |
| `regencool/boiling.py` | 포화 물성 (내장 IPA 또는 CSV), 끓음 여유 판정 |
| `regencool/radiation.py` | Leckner H2O/CO2 방사율, 복사 열유속 추정 |
| `regencool/injector.py` | 오리피스 질량유속 모델과 지름 |
| `regencool/fe2d.py`, `station2d.py` | 2D 열전도·탄성 유한요소, 1D 결과로 경계조건 구성 |
| `regencool/optimize.py`, `tolerance.py`, `manifold.py` | 설계 연구 |
| `regencool/plots.py` | matplotlib 그림 |

파이썬에서 직접 쓸 수도 있습니다.

```python
from regencool.case import Case

case = Case.load("examples/demo_engine.toml")
res = case.couple()            # Cantera 재생 결합이 꺼져 있으면 solve()와 같음
print(case.metrics(res))       # 최대 Twg, ΔP, 출구 압력 등
```

## 모델의 한계

결과를 쓰기 전에 꼭 읽어 주세요. 자세한 내용은 각 문서의 "한계" 절에 있습니다.

- **Bartz 상관식**은 보통 ±20–30 % 오차가 있습니다. 시험이나 CFD 자료로 `C_hg`를 보정하기 전에는
  절대값보다 설계안 간 비교에 쓰세요.
- **연소실 가스**는 분사기면부터 완전 연소한 균일 평형 가스로 봅니다. 연소 효율은 c*에만 `eta_cstar`로
  들어가고, 벽 근처 혼합비 변화와 막냉각은 반영하지 않습니다.
- **노즐 팽창**은 1D 등엔트로피입니다. 충격파와 유동 박리는 계산하지 않고, 출구 압력으로 Summerfield와
  Schmucker 기준의 박리 여부만 판정합니다. 해면 시험에서 과팽창 노즐이면 박리 위치 이후 결과는 쓰지 마세요.
- **냉각제 물성표**는 압력 하나에서 만든 단상 물성입니다. 끓음은 판정만 하고 열전달로 계산하지 않으므로
  채널 출구 압력이 냉각제 임계압보다 높아야 합니다 (`P_min_bar` 기본값 = 임계압 + 1 bar). 내장 IPA 물성은 액체에서 2–5 %,
  유사임계점 근처에서 15–30 % 정도 오차가 있습니다.
- **냉각제 물성 보정의 전환.** 냉각제 온도가 유사임계온도의 0.9배를 넘으면 벽면 물성 보정이 점도비에서
  Jackson 식으로 바뀌고, 그 지점에서 Twg가 계단처럼 변할 수 있습니다.
- **윤곽 원호 경계**에서는 곡률이 불연속이라 Dean 곡률 인자 때문에 hc에 작은 계단이 생깁니다.
- **1D 해석**은 국소 열저항 직렬 모델입니다. 축방향 벽 열전도는 선택 사항이고(기본 꺼짐), 채널 사이 고체는
  평균 핀 효율로 단순화합니다. 원주방향 불균일 가열은 무시합니다.
- **복사**는 그을음이 없는 H2O/CO2 회색 가스 추정이고, 냉각 해석과 결합하지 않습니다.
- **분사기**는 1차 크기 결정입니다. Cd는 수류 시험으로 보정해야 합니다.
- **과도(시동·정지) 해석**은 아직 없습니다. [docs/6_transient_roadmap.md](docs/6_transient_roadmap.md)를 보세요.
- **2D 단면**은 x를 고정한 횡단면입니다. 벽 기울기 |dr/dx| > 0.25 (수렴부·확대부 일부)에서는 근사이며
  경고를 냅니다. 재료는 선형 탄성이라 소성, 크리프, 저주기 피로는 별도로 평가해야 합니다.
- **매니폴드**는 1D 헤더 모델(Bernoulli + 마찰 + 분기 손실/운동량 회복)입니다. 3D 유동은 CFD로 확인하세요.

## 테스트

```bash
pytest
```

윤곽 생성, 면적비-Mach 역산, 1D 에너지 보존, 자유 열팽창 무응력 조건, CLI 실행, Cantera 연소실 상태,
축방향 전도 결합의 에너지 보존과 수렴, IPA 포화 물성, 저압 끓음 판정, 오리피스 모델(차가운 액체에서 SPI와 일치),
복사 추정을 확인합니다.

## 라이선스

[MIT License](LICENSE)
