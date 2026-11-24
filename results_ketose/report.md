# 과당 / 알룰로스 영역 선택 결과 (4000–16000 cm⁻¹, 호변이성질체 모델)

- 반복 5회, 검량 150개 / 테스트 100개 (독립)
- 시뮬레이션 비용: 스펙트럼당 10 ms
- 물 모델: real (문헌 흡수 곡선), 광로 1.0 mm
- 시료: 총 당 5–40 g/100 mL, 알룰로스 비율 0–100 %, 온도 25±2 °C (호변이성질체 비율이 온도에 따라 변함)

## 호변이성질체 비율 (25 °C, 문헌 NMR 근사)

| | β-피라노스 | α-피라노스 | β-퓨라노스 | α-퓨라노스 | 케토 |
|---|---|---|---|---|---|
| fructose | 68.3 % | 2.7 % | 22.4 % | 6.2 % | 0.5 % |
| allulose | 24.0 % | 22.0 % | 15.0 % | 38.7 % | 0.3 % |

## 1-b. 온도 이동 테스트 (검량 25 °C → 테스트 35 °C, RMSEP / R²)

| 방법 | 과당 (g/100 mL) | 알룰로스 (g/100 mL) | 총 당 (g/100 mL) | 알룰로스 비율 % (두 모델 예측값으로 계산) |
|---|---|---|---|---|
| full_4000_16000 | 7.912 / 0.038 | 7.015 / 0.356 | 2.181 / 0.948 | 116023793.537 / -42525569614057.891 |
| nir_4000_10000 | 7.935 / 0.036 | 7.037 / 0.352 | 2.184 / 0.948 | 115822920.262 / -42545210385333.906 |
| theory_tautomer | 1.320 / 0.972 | 0.926 / 0.988 | 0.368 / 0.999 | 9.145 / 0.861 |
| theory_dominant_only | 1.701 / 0.955 | 1.695 / 0.963 | 0.422 / 0.998 | 20.726 / 0.356 |
| VIP | 8.166 / 0.155 | 8.650 / 0.029 | 3.967 / 0.841 | 32312526.246 / -5657908082692.341 |
| iPLS | 1.166 / 0.979 | 1.254 / 0.977 | 0.461 / 0.998 | 9.047 / 0.882 |
| hybrid | 1.228 / 0.978 | 0.895 / 0.988 | 0.431 / 0.998 | 8.421 / 0.907 |

## 1. 성능 (25 °C 테스트 RMSEP ± SD / R² / 사용 포인트 수)

| 방법 | 과당 (g/100 mL) | 알룰로스 (g/100 mL) | 총 당 (g/100 mL) | 알룰로스 비율 % (두 모델 예측값으로 계산) |
|---|---|---|---|---|
| full_4000_16000 | 5.656±0.978 / 0.595 / 1501 | 5.187±0.847 / 0.658 / 1501 | 1.186±0.176 / 0.986 / 1501 | 41.994±8.341 / -1.149 / 1501 |
| nir_4000_10000 | 5.627±1.019 / 0.598 / 751 | 5.160±0.888 / 0.661 / 751 | 1.186±0.175 / 0.986 / 751 | 41.660±7.673 / -1.104 / 751 |
| theory_tautomer | 0.484±0.034 / 0.997 / 223 | 0.440±0.043 / 0.998 / 221 | 0.272±0.016 / 0.999 / 219 | 3.225±0.950 / 0.987 / 223 |
| theory_dominant_only | 0.480±0.049 / 0.997 / 224 | 0.461±0.108 / 0.997 / 224 | 0.271±0.009 / 0.999 / 219 | 3.116±0.423 / 0.989 / 226 |
| VIP | 7.791±0.608 / 0.256 / 115 | 7.978±0.268 / 0.197 / 127 | 1.936±0.185 / 0.963 / 117 | 30135554.687±57047407.378 / -4824737567115.205 / 127 |
| iPLS | 0.451±0.039 / 0.997 / 240 | 0.412±0.023 / 0.998 / 216 | 0.279±0.010 / 0.999 / 202 | 2.684±0.466 / 0.991 / 246 |
| hybrid | 0.458±0.021 / 0.997 / 166 | 0.432±0.039 / 0.998 / 159 | 0.272±0.018 / 0.999 / 141 | 2.942±0.684 / 0.990 / 173 |

## 2. 추천 영역 (hybrid, 반복 60 % 이상에서 선택)

### 과당 (g/100 mL)
- **4480–4512 cm⁻¹** (2216–2232 nm, 물 A 0.73–0.75): fructose b-pyr: OH(ax) str + CO str
- **5472–5768 cm⁻¹** (1734–1827 nm, 물 A 0.26–0.43): allulose a-fur: CH2 1st overtone; fructose b-pyr: CH2 1st overtone; fructose b-pyr: CH2_ring 1st overtone; fructose b-pyr: CH_ax 1st overtone
- **6096–6360 cm⁻¹** (1572–1640 nm, 물 A 0.24–0.34): allulose b-pyr: OH_hb 1st overtone (6508 cm⁻¹) 밴드의 경사면
- **6440–6632 cm⁻¹** (1508–1553 nm, 물 A 0.40–0.83): fructose b-pyr: OH_eq 1st overtone; allulose a-fur: OH_hb 1st overtone; allulose a-fur: OH_prim 1st overtone; fructose b-pyr: OH_anomK 1st overtone
- **9304–9344 cm⁻¹** (1070–1075 nm, 물 A 0.01–0.01): allulose b-pyr: OH_hb 2nd overtone (9485 cm⁻¹) 밴드의 경사면
- **9424–9664 cm⁻¹** (1035–1061 nm, 물 A 0.01–0.01): fructose b-pyr: OH_eq 2nd overtone; allulose a-fur: OH_hb 2nd overtone; allulose a-fur: OH_prim 2nd overtone; fructose b-pyr: OH_anomK 2nd overtone

### 알룰로스 (g/100 mL)
- **5504–5768 cm⁻¹** (1734–1817 nm, 물 A 0.26–0.39): allulose a-fur: CH2 1st overtone; fructose b-pyr: CH2 1st overtone; fructose b-pyr: CH2_ring 1st overtone; fructose b-pyr: CH_ax 1st overtone
- **5808–5848 cm⁻¹** (1710–1722 nm, 물 A 0.24–0.25): allulose b-fur: CH_fur 1st overtone (5705 cm⁻¹) 밴드의 경사면
- **6096–6632 cm⁻¹** (1508–1640 nm, 물 A 0.24–0.83): fructose b-pyr: OH_eq 1st overtone; allulose a-fur: OH_hb 1st overtone; allulose a-fur: OH_prim 1st overtone; fructose b-pyr: OH_anomK 1st overtone
- **9448–9656 cm⁻¹** (1036–1058 nm, 물 A 0.01–0.01): fructose b-pyr: OH_eq 2nd overtone; allulose a-fur: OH_hb 2nd overtone; allulose a-fur: OH_prim 2nd overtone; fructose b-pyr: OH_anomK 2nd overtone

### 총 당 (g/100 mL)
- **4576–4656 cm⁻¹** (2148–2185 nm, 물 A 0.74–0.92): fructose b-pyr: CH2OH OH str + bend (4740 cm⁻¹) 밴드의 경사면
- **5496–5832 cm⁻¹** (1715–1820 nm, 물 A 0.25–0.40): allulose a-fur: CH2 1st overtone; fructose b-pyr: CH2 1st overtone; fructose b-pyr: CH2_ring 1st overtone; fructose b-pyr: CH_ax 1st overtone
- **6240–6432 cm⁻¹** (1555–1603 nm, 물 A 0.29–0.39): allulose b-pyr: OH_hb 1st overtone (6508 cm⁻¹) 밴드의 경사면
- **6520–6584 cm⁻¹** (1519–1534 nm, 물 A 0.52–0.69): fructose b-pyr: OH_eq 1st overtone; allulose a-fur: OH_hb 1st overtone; fructose b-pyr: OH_anomK 1st overtone; allulose b-pyr: OH_hb 1st overtone
- **6720–6744 cm⁻¹** (1483–1488 nm, 물 A 1.03–1.08): fructose b-pyr: OH_ax 1st overtone (6662 cm⁻¹) 밴드의 경사면
- **9312–9560 cm⁻¹** (1046–1074 nm, 물 A 0.01–0.01): allulose a-fur: OH_hb 2nd overtone; allulose b-pyr: OH_hb 2nd overtone; allulose a-pyr: OH_hb 2nd overtone; fructose b-fur: OH_eq 2nd overtone
- **15776–15832 cm⁻¹** (632–634 nm, 물 A 0.00–0.00): fructose b-pyr: OH_ax 2nd overtone (9728 cm⁻¹) 밴드의 경사면

### 공통 합집합
4480–4512, 4576–4656, 5472–5848, 6096–6632, 6720–6744, 9304–9664, 15776–15832 cm⁻¹ (전체 포인트의 13 %)

## 3. 방법별 영역 (반복 0, cm⁻¹)

**과당 (g/100 mL)**
- theory_tautomer: 4480–4512, 5488–5864, 6040–6640, 8304–8392, 9000–9328, 9392–9680
- theory_dominant_only: 4520–4608, 5472–5880, 6096–6672, 8136–8184, 8264–8384, 8448–8496, 9120–9304, 9448–9736
- VIP: 4000–4272, 4840–5408, 6496–6616
- iPLS: 5504–5792, 6104–6696, 9400–9696, 14504–14792
- hybrid: 4480–4512, 5488–5784, 6040–6336, 6440–6632, 9000–9328, 9392–9680

**알룰로스 (g/100 mL)**
- theory_tautomer: 4472–4512, 5504–5864, 6040–6640, 8312–8392, 9000–9328, 9392–9680, 15936–15992
- theory_dominant_only: 4512–4624, 5480–5880, 6088–6680, 8264–8376, 8448–8496, 9104–9328, 9464–9752
- VIP: 4000–4312, 4776–5400, 6520–6632
- iPLS: 5504–5792, 6104–6696, 9400–9696
- hybrid: 5504–5696, 5808–5864, 6040–6632, 9000–9096, 9200–9296, 9496–9680

**총 당 (g/100 mL)**
- theory_tautomer: 4352–4384, 4576–4712, 5448–5832, 6216–6608, 8072–8184, 9336–9616, 9848–9992, 15736–16000
- theory_dominant_only: 4352–4408, 4544–4704, 5456–5824, 6168–6568, 8088–8208, 9248–9536, 9800–9880, 15728–16000
- VIP: 4000–4320, 4760–5408
- iPLS: 4304–4592, 5504–5792, 6400–6696, 7304–7592, 14504–14792
- hybrid: 4352–4384, 4576–4672, 5448–5832, 6216–6408, 6520–6608, 8072–8168, 9336–9616, 9848–9944, 15936–16000
