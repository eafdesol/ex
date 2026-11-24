# 과당 / 알룰로스 영역 선택 결과 (4000–16000 cm⁻¹, 호변이성질체 모델)

- 반복 5회, 검량 150개 / 테스트 100개 (독립)
- 시뮬레이션 비용: 스펙트럼당 9 ms
- 시료: 총 당 5–40 g/100 mL, 알룰로스 비율 0–100 %, 온도 25±2 °C (호변이성질체 비율이 온도에 따라 변함)

## 호변이성질체 비율 (25 °C, 문헌 NMR 근사)

| | β-피라노스 | α-피라노스 | β-퓨라노스 | α-퓨라노스 | 케토 |
|---|---|---|---|---|---|
| fructose | 68.3 % | 2.7 % | 22.4 % | 6.2 % | 0.5 % |
| allulose | 24.0 % | 22.0 % | 15.0 % | 38.7 % | 0.3 % |

## 1-b. 온도 이동 테스트 (검량 25 °C → 테스트 35 °C, RMSEP / R²)

| 방법 | 과당 (g/100 mL) | 알룰로스 (g/100 mL) | 총 당 (g/100 mL) | 알룰로스 비율 % (두 모델 예측값으로 계산) |
|---|---|---|---|---|
| full_4000_16000 | 0.760 / 0.992 | 0.770 / 0.992 | 0.274 / 0.999 | 7.266 / 0.935 |
| nir_4000_10000 | 0.669 / 0.994 | 0.734 / 0.993 | 0.272 / 0.999 | 6.611 / 0.946 |
| theory_tautomer | 0.703 / 0.993 | 0.440 / 0.998 | 0.283 / 0.999 | 3.228 / 0.985 |
| theory_dominant_only | 0.601 / 0.995 | 0.655 / 0.995 | 0.278 / 0.999 | 4.368 / 0.975 |
| VIP | 0.649 / 0.994 | 0.718 / 0.993 | 0.438 / 0.998 | 6.503 / 0.947 |
| iPLS | 0.477 / 0.997 | 0.709 / 0.992 | 0.254 / 0.999 | 3.368 / 0.986 |
| hybrid | 0.709 / 0.993 | 0.591 / 0.995 | 0.259 / 0.999 | 4.137 / 0.976 |

## 1. 성능 (25 °C 테스트 RMSEP ± SD / R² / 사용 포인트 수)

| 방법 | 과당 (g/100 mL) | 알룰로스 (g/100 mL) | 총 당 (g/100 mL) | 알룰로스 비율 % (두 모델 예측값으로 계산) |
|---|---|---|---|---|
| full_4000_16000 | 0.356±0.013 / 0.998 / 1501 | 0.338±0.025 / 0.999 / 1501 | 0.231±0.012 / 0.999 / 1501 | 2.497±0.354 / 0.993 / 1501 |
| nir_4000_10000 | 0.355±0.021 / 0.998 / 751 | 0.341±0.023 / 0.999 / 751 | 0.230±0.014 / 0.999 / 751 | 2.569±0.434 / 0.992 / 751 |
| theory_tautomer | 0.309±0.026 / 0.999 / 229 | 0.293±0.024 / 0.999 / 230 | 0.226±0.013 / 1.000 / 219 | 1.910±0.438 / 0.996 / 230 |
| theory_dominant_only | 0.314±0.052 / 0.999 / 226 | 0.285±0.038 / 0.999 / 230 | 0.222±0.012 / 1.000 / 221 | 1.865±0.392 / 0.996 / 230 |
| VIP | 0.362±0.021 / 0.998 / 268 | 0.337±0.032 / 0.999 / 299 | 0.242±0.014 / 0.999 / 218 | 2.538±0.396 / 0.992 / 299 |
| iPLS | 0.305±0.025 / 0.999 / 231 | 0.300±0.025 / 0.999 / 233 | 0.229±0.014 / 0.999 / 119 | 1.880±0.244 / 0.996 / 263 |
| hybrid | 0.325±0.037 / 0.999 / 184 | 0.286±0.032 / 0.999 / 188 | 0.225±0.013 / 1.000 / 125 | 1.797±0.288 / 0.996 / 201 |

## 2. 추천 영역 (hybrid, 반복 60 % 이상에서 선택)

### 과당 (g/100 mL)
- **4192–4288 cm⁻¹** (2332–2385 nm): fructose b-pyr: axial CH str + bend; allulose a-fur: furanose CH str + bend; allulose a-fur: CH2OH str + wag; fructose b-pyr: CH2OH str + wag
- **4416–4552 cm⁻¹** (2197–2264 nm): fructose b-pyr: OH str + CO str; allulose a-fur: H-bonded OH str + CO str; fructose b-pyr: hemiketal OH + C2-O str; allulose a-fur: CH2OH OH str + CO str
- **4656–4848 cm⁻¹** (2063–2148 nm): fructose b-pyr: OH str + COH bend; allulose a-fur: intramol. H-bonded OH str + bend; allulose a-fur: CH2OH OH str + bend; fructose b-pyr: hemiketal OH str + bend
- **4896–4936 cm⁻¹** (2026–2042 nm): fructose b-pyr: OH(ax) str + COH bend (4799 cm⁻¹) 밴드의 경사면
- **5496–5840 cm⁻¹** (1712–1820 nm): allulose a-fur: CH2 1st overtone; fructose b-pyr: CH2 1st overtone; fructose b-pyr: CH2_ring 1st overtone; fructose b-pyr: CH_ax 1st overtone
- **6160–6624 cm⁻¹** (1510–1623 nm): fructose b-pyr: OH_eq 1st overtone; allulose a-fur: OH_hb 1st overtone; allulose a-fur: OH_prim 1st overtone; fructose b-pyr: OH_anomK 1st overtone
- **9448–9616 cm⁻¹** (1040–1058 nm): fructose b-pyr: OH_eq 2nd overtone; allulose a-fur: OH_hb 2nd overtone; fructose b-pyr: OH_anomK 2nd overtone; allulose b-pyr: OH_hb 2nd overtone

### 알룰로스 (g/100 mL)
- **4216–4312 cm⁻¹** (2319–2372 nm): allulose a-fur: CH2OH str + scissor; fructose b-pyr: CH2OH str + scissor; fructose b-pyr: ring O-CH2 str + scissor; fructose b-pyr: axial CH str + bend
- **4352–4816 cm⁻¹** (2076–2298 nm): fructose b-pyr: OH str + COH bend; allulose a-fur: intramol. H-bonded OH str + bend; allulose a-fur: CH2OH OH str + bend; fructose b-pyr: hemiketal OH str + bend
- **5512–5832 cm⁻¹** (1715–1814 nm): allulose a-fur: CH2 1st overtone; fructose b-pyr: CH2 1st overtone; fructose b-pyr: CH2_ring 1st overtone; fructose b-pyr: CH_ax 1st overtone
- **6112–6624 cm⁻¹** (1510–1636 nm): fructose b-pyr: OH_eq 1st overtone; allulose a-fur: OH_hb 1st overtone; allulose a-fur: OH_prim 1st overtone; fructose b-pyr: OH_anomK 1st overtone
- **9456–9624 cm⁻¹** (1039–1058 nm): fructose b-pyr: OH_eq 2nd overtone; allulose a-fur: OH_hb 2nd overtone; allulose a-fur: OH_prim 2nd overtone; fructose b-pyr: OH_anomK 2nd overtone

### 총 당 (g/100 mL)
- **4112–4216 cm⁻¹** (2372–2432 nm): allulose a-fur: CH2OH str + wag; fructose b-pyr: CH2OH str + wag; fructose b-pyr: ring O-CH2 str + wag; fructose b-fur: CH2OH str + wag
- **4384–4656 cm⁻¹** (2148–2281 nm): fructose b-pyr: OH str + CO str; allulose a-fur: H-bonded OH str + CO str; fructose b-pyr: hemiketal OH + C2-O str; allulose a-fur: CH2OH OH str + CO str
- **4808–4952 cm⁻¹** (2019–2080 nm): fructose b-pyr: OH(ax) str + COH bend
- **5504–5832 cm⁻¹** (1715–1817 nm): allulose a-fur: CH2 1st overtone; fructose b-pyr: CH2 1st overtone; fructose b-pyr: CH2_ring 1st overtone; fructose b-pyr: CH_ax 1st overtone
- **6424–6512 cm⁻¹** (1536–1557 nm): allulose a-fur: OH_hb 1st overtone; allulose b-pyr: OH_hb 1st overtone; allulose a-pyr: OH_hb 1st overtone; allulose b-fur: OH_hb 1st overtone
- **6696–6744 cm⁻¹** (1483–1493 nm): fructose b-pyr: OH_ax 1st overtone (6662 cm⁻¹) 밴드의 경사면

### 공통 합집합
4112–4312, 4352–4952, 5496–5840, 6112–6624, 6696–6744, 9448–9624 cm⁻¹ (전체 포인트의 16 %)

## 3. 방법별 영역 (반복 0, cm⁻¹)

**과당 (g/100 mL)**
- theory_tautomer: 4240–4376, 4432–4536, 4672–4848, 4896–4936, 5504–5864, 6064–6632, 9120–9272, 9424–9632
- theory_dominant_only: 4080–4848, 5488–5872, 6112–6672, 8296–8344, 9552–9608
- VIP: 4000–4168, 4256–4304, 4528–4608, 4672–5576, 5624–5704, 6232–6392, 6456–7080
- iPLS: 4600–4896, 5504–5792, 6104–6392, 9104–9392, 10904–11192
- hybrid: 4344–4376, 4432–4528, 4672–4848, 4896–4936, 5504–5864, 6168–6560, 9120–9216, 9424–9616

**알룰로스 (g/100 mL)**
- theory_tautomer: 4216–4360, 4432–4536, 4656–4840, 4880–4952, 5520–5864, 6064–6624, 9128–9272, 9432–9624
- theory_dominant_only: 4080–4872, 5496–5872, 6104–6680, 8296–8344, 9576–9616
- VIP: 4000–4320, 4368–4432, 4536–5704, 6240–6840, 6904–7136
- iPLS: 4600–4896, 5504–5792, 6104–6696, 8200–8496, 8800–9096, 10904–11192
- hybrid: 4216–4312, 4432–4528, 4656–4752, 5520–5864, 6064–6624, 9128–9272, 9432–9624

**총 당 (g/100 mL)**
- theory_tautomer: 4112–4472, 4560–4968, 5472–5832, 6216–6632, 15768–16000
- theory_dominant_only: 4104–4448, 4520–4960, 5480–5832, 6176–6584, 15776–16000
- VIP: 4000–4304, 4592–4704, 4776–5600, 6568–6840, 7016–7256
- iPLS: 4600–4896, 5504–5792
- hybrid: 4112–4408, 4560–4656, 4864–4952, 5576–5832, 6216–6312, 6416–6512, 15968–16000
