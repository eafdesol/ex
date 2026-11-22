# 당류 판별용 파장 영역 선택 결과

- 반복 5회, 검량 150개 / 테스트 100개 (독립 세트)
- 시뮬레이션 비용: 스펙트럼 1개당 44.4 ms (Morse 슈뢰딩거 수치해 + 밴드 합성)
- 농도: 4종 각 0–10 g/100 mL 독립 랜덤, 수용액 (물 매트릭스, 온도 ±1.5 °C, 색소 간섭, 산란/광로 변동, 이분산 노이즈)

## 1. 성능 비교 (RMSEP 평균 ± SD / R² / 파장 수; 단위는 타깃별)

| 방법 | 포도당(DP1) | 말토오스(DP2) | 말토트리오스(DP3) | 말토테트라오스(DP4) | DP3+DP4 합 | 총 당량 (g/100 mL) | 환원 말단 (mmol/100 mL) | 평균 중합도 DP |
|---|---|---|---|---|---|---|---|---|
| full_400_2500 | 1.258±0.064 / 0.810 / 1051 | 2.567±0.120 / 0.181 / 1051 | 2.520±0.150 / 0.206 / 1051 | 2.412±0.176 / 0.351 / 1051 | 2.158±0.143 / 0.723 / 1051 | 0.187±0.012 / 0.999 / 1051 | 2.445±0.267 / 0.985 / 1051 | 0.180±0.010 / 0.746 / 1051 |
| nir_1000_2500 | 1.144±0.022 / 0.843 / 751 | 2.557±0.104 / 0.187 / 751 | 2.517±0.145 / 0.207 / 751 | 2.304±0.165 / 0.408 / 751 | 2.091±0.117 / 0.740 / 751 | 0.159±0.009 / 0.999 / 751 | 1.311±0.212 / 0.995 / 751 | 0.177±0.010 / 0.754 / 751 |
| theory_NAS | 0.997±0.052 / 0.881 / 155 | 2.503±0.067 / 0.224 / 156 | 2.423±0.188 / 0.264 / 158 | 2.214±0.121 / 0.455 / 158 | 1.931±0.107 / 0.779 / 160 | 0.134±0.014 / 0.999 / 159 | 1.870±0.127 / 0.991 / 158 | 0.150±0.021 / 0.825 / 158 |
| VIP | 1.424±0.133 / 0.756 / 370 | 2.575±0.119 / 0.176 / 247 | 2.555±0.154 / 0.182 / 248 | 2.623±0.166 / 0.233 / 348 | 2.501±0.170 / 0.627 / 350 | 0.235±0.016 / 0.998 / 279 | 4.454±0.752 / 0.948 / 328 | 0.231±0.061 / 0.568 / 262 |
| iPLS | 1.021±0.049 / 0.875 / 175 | 2.620±0.097 / 0.147 / 50 | 2.455±0.134 / 0.245 / 70 | 2.306±0.152 / 0.407 / 60 | 1.928±0.254 / 0.775 / 90 | 0.138±0.014 / 0.999 / 185 | 1.255±0.155 / 0.996 / 200 | 0.143±0.011 / 0.840 / 140 |
| hybrid | 1.026±0.037 / 0.874 / 102 | 2.538±0.093 / 0.201 / 58 | 2.383±0.164 / 0.287 / 91 | 2.269±0.144 / 0.426 / 95 | 1.896±0.104 / 0.786 / 100 | 0.139±0.012 / 0.999 / 130 | 1.915±0.213 / 0.991 / 117 | 0.149±0.013 / 0.826 / 110 |

## 2. 당류별 추천 영역 (hybrid, 반복의 60 % 이상에서 선택된 안정 영역)

### 포도당(DP1) — glucose
- **1528–1528 nm**: OH_anom 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1546–1582 nm**: OH_anom 1st overtone
- **2068–2092 nm**: anomeric OH str + bend
- **2206–2266 nm**: OH str + CO str; anomeric OH + C1-O str
- **2282–2312 nm**: CH str + CH bend; CH2 str + scissor; OH str + CO str; anomeric CH str + bend; glycosidic CH str + bend
- **2372–2372 nm**: CH2 str + wag

### 말토오스(DP2) — maltose
- **2102–2110 nm**: OH str + COH bend; anomeric OH str + bend
- **2154–2182 nm**: OH str + COH bend (2115 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2308–2308 nm**: CH str + CH bend; CH2 str + scissor; anomeric CH str + bend

### 말토트리오스(DP3) — maltotriose
- **1552–1580 nm**: OH_anom 1st overtone (1533 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **1788–1794 nm**: CH2 1st overtone
- **2062–2082 nm**: anomeric OH str + bend (2104 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2158–2184 nm**: OH str + COH bend (2115 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2256–2260 nm**: OH str + CO str; anomeric OH + C1-O str
- **2282–2310 nm**: CH str + CH bend; CH2 str + scissor; OH str + CO str; anomeric CH str + bend; glycosidic CH str + bend
- **2338–2366 nm**: C4-H(linked) str + bend; glycosidic CH str + bend

### 말토테트라오스(DP4) — maltotetraose
- **1526–1532 nm**: OH_anom 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1552–1576 nm**: OH_anom 1st overtone (1533 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2068–2110 nm**: OH str + COH bend; anomeric OH str + bend
- **2164–2172 nm**: OH str + COH bend (2115 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2238–2240 nm**: anomeric OH + C1-O str (2261 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2282–2310 nm**: CH str + CH bend; CH2 str + scissor; OH str + CO str; anomeric CH str + bend; glycosidic CH str + bend

### DP3+DP4 합 — DP3_4_sum
- **1558–1588 nm**: OH_anom 1st overtone (1533 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **1710–1734 nm**: CH_anom 1st overtone; CH_glyc 1st overtone
- **1766–1780 nm**: CH2 1st overtone; CH_link 1st overtone; CH_ring 1st overtone
- **2116–2138 nm**: OH str + COH bend; anomeric OH str + bend
- **2178–2188 nm**: OH str + COH bend (2115 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2216–2236 nm**: anomeric OH + C1-O str (2261 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2290–2312 nm**: CH str + CH bend; CH2 str + scissor; anomeric CH str + bend; glycosidic CH str + bend

### 총 당량 (g/100 mL) — total
- **1536–1582 nm**: OH_anom 1st overtone; OH_sec 1st overtone
- **1728–1760 nm**: CH_anom 1st overtone; CH_glyc 1st overtone; CH_link 1st overtone; CH_ring 1st overtone
- **1780–1804 nm**: CH2 1st overtone; CH_ring 1st overtone
- **2048–2086 nm**: anomeric OH str + bend (2104 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2130–2158 nm**: OH str + COH bend (2115 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2274–2306 nm**: CH str + CH bend; CH2 str + scissor; OH str + CO str; anomeric CH str + bend; anomeric OH + C1-O str
- **2322–2374 nm**: C4-H(linked) str + bend; CH str + CH bend; CH2 str + scissor; CH2 str + wag; anomeric CH str + bend; glycosidic CH str + bend

### 환원 말단 (mmol/100 mL) — reducing_ends
- **1520–1544 nm**: OH_anom 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1710–1738 nm**: CH_anom 1st overtone; CH_glyc 1st overtone; CH_link 1st overtone
- **1762–1796 nm**: CH2 1st overtone; CH_link 1st overtone; CH_ring 1st overtone
- **2092–2120 nm**: OH str + COH bend; anomeric OH str + bend
- **2188–2222 nm**: anomeric OH + C1-O str (2261 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2256–2344 nm**: C4-H(linked) str + bend; CH str + CH bend; CH2 str + scissor; OH str + CO str; anomeric CH str + bend; anomeric OH + C1-O str; glycosidic CH str + bend

### 평균 중합도 DP — mean_DP
- **1520–1546 nm**: OH_anom 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1718–1738 nm**: CH_anom 1st overtone; CH_glyc 1st overtone; CH_link 1st overtone
- **1762–1772 nm**: CH2 1st overtone; CH_link 1st overtone; CH_ring 1st overtone
- **2092–2118 nm**: OH str + COH bend; anomeric OH str + bend
- **2216–2222 nm**: anomeric OH + C1-O str (2261 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2254–2372 nm**: C4-H(linked) str + bend; CH str + CH bend; CH2 str + scissor; CH2 str + wag; OH str + CO str; anomeric CH str + bend; anomeric OH + C1-O str; glycosidic CH str + bend

### 공통 모델(PLS2 등)용 합집합
1520–1588, 1710–1804, 2048–2374  (총 247개 파장, 전체의 24 %)

## 3. 방법별 선택 영역 (반복 0)

**포도당(DP1)**
- theory_NAS: 1518–1582, 2044–2092, 2106–2178, 2192–2268, 2288–2304
- VIP: 400–674, 1454–1538, 1784–1794, 1844–2078, 2274–2310, 2330–2362, 2452–2500
- iPLS: 800–848, 1500–1548, 1650–1698, 2200–2348
- hybrid: 1518–1582, 2074–2092, 2166–2178, 2192–2268, 2288–2304

**말토오스(DP2)**
- theory_NAS: 1044–1054, 1516–1574, 2046–2094, 2108–2178, 2194–2266, 2284–2308
- VIP: 1368–1432, 1464–1532, 1786–1920, 1946–2088, 2332–2354, 2456–2500
- iPLS: 1700–1748
- hybrid: 1044–1054, 2108–2136, 2284–2308

**말토트리오스(DP3)**
- theory_NAS: 1042–1054, 1514–1564, 1732–1748, 2046–2096, 2112–2180, 2196–2266, 2284–2310
- VIP: 1366–1432, 1466–1534, 1816–1920, 1946–2082, 2338–2348, 2454–2500
- iPLS: 2150–2198, 2300–2348
- hybrid: 1544–1564, 2046–2074, 2196–2224, 2256–2266

**말토테트라오스(DP4)**
- theory_NAS: 1040–1054, 1514–1556, 1728–1746, 1772–1786, 2048–2100, 2116–2182, 2200–2264, 2282–2312
- VIP: 400–650, 1382–1426, 1468–1538, 1828–2076, 2458–2500
- iPLS: 2100–2148, 2300–2348
- hybrid: 1040–1054, 1514–1556, 1728–1746, 1772–1786, 2048–2076, 2282–2310

**DP3+DP4 합**
- theory_NAS: 1506–1536, 1552–1612, 1710–1740, 1768–1794, 2082–2150, 2172–2244, 2292–2308
- VIP: 402–618, 1376–1434, 1460–1556, 1826–2074, 2456–2500
- iPLS: 1750–1798, 2200–2348
- hybrid: 1506–1534, 1552–1580, 1710–1738, 1768–1794, 2172–2244, 2292–2308

**총 당량 (g/100 mL)**
- theory_NAS: 1540–1584, 1718–1752, 1766–1798, 2052–2086, 2128–2186, 2256–2296, 2314–2374
- VIP: 426–486, 522–584, 1368–1430, 1466–1530, 1788–1802, 1820–1918, 1942–2086, 2332–2354, 2454–2500
- iPLS: 900–948, 1000–1048, 1500–1548, 1650–1698, 2150–2198, 2250–2348
- hybrid: 1540–1584, 1718–1746, 1766–1794, 2052–2086, 2256–2296, 2314–2372

**환원 말단 (mmol/100 mL)**
- theory_NAS: 1528–1562, 1704–1746, 1762–1802, 2098–2158, 2252–2320, 2340–2386
- VIP: 402–632, 1374–1428, 1464–1534, 1826–2082, 2282–2300, 2334–2354, 2458–2500
- iPLS: 1550–1598, 1750–1798, 2050–2198, 2250–2298, 2350–2448
- hybrid: 1558–1562, 1704–1746, 1762–1802, 2098–2156, 2252–2280, 2312–2320, 2340–2386

**평균 중합도 DP**
- theory_NAS: 1528–1562, 1704–1746, 1762–1802, 2098–2158, 2252–2320, 2340–2386
- VIP: 400–670, 1466–1480, 1524–1536, 1838–1992, 2278–2306, 2334–2358, 2458–2500
- iPLS: 1250–1298, 2150–2248, 2300–2398
- hybrid: 1528–1562, 1762–1790, 2098–2156, 2252–2310, 2340–2386

## 그림
- `spectra_selectivity_regions.png`
- `rmsep_comparison.png`
- `selection_frequency.png`