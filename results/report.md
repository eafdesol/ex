# 당류 판별용 파장 영역 선택 결과

- 반복 5회, 검량 150개 / 테스트 100개 (독립 세트)
- 시뮬레이션 비용: 스펙트럼 1개당 33.5 ms (Morse 슈뢰딩거 수치해 + 밴드 합성)
- 농도: 4종 각 0–10 g/100 mL 독립 랜덤, 수용액 (물 매트릭스, 온도 ±1.5 °C, 색소 간섭, 산란/광로 변동, 이분산 노이즈)

## 1. 성능 비교 (RMSEP 평균 ± SD, g/100 mL / R² / 파장 수)

| 방법 | 포도당(단당류) | 수크로스(이당류) | 라피노스(삼당류) | 스타키오스(사당류) | 삼+사당류 합(RFO) |
|---|---|---|---|---|---|
| full_400_2500 | 0.948±0.053 / 0.892 / 1051 | 2.478±0.163 / 0.230 / 1051 | 2.520±0.150 / 0.206 / 1051 | 2.335±0.207 / 0.389 / 1051 | 1.958±0.093 / 0.772 / 1051 |
| nir_1000_2500 | 0.530±0.038 / 0.966 / 751 | 1.431±0.168 / 0.741 / 751 | 2.518±0.145 / 0.207 / 751 | 2.171±0.164 / 0.474 / 751 | 1.193±0.142 / 0.914 / 751 |
| theory_NAS | 0.626±0.054 / 0.953 / 163 | 1.302±0.057 / 0.790 / 159 | 2.374±0.164 / 0.292 / 160 | 1.794±0.102 / 0.641 / 161 | 1.087±0.154 / 0.926 / 158 |
| VIP | 1.237±0.106 / 0.816 / 362 | 2.584±0.133 / 0.171 / 290 | 2.562±0.159 / 0.178 / 251 | 2.572±0.209 / 0.264 / 350 | 2.439±0.337 / 0.644 / 353 |
| iPLS | 0.421±0.018 / 0.979 / 200 | 1.208±0.135 / 0.819 / 190 | 2.471±0.160 / 0.233 / 70 | 1.839±0.092 / 0.622 / 145 | 1.059±0.108 / 0.932 / 190 |
| hybrid | 0.643±0.066 / 0.950 / 148 | 1.347±0.049 / 0.775 / 132 | 2.386±0.135 / 0.287 / 98 | 1.824±0.067 / 0.628 / 129 | 1.074±0.154 / 0.928 / 130 |

## 2. 당류별 추천 영역 (hybrid, 반복의 60 % 이상에서 선택된 안정 영역)

### 포도당(단당류) — glucose
- **1562–1604 nm**: OH_anom 1st overtone (1528 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **1732–1808 nm**: CH2 1st overtone; CH_anom 1st overtone; CH_furan 1st overtone; CH_glyc 1st overtone; CH_ring 1st overtone
- **2048–2068 nm**: OH(ax) str + COH bend (2087 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2104–2158 nm**: OH str + COH bend; anomeric OH str + bend
- **2268–2368 nm**: CH str + CH bend; CH2 str + scissor; CH2 str + wag; OH str + CO str; anomeric CH str + bend; anomeric OH + C1-O str; furanose CH str + bend; glycosidic CH str + bend

### 수크로스(이당류) — sucrose
- **1510–1614 nm**: OH_anom 1st overtone; OH_ax 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1746–1798 nm**: CH2 1st overtone; CH_anom 1st overtone; CH_furan 1st overtone; CH_glyc 1st overtone; CH_ring 1st overtone
- **2088–2128 nm**: OH str + COH bend; OH(ax) str + COH bend; anomeric OH str + bend
- **2154–2208 nm**: OH(ax) str + CO str (2244 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2340–2352 nm**: glycosidic CH str + bend (2324 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호

### 라피노스(삼당류) — raffinose
- **1558–1588 nm**: OH_anom 1st overtone (1528 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **1746–1758 nm**: CH_anom 1st overtone; CH_furan 1st overtone; CH_glyc 1st overtone; CH_ring 1st overtone
- **1778–1796 nm**: CH2 1st overtone; CH_ring 1st overtone
- **2066–2124 nm**: OH str + COH bend; OH(ax) str + COH bend; anomeric OH str + bend
- **2282–2286 nm**: OH str + CO str
- **2344–2350 nm**: glycosidic CH str + bend (2324 nm) 밴드의 장파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호

### 스타키오스(사당류) — stachyose
- **1510–1608 nm**: OH_anom 1st overtone; OH_ax 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1740–1754 nm**: CH_anom 1st overtone; CH_furan 1st overtone; CH_glyc 1st overtone
- **1786–1798 nm**: CH2 1st overtone
- **2066–2124 nm**: OH str + COH bend; OH(ax) str + COH bend; anomeric OH str + bend
- **2168–2202 nm**: OH(ax) str + CO str (2244 nm) 밴드의 단파장 쪽 경사면 — 물 밴드와 겹치지 않는 미분 신호
- **2280–2286 nm**: OH str + CO str

### 삼+사당류 합(RFO) — RFO_sum
- **1522–1550 nm**: OH_anom 1st overtone; OH_prim 1st overtone; OH_sec 1st overtone
- **1750–1788 nm**: CH2 1st overtone; CH_anom 1st overtone; CH_furan 1st overtone; CH_glyc 1st overtone; CH_ring 1st overtone
- **2066–2126 nm**: OH str + COH bend; OH(ax) str + COH bend; anomeric OH str + bend
- **2268–2382 nm**: CH str + CH bend; CH2 str + scissor; CH2 str + wag; OH str + CO str; anomeric CH str + bend; anomeric OH + C1-O str; furanose CH str + bend; glycosidic CH str + bend

### 공통 모델(PLS2 등)용 합집합
1510–1614, 1732–1808, 2048–2208, 2268–2382  (총 231개 파장, 전체의 22 %)

## 3. 방법별 선택 영역 (반복 0)

**포도당(단당류)**
- theory_NAS: 1546–1606, 1744–1816, 2048–2084, 2124–2176, 2268–2360
- VIP: 400–676, 1460–1532, 1744–1756, 1782–1800, 1848–2062, 2270–2306, 2328–2352, 2462–2500
- iPLS: 1500–1598, 1750–1798, 2150–2398
- hybrid: 1546–1604, 1744–1816, 2048–2076, 2124–2152, 2268–2356

**수크로스(이당류)**
- theory_NAS: 1514–1540, 1560–1614, 1712–1738, 1784–1796, 2058–2122, 2154–2208, 2268–2308, 2328–2352
- VIP: 400–640, 1376–1428, 1468–1528, 1786–1800, 1828–1918, 1942–2094, 2280–2296, 2328–2348, 2460–2500
- iPLS: 1500–1598, 1750–1798, 2100–2248, 2300–2398
- hybrid: 1514–1540, 1560–1614, 1712–1738, 1784–1796, 2088–2116, 2154–2208, 2268–2308, 2328–2352

**라피노스(삼당류)**
- theory_NAS: 1512–1542, 1560–1614, 1714–1736, 1784–1800, 2058–2122, 2156–2210, 2270–2310, 2330–2350
- VIP: 1366–1432, 1466–1534, 1788–1920, 1946–2082, 2450–2500
- iPLS: 2150–2198, 2300–2348
- hybrid: 1560–1588, 2058–2086, 2330–2350

**스타키오스(사당류)**
- theory_NAS: 1512–1542, 1558–1616, 1716–1734, 1786–1800, 2058–2124, 2156–2212, 2272–2312, 2332–2350
- VIP: 400–650, 1386–1428, 1468–1534, 1830–2078, 2458–2500
- iPLS: 1000–1048, 2050–2098, 2300–2398
- hybrid: 1512–1540, 1588–1616, 1786–1800, 2058–2124, 2186–2212, 2272–2312, 2332–2350

**삼+사당류 합(RFO)**
- theory_NAS: 1518–1548, 1704–1788, 2052–2074, 2114–2162, 2248–2370
- VIP: 400–614, 1376–1434, 1460–1560, 1826–2076, 2456–2500
- iPLS: 1050–1098, 1250–1298, 1550–1598, 2150–2248, 2300–2448
- hybrid: 1518–1546, 1764–1788, 2052–2074, 2114–2142, 2248–2366

## 그림
- `spectra_selectivity_regions.png`
- `rmsep_comparison.png`
- `selection_frequency.png`