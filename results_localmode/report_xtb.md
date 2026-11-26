# C-H local-mode 배음 계산 (GFN2-xTB, 물 용매 모델)

포도당 전체 C-H 1차 배음 평균 (계산 원값): 5419 cm⁻¹ (1845 nm)

Δ = 해당 C-H의 1차 배음 − 포도당 평균 (cm⁻¹). 절대값보다 Δ가 더 믿을 만함.

## 1. C-H별 결과

| 분자 | C-H | 유형 | 방향 | ωe | ωexe | 1차 배음 cm⁻¹ | Δ | 2차 배음 cm⁻¹ |
|---|---|---|---|---|---|---|---|---|
| a-glucose | C1-H | anomeric H-C(O5)(OH) | eq | 2917 | 68.6 | 5423 | +4 | 7945 |
| a-glucose | C2-H | H-C-O-H | ax | 2984 | 65.2 | 5577 | +159 | 8179 |
| a-glucose | C3-H | H-C-O-H | ax | 2891 | 69.7 | 5364 | -55 | 7850 |
| a-glucose | C4-H | H-C-O-H | ax | 2887 | 72.5 | 5339 | -80 | 7804 |
| a-glucose | C5-H | H-C-O-C (ring O5) | ax | 2914 | 69.6 | 5411 | -7 | 7924 |
| a-glucose | C6-Ha | CH2-OH | CH2 | 2945 | 68.4 | 5479 | +60 | 8027 |
| a-glucose | C6-Hb | CH2-OH | CH2 | 2980 | 65.8 | 5566 | +147 | 8165 |
| b-glucose | C1-H | anomeric H-C(O5)(OH) | ax | 2810 | 75.4 | 5168 | -251 | 7545 |
| b-glucose | C2-H | H-C-O-H | ax | 2904 | 69.1 | 5393 | -26 | 7894 |
| b-glucose | C3-H | H-C-O-H | ax | 2892 | 70.8 | 5359 | -60 | 7837 |
| b-glucose | C4-H | H-C-O-H | ax | 2894 | 71.3 | 5360 | -59 | 7840 |
| b-glucose | C5-H | H-C-O-C (ring O5) | ax | 2888 | 71.3 | 5349 | -70 | 7826 |
| b-glucose | C6-Ha | CH2-OH | CH2 | 2966 | 66.4 | 5534 | +116 | 8117 |
| b-glucose | C6-Hb | CH2-OH | CH2 | 2970 | 66.3 | 5541 | +122 | 8127 |
| b-maltose | non-red. C1-H | anomeric H-C(O5)(O-C(glyc)) | eq | 2941 | 71.3 | 5454 | +35 | 7991 |
| b-maltose | non-red. C2-H | H-C-O-H | ax | 2927 | 70.1 | 5433 | +14 | 7952 |
| b-maltose | non-red. C3-H | H-C-O-H | ax | 2904 | 71.2 | 5380 | -39 | 7869 |
| b-maltose | non-red. C4-H | H-C-O-H | ax | 2883 | 72.9 | 5328 | -91 | 7788 |
| b-maltose | non-red. C5-H | H-C-O-C (ring O5) | ax | 2925 | 71.8 | 5420 | +2 | 7930 |
| b-maltose | non-red. C6-Ha | CH2-OH | CH2 | 2959 | 67.1 | 5516 | +97 | 8088 |
| b-maltose | non-red. C6-Hb | CH2-OH | CH2 | 2969 | 66.8 | 5536 | +117 | 8118 |
| b-maltose | red. C1-H | anomeric H-C(O5)(OH) | ax | 2828 | 74.3 | 5210 | -209 | 7612 |
| b-maltose | red. C2-H | H-C-O-H | ax | 2903 | 69.2 | 5391 | -28 | 7892 |
| b-maltose | red. C3-H | H-C-O-H | ax | 2909 | 71.5 | 5389 | -30 | 7881 |
| b-maltose | red. C4-H | H-C-O-C(glyc) | ax | 2949 | 71.6 | 5468 | +49 | 8003 |
| b-maltose | red. C5-H | H-C-O-C (ring O5) | ax | 2904 | 72.2 | 5375 | -43 | 7862 |
| b-maltose | red. C6-Ha | CH2-OH | CH2 | 2968 | 67.3 | 5532 | +114 | 8112 |
| b-maltose | red. C6-Hb | CH2-OH | CH2 | 2988 | 66.7 | 5577 | +158 | 8179 |
| a-maltose | non-red. C1-H | anomeric H-C(O5)(O-C(glyc)) | eq | 2944 | 71.3 | 5459 | +41 | 7998 |
| a-maltose | non-red. C2-H | H-C-O-H | ax | 2881 | 73.4 | 5321 | -97 | 7777 |
| a-maltose | non-red. C3-H | H-C-O-H | ax | 2910 | 70.7 | 5395 | -24 | 7892 |
| a-maltose | non-red. C4-H | H-C-O-H | ax | 2887 | 72.4 | 5339 | -80 | 7805 |
| a-maltose | non-red. C5-H | H-C-O-C (ring O5) | ax | 2912 | 73.1 | 5386 | -32 | 7877 |
| a-maltose | non-red. C6-Ha | CH2-OH | CH2 | 2998 | 65.1 | 5605 | +186 | 8226 |
| a-maltose | non-red. C6-Hb | CH2-OH | CH2 | 2915 | 71.1 | 5403 | -15 | 7907 |
| a-maltose | red. C1-H | anomeric H-C(O5)(OH) | eq | 2920 | 68.0 | 5431 | +13 | 7960 |
| a-maltose | red. C2-H | H-C-O-H | ax | 2984 | 65.4 | 5576 | +157 | 8177 |
| a-maltose | red. C3-H | H-C-O-H | ax | 2914 | 69.8 | 5408 | -10 | 7917 |
| a-maltose | red. C4-H | H-C-O-C(glyc) | ax | 2938 | 72.3 | 5442 | +23 | 7962 |
| a-maltose | red. C5-H | H-C-O-C (ring O5) | ax | 2935 | 70.0 | 5449 | +30 | 7977 |
| a-maltose | red. C6-Ha | CH2-OH | CH2 | 2935 | 69.3 | 5454 | +36 | 7990 |
| a-maltose | red. C6-Hb | CH2-OH | CH2 | 3025 | 65.2 | 5658 | +239 | 8304 |

## 2. 유형별 평균 (모든 분자, Δ cm⁻¹, ± 표준편차, 개수)

| 유형 | 축 C-H | 적도 C-H | 전체 |
|---|---|---|---|
| CH2-OH | – | – | +115 ± 65 (n=12) |
| H-C-O-C (ring O5) | -20 ± 33 (n=6) | – | -20 ± 33 (n=6) |
| H-C-O-C(glyc) | +36 ± 13 (n=2) | – | +36 ± 13 (n=2) |
| H-C-O-H | -22 ± 74 (n=16) | – | -22 ± 74 (n=16) |
| anomeric H-C(O5)(O-C(glyc)) | – | +38 ± 3 (n=2) | +38 ± 3 (n=2) |
| anomeric H-C(O5)(OH) | -230 ± 21 (n=2) | +8 ± 4 (n=2) | -111 ± 120 (n=4) |
| (방향만) ax | | | -33 ± 84 (n=26) |
| (방향만) eq | | | +23 ± 15 (n=4) |
| (방향만) CH2 | | | +115 ± 65 (n=12) |

## 3. 포도당 vs 말토오스 (그램당, 아노머 비율 반영, C-H당 FWHM 60 cm⁻¹ 가정)

- 말토오스가 더 큰 구간 (차이 > 최대의 50 %): 5397–5464 cm⁻¹ (Δ -22~+45)
- 포도당이 더 큰 구간: 5141–5184 cm⁻¹ (Δ -278~-235), 5338–5370 cm⁻¹ (Δ -81~-49), 5526–5571 cm⁻¹ (Δ +107~+152)
- 말토오스 글리코시드 C-H(비환원 C1-H, 환원 C4-H) 밴드 중심: b-maltose non-red. C1-H 5454 (Δ +35), b-maltose red. C4-H 5468 (Δ +49), a-maltose non-red. C1-H 5459 (Δ +41), a-maltose red. C4-H 5442 (Δ +23)
- 최대 상대 차이 |Δ스펙트럼|/최대흡광: 0.548