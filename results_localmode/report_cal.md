# C-H local-mode 배음 계산 (GFN2-xTB → B3LYP/6-31G** 보정, 물 용매 모델)

포도당 전체 C-H 1차 배음 평균 (계산 원값): 5750 cm⁻¹ (1739 nm)

Δ = 해당 C-H의 1차 배음 − 포도당 평균 (cm⁻¹). 절대값보다 Δ가 더 믿을 만함.

## 1. C-H별 결과

| 분자 | C-H | 유형 | 방향 | ωe | ωexe | 1차 배음 cm⁻¹ | Δ | 2차 배음 cm⁻¹ |
|---|---|---|---|---|---|---|---|---|
| a-glucose | C1-H | anomeric H-C(O5)(OH) | eq | 2917 | 68.6 | 5753 | +3 | 8484 |
| a-glucose | C2-H | H-C-O-H | ax | 2984 | 65.2 | 5857 | +107 | 8644 |
| a-glucose | C3-H | H-C-O-H | ax | 2891 | 69.7 | 5713 | -37 | 8419 |
| a-glucose | C4-H | H-C-O-H | ax | 2887 | 72.5 | 5696 | -54 | 8387 |
| a-glucose | C5-H | H-C-O-C (ring O5) | ax | 2914 | 69.6 | 5745 | -5 | 8469 |
| a-glucose | C6-Ha | CH2-OH | CH2 | 2945 | 68.4 | 5790 | +40 | 8540 |
| a-glucose | C6-Hb | CH2-OH | CH2 | 2980 | 65.8 | 5849 | +99 | 8634 |
| b-glucose | C1-H | anomeric H-C(O5)(OH) | ax | 2810 | 75.4 | 5581 | -169 | 8210 |
| b-glucose | C2-H | H-C-O-H | ax | 2904 | 69.1 | 5732 | -18 | 8449 |
| b-glucose | C3-H | H-C-O-H | ax | 2892 | 70.8 | 5710 | -40 | 8410 |
| b-glucose | C4-H | H-C-O-H | ax | 2894 | 71.3 | 5710 | -39 | 8412 |
| b-glucose | C5-H | H-C-O-C (ring O5) | ax | 2888 | 71.3 | 5703 | -47 | 8402 |
| b-glucose | C6-Ha | CH2-OH | CH2 | 2966 | 66.4 | 5828 | +78 | 8601 |
| b-glucose | C6-Hb | CH2-OH | CH2 | 2970 | 66.3 | 5832 | +82 | 8608 |
| b-maltose | non-red. C1-H | anomeric H-C(O5)(O-C(glyc)) | eq | 2941 | 71.3 | 5774 | +24 | 8515 |
| b-maltose | non-red. C2-H | H-C-O-H | ax | 2927 | 70.1 | 5759 | +9 | 8488 |
| b-maltose | non-red. C3-H | H-C-O-H | ax | 2904 | 71.2 | 5724 | -26 | 8432 |
| b-maltose | non-red. C4-H | H-C-O-H | ax | 2883 | 72.9 | 5689 | -61 | 8376 |
| b-maltose | non-red. C5-H | H-C-O-C (ring O5) | ax | 2925 | 71.8 | 5751 | +1 | 8474 |
| b-maltose | non-red. C6-Ha | CH2-OH | CH2 | 2959 | 67.1 | 5815 | +65 | 8581 |
| b-maltose | non-red. C6-Hb | CH2-OH | CH2 | 2969 | 66.8 | 5829 | +79 | 8602 |
| b-maltose | red. C1-H | anomeric H-C(O5)(OH) | ax | 2828 | 74.3 | 5610 | -140 | 8256 |
| b-maltose | red. C2-H | H-C-O-H | ax | 2903 | 69.2 | 5731 | -19 | 8447 |
| b-maltose | red. C3-H | H-C-O-H | ax | 2909 | 71.5 | 5730 | -20 | 8440 |
| b-maltose | red. C4-H | H-C-O-C(glyc) | ax | 2949 | 71.6 | 5783 | +33 | 8524 |
| b-maltose | red. C5-H | H-C-O-C (ring O5) | ax | 2904 | 72.2 | 5721 | -29 | 8427 |
| b-maltose | red. C6-Ha | CH2-OH | CH2 | 2968 | 67.3 | 5826 | +76 | 8598 |
| b-maltose | red. C6-Hb | CH2-OH | CH2 | 2988 | 66.7 | 5856 | +106 | 8644 |
| a-maltose | non-red. C1-H | anomeric H-C(O5)(O-C(glyc)) | eq | 2944 | 71.3 | 5777 | +27 | 8520 |
| a-maltose | non-red. C2-H | H-C-O-H | ax | 2881 | 73.4 | 5685 | -65 | 8369 |
| a-maltose | non-red. C3-H | H-C-O-H | ax | 2910 | 70.7 | 5734 | -16 | 8447 |
| a-maltose | non-red. C4-H | H-C-O-H | ax | 2887 | 72.4 | 5696 | -54 | 8388 |
| a-maltose | non-red. C5-H | H-C-O-C (ring O5) | ax | 2912 | 73.1 | 5728 | -22 | 8437 |
| a-maltose | non-red. C6-Ha | CH2-OH | CH2 | 2998 | 65.1 | 5875 | +125 | 8676 |
| a-maltose | non-red. C6-Hb | CH2-OH | CH2 | 2915 | 71.1 | 5740 | -10 | 8458 |
| a-maltose | red. C1-H | anomeric H-C(O5)(OH) | eq | 2920 | 68.0 | 5758 | +8 | 8494 |
| a-maltose | red. C2-H | H-C-O-H | ax | 2984 | 65.4 | 5856 | +106 | 8642 |
| a-maltose | red. C3-H | H-C-O-H | ax | 2914 | 69.8 | 5743 | -7 | 8465 |
| a-maltose | red. C4-H | H-C-O-C(glyc) | ax | 2938 | 72.3 | 5765 | +15 | 8495 |
| a-maltose | red. C5-H | H-C-O-C (ring O5) | ax | 2935 | 70.0 | 5770 | +20 | 8506 |
| a-maltose | red. C6-Ha | CH2-OH | CH2 | 2935 | 69.3 | 5774 | +24 | 8515 |
| a-maltose | red. C6-Hb | CH2-OH | CH2 | 3025 | 65.2 | 5911 | +161 | 8729 |

## 2. 유형별 평균 (모든 분자, Δ cm⁻¹, ± 표준편차, 개수)

| 유형 | 축 C-H | 적도 C-H | 전체 |
|---|---|---|---|
| CH2-OH | – | – | +77 ± 43 (n=12) |
| H-C-O-C (ring O5) | -14 ± 22 (n=6) | – | -14 ± 22 (n=6) |
| H-C-O-C(glyc) | +24 ± 9 (n=2) | – | +24 ± 9 (n=2) |
| H-C-O-H | -15 ± 50 (n=16) | – | -15 ± 50 (n=16) |
| anomeric H-C(O5)(O-C(glyc)) | – | +26 ± 2 (n=2) | +26 ± 2 (n=2) |
| anomeric H-C(O5)(OH) | -154 ± 14 (n=2) | +6 ± 3 (n=2) | -74 ± 81 (n=4) |
| (방향만) ax | | | -22 ± 57 (n=26) |
| (방향만) eq | | | +16 ± 10 (n=4) |
| (방향만) CH2 | | | +77 ± 43 (n=12) |

## 3. 포도당 vs 말토오스 (그램당, 아노머 비율 반영, C-H당 FWHM 60 cm⁻¹ 가정)

- 말토오스가 더 큰 구간 (차이 > 최대의 50 %): 5737–5784 cm⁻¹ (Δ -13~+34)
- 포도당이 더 큰 구간: 5557–5592 cm⁻¹ (Δ -193~-158), 5817–5854 cm⁻¹ (Δ +67~+104)
- 말토오스 글리코시드 C-H(비환원 C1-H, 환원 C4-H) 밴드 중심: b-maltose non-red. C1-H 5774 (Δ +24), b-maltose red. C4-H 5783 (Δ +33), a-maltose non-red. C1-H 5777 (Δ +27), a-maltose red. C4-H 5765 (Δ +15)
- 최대 상대 차이 |Δ스펙트럼|/최대흡광: 0.432

## 4. xTB → DFT 보정 (1차 배음)

DFT = 2109 + 0.672 × xTB  (쌍 8개, 잔차 RMSE 28 cm⁻¹)

| 분자 | C-H | 유형 | xTB | DFT | 보정 xTB |
|---|---|---|---|---|---|
| b-glucose | C1-H | anomeric H-C(O5)(OH) | 5171 | 5595 | 5584 |
| b-glucose | C2-H | H-C-O-H | 5484 | 5754 | 5794 |
| b-glucose | C4-H | H-C-O-H | 5437 | 5727 | 5762 |
| b-glucose | C5-H | H-C-O-C (ring O5) | 5400 | 5723 | 5737 |
| b-glucose | C6-Ha | CH2-OH | 5689 | 5917 | 5932 |
| a-glucose | C1-H | anomeric H-C(O5)(OH) | 5499 | 5836 | 5804 |
| a-glucose | C2-H | H-C-O-H | 5650 | 5932 | 5905 |
| me-a-glucoside | C1-H | anomeric H-C(O5)(O-C(glyc)) | 5447 | 5802 | 5769 |

## 5. 파장으로 환산한 결과 (보정값)

- 말토오스가 더 큰 구간: 5737–5784 cm⁻¹ = 1729–1743 nm
- 포도당이 더 큰 구간: 5557–5592 cm⁻¹ = 1788–1800 nm, 5817–5854 cm⁻¹ = 1708–1719 nm