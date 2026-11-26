# C-H local-mode 배음 계산 (GFN2-xTB → B3LYP/6-31G** 보정, 물 용매 모델)

포도당 전체 C-H 1차 배음 평균 (계산 원값): 5752 cm⁻¹ (1738 nm)

Δ = 해당 C-H의 1차 배음 − 포도당 평균 (cm⁻¹). 절대값보다 Δ가 더 믿을 만함.

## 1. C-H별 결과

| 분자 | C-H | 유형 | 방향 | ωe | ωexe | 1차 배음 cm⁻¹ | Δ | 2차 배음 cm⁻¹ |
|---|---|---|---|---|---|---|---|---|
| a-glucose | C1-H | anomeric H-C(O5)(OH) | eq | 2917 | 68.6 | 5755 | +3 | 8488 |
| a-glucose | C2-H | H-C-O-H | ax | 2984 | 65.2 | 5859 | +107 | 8648 |
| a-glucose | C3-H | H-C-O-H | ax | 2891 | 69.7 | 5715 | -37 | 8423 |
| a-glucose | C4-H | H-C-O-H | ax | 2887 | 72.5 | 5698 | -54 | 8391 |
| a-glucose | C5-H | H-C-O-C (ring O5) | ax | 2914 | 69.6 | 5747 | -5 | 8473 |
| a-glucose | C6-Ha | CH2-OH | CH2 | 2945 | 68.4 | 5793 | +40 | 8544 |
| a-glucose | C6-Hb | CH2-OH | CH2 | 2980 | 65.8 | 5851 | +99 | 8639 |
| b-glucose | C1-H | anomeric H-C(O5)(OH) | ax | 2810 | 75.4 | 5583 | -169 | 8214 |
| b-glucose | C2-H | H-C-O-H | ax | 2904 | 69.1 | 5735 | -18 | 8453 |
| b-glucose | C3-H | H-C-O-H | ax | 2892 | 70.8 | 5712 | -40 | 8414 |
| b-glucose | C4-H | H-C-O-H | ax | 2894 | 71.3 | 5713 | -40 | 8416 |
| b-glucose | C5-H | H-C-O-C (ring O5) | ax | 2888 | 71.3 | 5705 | -47 | 8406 |
| b-glucose | C6-Ha | CH2-OH | CH2 | 2966 | 66.4 | 5830 | +78 | 8605 |
| b-glucose | C6-Hb | CH2-OH | CH2 | 2970 | 66.3 | 5835 | +83 | 8613 |
| b-maltose | non-red. C1-H | anomeric H-C(O5)(O-C(glyc)) | eq | 2941 | 71.3 | 5776 | +24 | 8520 |
| b-maltose | non-red. C2-H | H-C-O-H | ax | 2927 | 70.1 | 5762 | +9 | 8493 |
| b-maltose | non-red. C3-H | H-C-O-H | ax | 2904 | 71.2 | 5726 | -26 | 8436 |
| b-maltose | non-red. C4-H | H-C-O-H | ax | 2883 | 72.9 | 5691 | -61 | 8380 |
| b-maltose | non-red. C5-H | H-C-O-C (ring O5) | ax | 2925 | 71.8 | 5753 | +1 | 8478 |
| b-maltose | non-red. C6-Ha | CH2-OH | CH2 | 2959 | 67.1 | 5818 | +66 | 8586 |
| b-maltose | non-red. C6-Hb | CH2-OH | CH2 | 2969 | 66.8 | 5831 | +79 | 8607 |
| b-maltose | red. C1-H | anomeric H-C(O5)(OH) | ax | 2828 | 74.3 | 5612 | -141 | 8260 |
| b-maltose | red. C2-H | H-C-O-H | ax | 2903 | 69.2 | 5734 | -19 | 8451 |
| b-maltose | red. C3-H | H-C-O-H | ax | 2909 | 71.5 | 5732 | -20 | 8444 |
| b-maltose | red. C4-H | H-C-O-C(glyc) | ax | 2949 | 71.6 | 5786 | +33 | 8528 |
| b-maltose | red. C5-H | H-C-O-C (ring O5) | ax | 2904 | 72.2 | 5723 | -29 | 8431 |
| b-maltose | red. C6-Ha | CH2-OH | CH2 | 2968 | 67.3 | 5829 | +77 | 8602 |
| b-maltose | red. C6-Hb | CH2-OH | CH2 | 2988 | 66.7 | 5859 | +106 | 8649 |
| a-maltose | non-red. C1-H | anomeric H-C(O5)(O-C(glyc)) | eq | 2944 | 71.3 | 5780 | +27 | 8524 |
| a-maltose | non-red. C2-H | H-C-O-H | ax | 2881 | 73.4 | 5687 | -66 | 8373 |
| a-maltose | non-red. C3-H | H-C-O-H | ax | 2910 | 70.7 | 5736 | -16 | 8451 |
| a-maltose | non-red. C4-H | H-C-O-H | ax | 2887 | 72.4 | 5699 | -54 | 8392 |
| a-maltose | non-red. C5-H | H-C-O-C (ring O5) | ax | 2912 | 73.1 | 5730 | -22 | 8441 |
| a-maltose | non-red. C6-Ha | CH2-OH | CH2 | 2998 | 65.1 | 5878 | +125 | 8681 |
| a-maltose | non-red. C6-Hb | CH2-OH | CH2 | 2915 | 71.1 | 5742 | -10 | 8462 |
| a-maltose | red. C1-H | anomeric H-C(O5)(OH) | eq | 2920 | 68.0 | 5761 | +8 | 8498 |
| a-maltose | red. C2-H | H-C-O-H | ax | 2984 | 65.4 | 5858 | +106 | 8647 |
| a-maltose | red. C3-H | H-C-O-H | ax | 2914 | 69.8 | 5745 | -7 | 8469 |
| a-maltose | red. C4-H | H-C-O-C(glyc) | ax | 2938 | 72.3 | 5768 | +15 | 8500 |
| a-maltose | red. C5-H | H-C-O-C (ring O5) | ax | 2935 | 70.0 | 5773 | +21 | 8510 |
| a-maltose | red. C6-Ha | CH2-OH | CH2 | 2935 | 69.3 | 5776 | +24 | 8519 |
| a-maltose | red. C6-Hb | CH2-OH | CH2 | 3025 | 65.2 | 5913 | +161 | 8734 |

## 2. 유형별 평균 (모든 분자, Δ cm⁻¹, ± 표준편차, 개수)

| 유형 | 축 C-H | 적도 C-H | 전체 |
|---|---|---|---|
| CH2-OH | – | – | +77 ± 44 (n=12) |
| H-C-O-C (ring O5) | -14 ± 22 (n=6) | – | -14 ± 22 (n=6) |
| H-C-O-C(glyc) | +24 ± 9 (n=2) | – | +24 ± 9 (n=2) |
| H-C-O-H | -15 ± 50 (n=16) | – | -15 ± 50 (n=16) |
| anomeric H-C(O5)(O-C(glyc)) | – | +26 ± 2 (n=2) | +26 ± 2 (n=2) |
| anomeric H-C(O5)(OH) | -155 ± 14 (n=2) | +6 ± 3 (n=2) | -75 ± 81 (n=4) |
| (방향만) ax | | | -22 ± 57 (n=26) |
| (방향만) eq | | | +16 ± 10 (n=4) |
| (방향만) CH2 | | | +77 ± 44 (n=12) |

## 3. 포도당 vs 말토오스 (그램당, 아노머 비율 반영, C-H당 FWHM 60 cm⁻¹ 가정)

- 말토오스가 더 큰 구간 (차이 > 최대의 50 %): 5739–5786 cm⁻¹ (Δ -13~+34)
- 포도당이 더 큰 구간: 5559–5594 cm⁻¹ (Δ -193~-158), 5820–5857 cm⁻¹ (Δ +68~+105)
- 말토오스 글리코시드 C-H(비환원 C1-H, 환원 C4-H) 밴드 중심: b-maltose non-red. C1-H 5776 (Δ +24), b-maltose red. C4-H 5786 (Δ +33), a-maltose non-red. C1-H 5780 (Δ +27), a-maltose red. C4-H 5768 (Δ +15)
- 최대 상대 차이 |Δ스펙트럼|/최대흡광: 0.433

## 4. xTB → DFT 보정 (1차 배음)

DFT = 2102 + 0.674 × xTB  (쌍 9개, 잔차 RMSE 27 cm⁻¹)

| 분자 | C-H | 유형 | xTB | DFT | 보정 xTB |
|---|---|---|---|---|---|
| b-glucose | C1-H | anomeric H-C(O5)(OH) | 5171 | 5595 | 5586 |
| b-glucose | C2-H | H-C-O-H | 5484 | 5754 | 5796 |
| b-glucose | C4-H | H-C-O-H | 5437 | 5727 | 5765 |
| b-glucose | C5-H | H-C-O-C (ring O5) | 5400 | 5723 | 5739 |
| b-glucose | C6-Ha | CH2-OH | 5689 | 5917 | 5934 |
| a-glucose | C1-H | anomeric H-C(O5)(OH) | 5499 | 5836 | 5806 |
| a-glucose | C2-H | H-C-O-H | 5650 | 5932 | 5908 |
| me-a-glucoside | C1-H | anomeric H-C(O5)(O-C(glyc)) | 5447 | 5802 | 5771 |
| 4-OMe-b-glucose | C4-H | H-C-O-C(glyc) | 5487 | 5818 | 5798 |

## 5. 파장으로 환산한 결과 (보정값)

- 말토오스가 더 큰 구간: 5739–5786 cm⁻¹ = 1728–1742 nm
- 포도당이 더 큰 구간: 5559–5594 cm⁻¹ = 1788–1799 nm, 5820–5857 cm⁻¹ = 1707–1718 nm