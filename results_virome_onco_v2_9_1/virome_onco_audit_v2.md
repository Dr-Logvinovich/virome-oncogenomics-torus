# Zero-RNG Toric Virome & Oncogenomics Audit (v2.9.1)

> **Methodology:** 100% deterministic, zero-RNG Klein $V_4$ & $V_4$-closed $496 \oplus 3600$ screen engine (894/1024 dimensions inside $W_{-+}$) · 100/100 startup self-tests passed · Human Reference Arrow: `/Volumes/Мои документы/Virusu/human_t2t_euch_ref_v2_8.npz`.

## 1. Theorem V1 — 33 Real NCBI RefSeq Viruses on the Codon Torus

| Virus | Class | Mode | Len (bp) | GC % | ‖W_-+‖²_unb | z_iid(W_-+) | PR2-ok (unb) % | 3600-D unb % (raw) | CpG-free % | cos_raw | cos_unb (z≥3) | Best Human Chr (cos) | Top-1 Motif |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| **HIV1_HXB2** | retrovirus_int | integrating | 9,181 | 42.1% | 6.81e-06 | 4.5 | 61.75% | 24.6% (74.7%) | 92.3% | +0.317 | +0.706 | chrX (+0.335) | `AAT-TTT` (2.9%) |
| **HIV2_BEN** | retrovirus_int | integrating | 10,359 | 45.7% | 8.95e-06 | 7.0 | 58.02% | 42.7% (75.1%) | 88.7% | +0.352 | +0.674 | chr11 (+0.365) | `TTC-CTG` (2.6%) |
| **HTLV1** | retrovirus_int | onco_int | 8,507 | 53.5% | 7.79e-06 | 5.1 | 71.36% | 23.9% (73.9%) | 80.4% | +0.273 | +0.593 | chr11 (+0.280) | `CAG-GAG` (2.1%) |
| **HTLV2** | retrovirus_int | integrating | 8,952 | 53.8% | 1.45e-05 | 9.6 | 89.21% | 47.9% (73.8%) | 78.4% | +0.238 | +0.405 | chr22 (+0.248) | `CTC-CCC` (2.2%) |
| **SIV** | retrovirus_int | integrating | 9,623 | 44.4% | 5.28e-06 | 3.8 | 46.39% | 26.3% (76.9%) | 88.2% | +0.272 | +0.657 | chrX (+0.286) | `TAT-CCC` (1.7%) |
| **MoMLV** | retrovirus_int | onco_int | 8,332 | 53.0% | 9.60e-06 | 6.1 | 68.54% | 25.6% (72.3%) | 79.8% | +0.195 | +0.396 | chr11 (+0.211) | `CCC-CTT` (1.5%) |
| **HPV16** | dna_onco_small | onco_int | 7,906 | 36.5% | 1.05e-05 | 5.3 | 100.00%* | 35.7% (74.3%) | 92.5% | +0.187 | +0.373 | chr13 (+0.206) | `TTT-TAA` (2.4%) |
| **HPV18** | dna_onco_small | onco_int | 7,857 | 40.4% | 1.45e-05 | 7.8 | 100.00%* | 17.1% (65.2%) | 90.2% | +0.214 | +0.382 | chr13 (+0.246) | `TTT-AAA` (2.1%) |
| **HPV6b** | dna_benign_ctrl | non_onco | 7,902 | 40.9% | 1.17e-05 | 6.2 | 100.00%* | 23.1% (69.9%) | 89.1% | +0.161 | +0.310 | chr13 (+0.177) | `TTT-AAA` (1.2%) |
| **HPV11** | dna_benign_ctrl | non_onco | 7,931 | 41.0% | 9.44e-06 | 5.1 | 100.00%* | 17.8% (71.1%) | 88.1% | +0.174 | +0.362 | chrY (+0.191) | `TTT-AAA` (1.3%) |
| **MCPyV** | dna_onco_small | onco_int | 5,387 | 40.7% | 1.32e-05 | 4.8 | 99.22% | 23.1% (72.9%) | 93.2% | +0.292 | +0.617 | chr11 (+0.309) | `CCT-CTG` (2.9%) |
| **SV40** | dna_onco_small | onco_model | 5,243 | 40.8% | 2.05e-05 | 7.1 | 100.00%* | 42.1% (73.6%) | 94.6% | +0.347 | +0.630 | chr13 (+0.368) | `TTT-AAA` (1.7%) |
| **JCV** | dna_benign_ctrl | non_onco | 5,130 | 40.4% | 1.21e-05 | 4.1 | 100.00%* | 9.2% (71.5%) | 93.1% | +0.371 | +0.825 | chr5 (+0.380) | `CAA-AAA` (2.7%) |
| **BKV** | dna_benign_ctrl | non_onco | 5,153 | 39.5% | 1.24e-05 | 4.2 | 89.51% | 10.8% (71.6%) | 95.0% | +0.342 | +0.754 | chr13 (+0.358) | `AAT-ATT` (1.1%) |
| **HBV** | hepatitis_hcc | onco_int | 3,182 | 48.5% | 9.27e-06 | 2.2 | 51.66% | 38.2% (82.1%) | 82.0% | +0.188 | unresolved | chrX (+0.201) | `TTT-GCT` (1.4%) |
| **HCV** | hepatitis_hcc | onco_rna | 9,646 | 58.2% | 7.05e-07 | 0.5 | 22.20% | 0.0% (77.9%) | 69.3% | +0.253 | unresolved | chr11 (+0.264) | `TCC-TGG` (2.0%) |
| **EBV_HHV4** | herpes_large | onco_episome | 171,823 | 59.5% | 1.36e-05 | 170.8 | 86.52% | 44.3% (48.4%) | 60.0% | +0.431 | +0.453 | chr20 (+0.454) | `GCC-CCC` (4.3%) |
| **KSHV_HHV8** | herpes_large | onco_episome | 137,969 | 53.8% | 5.36e-06 | 63.1 | 94.41% | 26.3% (41.6%) | 73.8% | +0.504 | +0.582 | chr21 (+0.514) | `CCC-AGG` (1.5%) |
| **HSV1** | herpes_large | non_onco | 152,222 | 68.3% | 4.56e-06 | 40.5 | 100.00%* | 50.3% (59.9%) | 49.7% | +0.055 | +0.064 | chr19 (+0.070) | `CCC-CCG` (3.2%) |
| **HCMV_HHV5** | herpes_large | non_onco | 235,646 | 57.5% | 6.00e-06 | 113.3 | 100.00%* | 37.1% (44.5%) | 51.7% | +0.081 | +0.088 | chr18 (+0.087) | `CGC-GCC` (2.6%) |
| **AAV2** | vector_dna | site_int | 4,679 | 53.8% | 1.15e-05 | 4.1 | 100.00%* | 58.5% (82.2%) | 69.9% | +0.211 | +0.499 | chr11 (+0.218) | `CAC-GTG` (1.3%) |
| **HAdV5** | vector_dna | onco_model | 35,937 | 55.2% | 6.02e-06 | 17.7 | 93.50% | 33.8% (62.3%) | 68.1% | +0.303 | +0.444 | chr21 (+0.313) | `GTG-CGC` (2.2%) |
| **SARS_CoV_2** | rna_acute_ctrl | non_int_rna | 29,903 | 38.0% | 8.69e-06 | 18.8 | 80.35% | 29.0% (57.4%) | 91.3% | +0.242 | +0.338 | chr13 (+0.270) | `TTG-TTA` (3.6%) |
| **SARS_CoV_1** | rna_acute_ctrl | non_int_rna | 29,751 | 40.8% | 6.02e-06 | 13.7 | 69.03% | 11.5% (55.5%) | 88.5% | +0.360 | +0.556 | chr13 (+0.392) | `ATT-GTT` (1.5%) |
| **MERS_CoV** | rna_acute_ctrl | non_int_rna | 30,119 | 41.2% | 4.18e-06 | 9.9 | 86.58% | 20.5% (64.8%) | 86.3% | +0.260 | +0.447 | chr13 (+0.279) | `TAT-TAA` (1.6%) |
| **Ebola_Zaire** | rna_acute_ctrl | non_int_rna | 18,959 | 41.1% | 4.82e-06 | 7.0 | 100.00%* | 33.8% (72.9%) | 87.7% | +0.216 | +0.416 | chrX (+0.241) | `AAG-AAA` (1.7%) |
| **Measles** | rna_acute_ctrl | non_int_rna | 15,894 | 47.4% | 4.91e-06 | 6.4 | 67.46% | 2.8% (67.1%) | 78.0% | +0.254 | +0.520 | chr11 (+0.273) | `GAT-ATC` (1.0%) |
| **Poliovirus** | rna_acute_ctrl | non_int_rna | 7,440 | 46.3% | 7.58e-06 | 4.4 | 62.63% | 22.4% (75.3%) | 83.6% | +0.212 | +0.492 | chrX (+0.223) | `CAA-ACA` (1.6%) |
| **Zika** | rna_acute_ctrl | non_int_rna | 10,794 | 50.9% | 1.23e-05 | 10.1 | 61.99% | 19.3% (63.5%) | 82.4% | +0.406 | +0.685 | chr11 (+0.417) | `TCC-ACA` (2.0%) |
| **Dengue2** | rna_acute_ctrl | non_int_rna | 10,723 | 45.8% | 1.31e-05 | 10.4 | 78.67% | 15.1% (61.2%) | 86.7% | +0.407 | +0.676 | chrY (+0.426) | `AGG-AAA` (1.5%) |
| **Phage_Lambda** | phage_outgroup | phage_int | 48,502 | 49.9% | 1.54e-05 | 61.9 | 95.87% | 29.4% (43.8%) | 72.8% | +0.181 | +0.208 | chrY (+0.200) | `GCC-GGA` (2.7%) |
| **Phage_T4** | phage_outgroup | phage_lytic | 168,903 | 35.3% | 9.69e-06 | 110.1 | 87.07% | 44.6% (50.2%) | 89.9% | +0.153 | +0.165 | chr4 (+0.186) | `TCT-TTA` (2.3%) |
| **Phage_PhiX174** | phage_outgroup | phage_ssdna | 5,386 | 44.8% | 9.48e-06 | 3.8 | 54.36% | 33.8% (78.2%) | 75.3% | +0.104 | +0.251 | chr13 (+0.117) | `AAG-CTG` (1.6%) |

## 2. Theorem V4 — 28 Real Human Cancer Driver, Fragile-Site & Control Loci (NCBI RefSeqGene)

| Locus | Chr | Role | Len (bp) | GC % | ‖W_-+‖²_unb | PR2-ok (unb) % | Non-Env 3600-D % | CpG-free % | cos(W_-+, Human) | Pal d=6 Share % |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **TP53** | chr17 | tumor_suppressor | 32,772 | 49.0% | 6.86e-05 | 85.51% | 73.2% | 88.7% | +0.728 | 10.27% |
| **TP63** | chr3 | tumor_suppressor | 272,853 | 37.1% | 9.84e-06 | 100.00% | 36.4% | 88.6% | +0.929 | 10.65% |
| **TP73** | chr1 | tumor_suppressor | 90,637 | 59.1% | 3.26e-05 | 88.97% | 39.7% | 79.4% | +0.692 | 3.60% |
| **BRCA1** | chr17 | tumor_suppressor | 193,689 | 45.1% | 3.88e-05 | 99.55% | 66.1% | 89.2% | +0.828 | 9.47% |
| **BRCA2** | chr13 | tumor_suppressor | 91,193 | 38.6% | 1.63e-05 | 100.00% | 50.8% | 89.7% | +0.882 | 8.96% |
| **PTEN** | chr10 | tumor_suppressor | 115,492 | 36.2% | 1.16e-05 | 96.98% | 39.1% | 91.2% | +0.861 | 9.19% |
| **RB1** | chr13 | tumor_suppressor | 185,144 | 37.0% | 1.07e-05 | 96.56% | 37.5% | 89.2% | +0.902 | 4.42% |
| **APC** | chr5 | tumor_suppressor | 163,719 | 37.8% | 1.39e-05 | 98.02% | 45.8% | 90.3% | +0.911 | 8.62% |
| **VHL** | chr3 | tumor_suppressor | 17,444 | 47.0% | 7.23e-05 | 96.17% | 70.5% | 89.7% | +0.690 | 7.83% |
| **MLH1** | chr3 | dna_repair | 79,540 | 41.7% | 2.05e-05 | 100.00% | 51.6% | 89.2% | +0.893 | 7.19% |
| **CDKN2A** | chr9 | tumor_suppressor | 33,740 | 41.8% | 8.13e-06 | 91.05% | 31.0% | 87.6% | +0.690 | 7.94% |
| **ATM** | chr11 | dna_repair | 153,268 | 37.5% | 1.31e-05 | 99.47% | 47.5% | 90.3% | +0.914 | 9.99% |
| **MYC** | chr8 | oncogene_viral | 14,518 | 48.4% | 9.88e-06 | 100.00% | 24.5% | 84.9% | +0.588 | 5.52% |
| **TERT** | chr5 | oncogene_viral | 48,881 | 58.2% | 4.32e-05 | 66.62% | 59.3% | 75.5% | +0.525 | 6.15% |
| **KRAS** | chr12 | oncogene | 53,215 | 37.1% | 1.37e-05 | 98.03% | 44.4% | 91.7% | +0.793 | 8.06% |
| **NRAS** | chr1 | oncogene | 19,438 | 39.0% | 1.35e-05 | 100.00% | 49.4% | 90.5% | +0.714 | 3.42% |
| **HRAS** | chr11 | oncogene | 10,309 | 65.9% | 5.96e-05 | 62.81% | 65.1% | 73.5% | +0.397 | 5.17% |
| **PIK3CA** | chr3 | oncogene_apobec | 98,571 | 35.8% | 1.25e-05 | 100.00% | 42.0% | 91.2% | +0.847 | 6.33% |
| **EGFR** | chr7 | oncogene_ecdna | 244,589 | 45.1% | 1.31e-05 | 100.00% | 23.8% | 84.2% | +0.893 | 6.85% |
| **BRAF** | chr7 | oncogene_uv | 212,438 | 38.1% | 1.72e-05 | 96.02% | 50.2% | 90.0% | +0.908 | 10.50% |
| **ERBB2** | chr17 | oncogene | 47,523 | 52.0% | 3.26e-05 | 100.00% | 53.9% | 88.5% | +0.808 | 8.22% |
| **BCL2** | chr18 | oncogene | 203,035 | 42.9% | 1.23e-05 | 98.76% | 31.4% | 86.4% | +0.939 | 9.37% |
| **ALK** | chr2 | oncogene | 735,793 | 43.5% | 1.21e-05 | 99.83% | 26.4% | 83.8% | +0.936 | 8.89% |
| **FHIT** | chr3 | fragile_viral | 1,511,182 | 38.9% | 9.94e-06 | 99.29% | 30.4% | 86.3% | +0.955 | 12.49% |
| **GAPDH** | chr12 | control_hk | 10,880 | 56.5% | 1.81e-05 | 97.49% | 37.4% | 80.5% | +0.541 | 4.22% |
| **ACTB** | chr7 | control_hk | 10,454 | 58.4% | 3.22e-05 | 98.28% | 52.5% | 84.1% | +0.586 | 4.02% |
| **HBB** | chr11 | control_nononco | 81,706 | 39.5% | 9.94e-06 | 100.00% | 33.2% | 88.4% | +0.815 | 2.53% |
| **APOE** | chr19 | control_gc19 | 10,612 | 58.0% | 4.06e-05 | 97.18% | 57.2% | 85.5% | +0.686 | 8.11% |

## 3. Theorems V3 & V3.3 — Unbiased COSMIC Mutagenesis Jacobian & Exact De Novo vs. Substrate-Coupled Split

| Locus | SBS1 Unb Chiral % | SBS1 De Novo χ_-- | SBS1 Sub χ_-- (×10⁴) | APOBEC-SBS2(+) Unb Chiral % | APOBEC De Novo χ_--^unb (×10⁴) | APOBEC De Novo % | ε_crit (1/4N) | UV-SBS7a(+) Unb Chiral % | Tobacco-SBS4(+) Unb Chiral % |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **TP53** | 18.39% | 2.43e-35 | 0.0821 | 48.73% | 0.163 | 75.9% | 0.594 | 46.21% | 37.25% |
| **TP63** | 3.60% | 6.84e-36 | 0.0005 | 48.70% | 0.115 | 97.8% | 0.276 | 45.65% | 41.69% |
| **TP73** | 8.84% | 1.12e-34 | 0.0202 | 49.17% | 0.111 | 93.4% | 0.474 | 43.66% | 41.04% |
| **BRCA1** | 2.12% | 7.17e-35 | 0.0032 | 47.67% | 0.139 | 97.4% | 0.305 | 44.06% | 34.92% |
| **BRCA2** | 5.83% | 3.58e-35 | 0.0027 | 47.99% | 0.119 | 94.1% | 0.490 | 44.97% | 36.33% |
| **PTEN** | 8.74% | 1.45e-35 | 0.0018 | 48.60% | 0.105 | 92.6% | 0.450 | 46.80% | 42.57% |
| **RB1** | 11.44% | 1.57e-35 | 0.0025 | 49.13% | 0.113 | 95.7% | 0.335 | 45.52% | 40.00% |
| **APC** | 9.15% | 1.28e-35 | 0.0025 | 48.48% | 0.108 | 95.9% | 0.363 | 46.78% | 36.68% |
| **VHL** | 16.47% | 1.92e-34 | 0.0710 | 46.56% | 0.163 | 83.0% | 0.955 | 44.51% | 36.31% |
| **MLH1** | 3.54% | 1.75e-35 | 0.0017 | 48.86% | 0.123 | 94.1% | 0.503 | 45.24% | 36.81% |
| **CDKN2A** | 24.79% | 8.09e-35 | 0.0332 | 48.34% | 0.090 | 86.3% | 0.890 | 48.22% | 49.31% |
| **ATM** | 6.52% | 1.18e-35 | 0.0015 | 48.77% | 0.116 | 96.6% | 0.381 | 43.75% | 37.01% |
| **MYC** | 25.43% | 1.96e-34 | 0.0774 | 49.05% | 0.100 | 79.1% | 1.308 | 48.34% | 43.63% |
| **TERT** | 30.93% | 3.65e-34 | 0.3498 | 47.45% | 0.158 | 66.8% | 0.530 | 44.26% | 45.84% |
| **KRAS** | 19.60% | 3.22e-35 | 0.0093 | 48.38% | 0.096 | 88.7% | 0.689 | 46.60% | 39.68% |
| **NRAS** | 21.23% | 1.96e-35 | 0.0131 | 49.36% | 0.110 | 81.4% | 1.043 | 46.41% | 39.75% |
| **HRAS** | 37.58% | 7.28e-34 | 0.7594 | 50.55% | 0.178 | 55.5% | 0.665 | 45.95% | 43.74% |
| **PIK3CA** | 7.72% | 1.26e-35 | 0.0021 | 48.43% | 0.113 | 94.6% | 0.478 | 44.64% | 38.54% |
| **EGFR** | 6.11% | 2.13e-35 | 0.0029 | 49.13% | 0.114 | 97.7% | 0.289 | 45.79% | 37.73% |
| **BRAF** | 6.09% | 2.13e-35 | 0.0022 | 48.13% | 0.113 | 96.1% | 0.319 | 44.30% | 37.53% |
| **ERBB2** | 12.29% | 4.10e-35 | 0.0113 | 48.05% | 0.117 | 91.2% | 0.666 | 41.81% | 36.20% |
| **BCL2** | 5.53% | 6.51e-36 | 0.0017 | 49.10% | 0.114 | 95.9% | 0.318 | 45.52% | 40.59% |
| **ALK** | 2.90% | 8.25e-36 | 0.0005 | 49.15% | 0.127 | 99.1% | 0.164 | 45.26% | 37.58% |
| **FHIT** | 1.83% | 4.73e-36 | 0.0002 | 49.21% | 0.117 | 99.1% | 0.119 | 45.86% | 41.07% |
| **GAPDH** | 31.03% | 1.65e-34 | 0.1293 | 47.20% | 0.117 | 78.0% | 1.452 | 48.84% | 44.60% |
| **ACTB** | 21.27% | 2.41e-34 | 0.1377 | 50.50% | 0.117 | 72.5% | 1.187 | 45.56% | 35.44% |
| **HBB** | 20.40% | 9.87e-36 | 0.0032 | 48.87% | 0.128 | 93.7% | 0.492 | 46.36% | 46.15% |
| **APOE** | 20.06% | 2.56e-34 | 0.0864 | 48.00% | 0.138 | 76.3% | 1.191 | 42.49% | 36.16% |

## 4. Theorem V2 — Deterministic In-Silico Concatenation Proviral Model into Real Human Onco-Loci

| Virus -> Host Locus | Virus % | Host PR2-ok % | +Strand PR2-ok % | +Strand E/floor | HH E/floor | cos(ΔW_-- raw) | cos(ΔB_-- count-space) | cos(ΔB_-+ count-space) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `HBV` -> `TERT` | 6.1% | 66.62% | 66.81% | 4.0x | 3.6x | -0.524 | -0.9955 | +0.9970 |
| `HBV` -> `MYC` | 18.0% | 100.00% | 100.00% | -0.1x | -0.4x | -0.701 | -0.9966 | +0.9970 |
| `HPV16` -> `MYC` | 35.3% | 100.00% | 100.00% | -0.1x | -0.5x | -0.328 | -0.9980 | +0.9992 |
| `HPV16` -> `FHIT` | 0.5% | 99.29% | 99.31% | 0.4x | 0.4x | -0.982 | -0.9975 | +0.9990 |
| `HPV18` -> `MYC` | 35.1% | 100.00% | 100.00% | -0.0x | -0.5x | -0.348 | -0.9981 | +0.9991 |
| `HPV6b` -> `MYC` | 35.2% | 100.00% | 100.00% | -0.1x | -0.5x | -0.344 | -0.9979 | +0.9990 |
| `HIV1_HXB2` -> `MYC` | 38.7% | 100.00% | 90.55% | 0.1x | -0.6x | -0.347 | -0.9983 | +0.9993 |
| `HIV1_HXB2` -> `TP53` | 21.9% | 85.51% | 86.03% | 1.2x | 0.6x | -0.232 | -0.9982 | +0.9993 |
| `HIV1_HXB2` -> `APOE` | 46.4% | 97.18% | 93.95% | 0.1x | -0.6x | -0.121 | -0.9985 | +0.9994 |
| `HTLV1` -> `TP53` | 20.6% | 85.51% | 85.39% | 1.3x | 0.7x | -0.249 | -0.9980 | +0.9991 |
| `MCPyV` -> `TP53` | 14.1% | 85.51% | 85.79% | 1.3x | 0.9x | -0.402 | -0.9970 | +0.9986 |
| `AAV2` -> `APOE` | 30.6% | 97.18% | 98.61% | 0.0x | -0.4x | -0.361 | -0.9941 | +1.0000 |

## 5. Theorem V5 — Exact T2T-CHM13v2.0 ±50 kb Chromosomal Neighborhoods Around Cancer & Viral Hotspots

| Locus | T2T Coord | Strand | Win GC % | ‖W_-+‖²_unb | PR2-ok % | Non-Env 3600-D % | Pal d=6 % | cos(vs Chr Euch) | Top Matching Virus (cos) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **TP53** | `chr17:7.58 Mb` | `-` | 51.3% | 4.83e-05 | 93.18% | 68.7% | 10.2% | +0.914 | `KSHV_HHV8` (+0.362) |
| **BRCA1** | `chr17:43.98 Mb` | `-` | 44.0% | 3.65e-05 | 99.56% | 68.5% | 9.7% | +0.908 | `KSHV_HHV8` (+0.375) |
| **ERBB2** | `chr17:40.57 Mb` | `+` | 54.5% | 3.14e-05 | 100.00% | 50.5% | 8.3% | +0.899 | `EBV_HHV4` (+0.489) |
| **TP63** | `chr3:192.58 Mb` | `+` | 37.1% | 9.87e-06 | 100.00% | 46.4% | 8.8% | +0.880 | `KSHV_HHV8` (+0.424) |
| **VHL** | `chr3:10.14 Mb` | `+` | 47.2% | 5.64e-05 | 100.00% | 70.3% | 9.5% | +0.722 | `KSHV_HHV8` (+0.358) |
| **MLH1** | `chr3:37.03 Mb` | `+` | 41.1% | 2.06e-05 | 100.00% | 56.1% | 8.2% | +0.878 | `KSHV_HHV8` (+0.429) |
| **PIK3CA** | `chr3:182.00 Mb` | `+` | 35.9% | 1.28e-05 | 100.00% | 49.8% | 6.3% | +0.858 | `KSHV_HHV8` (+0.439) |
| **FHIT** | `chr3:60.54 Mb` | `-` | 38.4% | 9.52e-06 | 93.58% | 45.1% | 7.8% | +0.882 | `KSHV_HHV8` (+0.436) |
| **TP73** | `chr1:3.20 Mb` | `+` | 58.8% | 3.35e-05 | 84.35% | 47.4% | 2.4% | +0.693 | `EBV_HHV4` (+0.516) |
| **NRAS** | `chr1:114.72 Mb` | `-` | 39.3% | 1.43e-05 | 99.77% | 51.6% | 7.9% | +0.891 | `KSHV_HHV8` (+0.446) |
| **BRCA2** | `chr13:31.57 Mb` | `+` | 38.6% | 1.65e-05 | 100.00% | 55.7% | 8.6% | +0.846 | `KSHV_HHV8` (+0.434) |
| **RB1** | `chr13:47.61 Mb` | `+` | 35.5% | 1.14e-05 | 93.45% | 48.1% | 3.2% | +0.868 | `KSHV_HHV8` (+0.427) |
| **PTEN** | `chr10:88.80 Mb` | `+` | 35.5% | 1.25e-05 | 96.79% | 50.0% | 8.9% | +0.850 | `KSHV_HHV8` (+0.450) |
| **APC** | `chr5:113.28 Mb` | `+` | 37.0% | 1.53e-05 | 97.91% | 53.7% | 8.5% | +0.867 | `KSHV_HHV8` (+0.428) |
| **TERT** | `chr5:1.18 Mb` | `-` | 55.9% | 2.86e-05 | 78.97% | 49.5% | 5.8% | +0.643 | `EBV_HHV4` (+0.481) |
| **CDKN2A** | `chr9:22.00 Mb` | `-` | 40.3% | 9.09e-06 | 98.39% | 44.1% | 5.6% | +0.892 | `KSHV_HHV8` (+0.468) |
| **ATM** | `chr11:108.30 Mb` | `+` | 36.8% | 1.27e-05 | 100.00% | 53.6% | 10.2% | +0.869 | `KSHV_HHV8` (+0.422) |
| **HRAS** | `chr11:0.58 Mb` | `-` | 59.7% | 3.91e-05 | 92.91% | 56.3% | 7.5% | +0.746 | `KSHV_HHV8` (+0.428) |
| **HBB** | `chr11:5.32 Mb` | `-` | 39.0% | 9.36e-06 | 100.00% | 44.1% | 3.1% | +0.850 | `KSHV_HHV8` (+0.414) |
| **MYC** | `chr8:128.87 Mb` | `+` | 44.5% | 2.55e-05 | 100.00% | 59.1% | 9.9% | +0.857 | `KSHV_HHV8` (+0.389) |
| **KRAS** | `chr12:25.10 Mb` | `-` | 38.0% | 1.44e-05 | 96.01% | 52.8% | 8.0% | +0.884 | `KSHV_HHV8` (+0.423) |
| **GAPDH** | `chr12:6.55 Mb` | `+` | 49.1% | 3.27e-05 | 100.00% | 64.0% | 9.8% | +0.861 | `KSHV_HHV8` (+0.397) |
| **EGFR** | `chr7:55.30 Mb` | `+` | 45.1% | 1.41e-05 | 99.07% | 34.2% | 4.4% | +0.829 | `KSHV_HHV8` (+0.478) |
| **BRAF** | `chr7:142.14 Mb` | `-` | 37.0% | 1.59e-05 | 98.47% | 54.1% | 9.2% | +0.864 | `KSHV_HHV8` (+0.435) |
| **ACTB** | `chr7:5.65 Mb` | `-` | 54.7% | 7.37e-05 | 97.51% | 70.6% | 10.3% | +0.752 | `KSHV_HHV8` (+0.348) |
| **BCL2** | `chr18:63.43 Mb` | `-` | 42.3% | 1.26e-05 | 90.10% | 43.7% | 5.0% | +0.857 | `KSHV_HHV8` (+0.458) |
| **ALK** | `chr2:29.60 Mb` | `-` | 43.0% | 1.16e-05 | 99.74% | 36.5% | 5.1% | +0.867 | `KSHV_HHV8` (+0.481) |
| **APOE** | `chr19:47.73 Mb` | `+` | 53.0% | 5.98e-05 | 100.00% | 68.6% | 10.0% | +0.963 | `KSHV_HHV8` (+0.376) |

## 6. Theorem V6 — Multi-Lag Toric Decay & Human Impedance Across Lags (d = 3, 6, 30, 171 bp)

| Virus | Class | d=3 ‖W_-+‖² (cos_Hum) | d=6 ‖W_-+‖² (Pal %) | d=30 ‖W_-+‖² (cos_Hum) | d=171 ‖W_-+‖² (cos_Hum) |
|---|---|---|---|---|---|
| **HIV1_HXB2** | retrovirus_int | 6.81e-06 (cos +0.32) | -2.67e-07 (Pal 3.3%) | -1.79e-08 (cos -0.00) | -2.43e-06 (cos +0.00) |
| **HIV2_BEN** | retrovirus_int | 8.95e-06 (cos +0.35) | 1.76e-06 (Pal 5.2%) | 5.92e-06 (cos +0.06) | 2.93e-06 (cos -0.03) |
| **HTLV1** | retrovirus_int | 7.79e-06 (cos +0.27) | 1.41e-06 (Pal 4.4%) | -8.63e-07 (cos -0.01) | -2.34e-06 (cos -0.01) |
| **HTLV2** | retrovirus_int | 1.45e-05 (cos +0.24) | 8.27e-06 (Pal 3.6%) | 4.53e-06 (cos -0.02) | 3.98e-06 (cos +0.03) |
| **SIV** | retrovirus_int | 5.28e-06 (cos +0.27) | 3.42e-06 (Pal 2.9%) | 3.40e-06 (cos -0.02) | 1.65e-06 (cos +0.04) |
| **MoMLV** | retrovirus_int | 9.60e-06 (cos +0.20) | 1.25e-06 (Pal 3.2%) | -1.26e-06 (cos -0.01) | -1.77e-06 (cos +0.02) |
| **HPV16** | dna_onco_small | 1.05e-05 (cos +0.19) | 1.60e-07 (Pal 2.6%) | -1.48e-06 (cos -0.00) | -2.49e-06 (cos -0.05) |
| **HPV18** | dna_onco_small | 1.45e-05 (cos +0.21) | -2.26e-06 (Pal 4.1%) | -6.04e-07 (cos +0.11) | -2.35e-06 (cos -0.01) |
| **HPV6b** | dna_benign_ctrl | 1.17e-05 (cos +0.16) | -9.47e-07 (Pal 3.4%) | 1.80e-06 (cos +0.08) | -4.44e-07 (cos -0.02) |
| **HPV11** | dna_benign_ctrl | 9.44e-06 (cos +0.17) | -1.36e-06 (Pal 3.5%) | -2.51e-06 (cos +0.07) | -1.50e-07 (cos -0.02) |
| **MCPyV** | dna_onco_small | 1.32e-05 (cos +0.29) | 1.09e-06 (Pal 6.0%) | -4.53e-07 (cos +0.05) | 3.97e-06 (cos +0.05) |
| **SV40** | dna_onco_small | 2.05e-05 (cos +0.35) | -7.30e-07 (Pal 3.0%) | -7.37e-08 (cos +0.01) | -3.54e-06 (cos +0.04) |
| **JCV** | dna_benign_ctrl | 1.21e-05 (cos +0.37) | -1.56e-06 (Pal 3.4%) | 4.14e-06 (cos -0.01) | 1.20e-06 (cos -0.03) |
| **BKV** | dna_benign_ctrl | 1.24e-05 (cos +0.34) | -4.91e-07 (Pal 3.9%) | -3.03e-06 (cos +0.03) | -2.82e-07 (cos +0.05) |
| **HBV** | hepatitis_hcc | 9.27e-06 (cos +0.19) | -5.51e-06 (Pal 3.9%) | -3.48e-06 (cos -0.06) | -2.32e-06 (cos +0.05) |
| **HCV** | hepatitis_hcc | 7.05e-07 (cos +0.25) | -1.22e-06 (Pal 4.0%) | -3.06e-06 (cos +0.07) | 1.01e-06 (cos +0.00) |
| **EBV_HHV4** | herpes_large | 1.36e-05 (cos +0.43) | 5.19e-06 (Pal 3.7%) | 1.65e-06 (cos -0.01) | 1.89e-06 (cos -0.04) |
| **KSHV_HHV8** | herpes_large | 5.36e-06 (cos +0.50) | 1.07e-06 (Pal 4.3%) | 2.89e-07 (cos +0.03) | 1.10e-07 (cos -0.00) |
| **HSV1** | herpes_large | 4.56e-06 (cos +0.06) | 2.75e-06 (Pal 2.2%) | 5.63e-07 (cos +0.08) | 3.67e-07 (cos +0.02) |
| **HCMV_HHV5** | herpes_large | 6.00e-06 (cos +0.08) | 1.16e-06 (Pal 5.1%) | 1.14e-07 (cos +0.01) | 6.62e-08 (cos -0.02) |
| **AAV2** | vector_dna | 1.15e-05 (cos +0.21) | 6.21e-06 (Pal 3.8%) | 6.25e-07 (cos +0.01) | 1.28e-06 (cos +0.02) |
| **HAdV5** | vector_dna | 6.02e-06 (cos +0.30) | 4.40e-07 (Pal 5.0%) | -2.00e-07 (cos +0.05) | -2.09e-07 (cos +0.01) |
| **SARS_CoV_2** | rna_acute_ctrl | 8.69e-06 (cos +0.24) | 5.05e-07 (Pal 2.7%) | -3.62e-07 (cos +0.04) | -7.43e-07 (cos -0.04) |
| **SARS_CoV_1** | rna_acute_ctrl | 6.02e-06 (cos +0.36) | 1.24e-07 (Pal 3.9%) | 2.88e-07 (cos -0.01) | -7.46e-07 (cos -0.04) |
| **MERS_CoV** | rna_acute_ctrl | 4.18e-06 (cos +0.26) | 4.96e-09 (Pal 4.2%) | -4.01e-07 (cos +0.02) | 1.32e-07 (cos +0.02) |
| **Ebola_Zaire** | rna_acute_ctrl | 4.82e-06 (cos +0.22) | -1.05e-06 (Pal 2.5%) | -1.32e-06 (cos +0.01) | -3.66e-07 (cos +0.07) |
| **Measles** | rna_acute_ctrl | 4.91e-06 (cos +0.25) | 6.38e-07 (Pal 4.4%) | -1.27e-07 (cos +0.01) | 5.60e-07 (cos -0.03) |
| **Poliovirus** | rna_acute_ctrl | 7.58e-06 (cos +0.21) | -9.96e-07 (Pal 3.6%) | 2.89e-06 (cos -0.00) | -3.79e-07 (cos -0.07) |
| **Zika** | rna_acute_ctrl | 1.23e-05 (cos +0.41) | -6.15e-07 (Pal 4.1%) | -6.52e-07 (cos -0.01) | 2.93e-07 (cos +0.01) |
| **Dengue2** | rna_acute_ctrl | 1.31e-05 (cos +0.41) | 1.06e-06 (Pal 2.4%) | -1.55e-06 (cos -0.06) | -5.71e-07 (cos -0.00) |
| **Phage_Lambda** | phage_outgroup | 1.54e-05 (cos +0.18) | 5.60e-07 (Pal 3.0%) | 8.83e-08 (cos +0.01) | 2.29e-08 (cos +0.06) |
| **Phage_T4** | phage_outgroup | 9.69e-06 (cos +0.15) | 1.25e-06 (Pal 4.6%) | -5.73e-08 (cos -0.02) | -2.10e-08 (cos +0.03) |
| **Phage_PhiX174** | phage_outgroup | 9.48e-06 (cos +0.10) | 7.78e-07 (Pal 3.7%) | -6.44e-06 (cos -0.02) | -1.93e-06 (cos -0.01) |

## 7. Theorem V7 — 1-kb Tiled COSMIC Mutagenesis Chiral-Susceptibility Hotspots (χ_-- × 10⁴)

| Locus | Tiles (1 kb) | APOBEC SBS2(+) Med / Peak χ_--×10⁴ (Hotspot) | APOBEC SBS13(+) Peak χ_--×10⁴ (Hotspot) | UV SBS7a(+) Peak χ_--×10⁴ | Tobacco SBS4(+) Peak χ_--×10⁴ |
|---|---:|---|---|---:|---:|
| **TP53** | 32 | 0.27 / **1.12** (`10-11 kb`) | **1.16** (`10-11 kb`) | 7.73 | 32.18 |
| **TP63** | 272 | 0.20 / **2.54** (`107-108 kb`) | **2.24** (`17-18 kb`) | 14.38 | 32.04 |
| **TP73** | 90 | 0.15 / **19.81** (`32-33 kb`) | **17.34** (`32-33 kb`) | 290.54 | 1149.36 |
| **BRCA1** | 193 | 0.21 / **2.40** (`142-143 kb`) | **1.72** (`142-143 kb`) | 11.59 | 17.95 |
| **BRCA2** | 91 | 0.18 / **0.91** (`40-41 kb`) | **1.06** (`40-41 kb`) | 4.69 | 12.45 |
| **PTEN** | 115 | 0.15 / **3.32** (`62-63 kb`) | **2.36** (`62-63 kb`) | 5.67 | 10.26 |
| **RB1** | 185 | 0.20 / **30.00** (`129-130 kb`) | **30.00** (`129-130 kb`) | 82.35 | 61.94 |
| **APC** | 163 | 0.19 / **1.07** (`150-151 kb`) | **1.44** (`150-151 kb`) | 8.27 | 9.87 |
| **VHL** | 17 | 0.16 / **0.56** (`9-10 kb`) | **0.70** (`9-10 kb`) | 3.98 | 11.17 |
| **MLH1** | 79 | 0.21 / **0.75** (`55-56 kb`) | **0.88** (`12-13 kb`) | 4.36 | 18.49 |
| **CDKN2A** | 33 | 0.14 / **1.40** (`22-23 kb`) | **0.91** (`22-23 kb`) | 6.61 | 12.97 |
| **ATM** | 153 | 0.20 / **2.18** (`23-24 kb`) | **2.18** (`23-24 kb`) | 4.53 | 13.35 |
| **MYC** | 14 | 0.12 / **0.66** (`1-2 kb`) | **0.51** (`2-3 kb`) | 6.81 | 8.10 |
| **TERT** | 48 | 0.24 / **19.54** (`40-41 kb`) | **19.79** (`40-41 kb`) | 269.57 | 379.65 |
| **KRAS** | 53 | 0.15 / **1.13** (`22-23 kb`) | **1.22** (`22-23 kb`) | 5.08 | 9.78 |
| **NRAS** | 19 | 0.18 / **1.12** (`11-12 kb`) | **0.75** (`11-12 kb`) | 4.31 | 8.68 |
| **HRAS** | 10 | 0.27 / **30.37** (`9-10 kb`) | **25.94** (`9-10 kb`) | 207.28 | 154.21 |
| **PIK3CA** | 98 | 0.19 / **1.03** (`57-58 kb`) | **1.23** (`32-33 kb`) | 3.24 | 7.44 |
| **EGFR** | 244 | 0.21 / **4.56** (`206-207 kb`) | **4.17** (`206-207 kb`) | 32.10 | 85.62 |
| **BRAF** | 212 | 0.19 / **1.26** (`65-66 kb`) | **1.32** (`65-66 kb`) | 8.42 | 26.60 |
| **ERBB2** | 47 | 0.17 / **1.22** (`26-27 kb`) | **0.92** (`26-27 kb`) | 6.32 | 11.18 |
| **BCL2** | 203 | 0.20 / **13.27** (`125-126 kb`) | **13.67** (`125-126 kb`) | 13.46 | 39.92 |
| **ALK** | 735 | 0.22 / **5.10** (`672-673 kb`) | **4.96** (`672-673 kb`) | 233.21 | 97.41 |
| **FHIT** | 1511 | 0.19 / **11.29** (`395-396 kb`) | **23.30** (`395-396 kb`) | 29.93 | 20.04 |
| **GAPDH** | 10 | 0.19 / **0.38** (`6-7 kb`) | **0.31** (`0-1 kb`) | 5.55 | 11.93 |
| **ACTB** | 10 | 0.25 / **5.57** (`9-10 kb`) | **5.78** (`9-10 kb`) | 7.68 | 8.85 |
| **HBB** | 81 | 0.17 / **5.54** (`3-4 kb`) | **5.36** (`3-4 kb`) | 10.53 | 23.20 |
| **APOE** | 10 | 0.16 / **1.21** (`7-8 kb`) | **0.98** (`7-8 kb`) | 8.78 | 9.91 |

## 8. Theorem V8 — 6-Frame Duplex Stop Barrier & Synonymous Mutation Toric Disruption

| Locus | Chr | Role | 6-Frame Stop % | Longest ORF (codons) | Synonymous Hits | Antisense (-0) Stop Flip % | Mean Δ‖B_-+^3600‖² / Syn (count²) |
|---|---|---|---:|---:|---:|---:|---:|
| **TP53** | chr17 | tumor_suppressor | 4.15% | 186 | 380 | 14.47% | 2.538 |
| **TP63** | chr3 | tumor_suppressor | 5.97% | 213 | 426 | 9.86% | 2.673 |
| **TP73** | chr1 | tumor_suppressor | 2.58% | 451 | 1,065 | 11.74% | 2.368 |
| **BRCA1** | chr17 | tumor_suppressor | 4.82% | 1,103 | 1,906 | 12.01% | 2.669 |
| **BRCA2** | chr13 | tumor_suppressor | 5.77% | 1,517 | 2,622 | 13.16% | 2.674 |
| **PTEN** | chr10 | tumor_suppressor | 6.14% | 159 | 276 | 13.04% | 2.542 |
| **RB1** | chr13 | tumor_suppressor | 5.88% | 646 | 1,145 | 12.93% | 2.673 |
| **APC** | chr5 | tumor_suppressor | 5.95% | 2,143 | 3,798 | 14.80% | 2.634 |
| **VHL** | chr3 | tumor_suppressor | 4.48% | 148 | 325 | 4.92% | 2.593 |
| **MLH1** | chr3 | dna_repair | 5.30% | 132 | 205 | 6.83% | 2.018 |
| **CDKN2A** | chr9 | tumor_suppressor | 5.23% | 160 | 301 | 4.98% | 2.221 |
| **ATM** | chr11 | dna_repair | 5.95% | 751 | 1,307 | 13.24% | 2.675 |
| **MYC** | chr8 | oncogene_viral | 4.28% | 355 | 710 | 11.55% | 2.558 |
| **TERT** | chr5 | oncogene_viral | 2.87% | 821 | 1,958 | 7.35% | 2.620 |
| **KRAS** | chr12 | oncogene | 6.15% | 308 | 780 | 5.77% | 2.426 |
| **NRAS** | chr1 | oncogene | 5.72% | 180 | 337 | 15.43% | 2.549 |
| **HRAS** | chr11 | oncogene | 1.90% | 330 | 817 | 10.28% | 2.735 |
| **PIK3CA** | chr3 | oncogene_apobec | 6.18% | 128 | 215 | 16.74% | 2.674 |
| **EGFR** | chr7 | oncogene_ecdna | 4.66% | 373 | 714 | 1.68% | 2.125 |
| **BRAF** | chr7 | oncogene_uv | 5.83% | 151 | 339 | 8.55% | 2.418 |
| **ERBB2** | chr17 | oncogene | 3.58% | 528 | 1,081 | 10.27% | 2.634 |
| **BCL2** | chr18 | oncogene | 4.92% | 393 | 822 | 13.38% | 1.637 |
| **ALK** | chr2 | oncogene | 5.01% | 522 | 872 | 13.42% | 2.674 |
| **FHIT** | chr3 | fragile_viral | 5.67% | 411 | 701 | 13.12% | 2.682 |
| **GAPDH** | chr12 | control_hk | 3.18% | 302 | 647 | 9.74% | 2.610 |
| **ACTB** | chr7 | control_hk | 2.91% | 167 | 406 | 6.16% | 2.498 |
| **HBB** | chr11 | control_nononco | 5.62% | 372 | 637 | 14.60% | 2.710 |
| **APOE** | chr19 | control_gc19 | 2.87% | 236 | 500 | 6.80% | 2.604 |

## 9. Theorem V9 — Nonparametric Mann–Whitney Group Discrimination (Zero-RNG)

| Contrast | Metric | n1 / n2 | Median 1 vs Median 2 | Rank-Biserial r | Two-Sided p |
|---|---|---:|---|---:|---:|
| Onco/Integrating Viruses vs Acute RNA Controls | `cos_Wmp_vs_human_grand` | 18 / 8 | 0.2625 vs 0.2568 | -0.208 | 4.29e-01 |
| Onco/Integrating Viruses vs Acute RNA Controls | `cos_Wmp_unb_vs_human_grand` | 16 / 8 | 0.4763 vs 0.5056 | -0.109 | 6.97e-01 |
| Onco/Integrating Viruses vs Acute RNA Controls | `cos_Wmp_nonenv3600_vs_human` | 18 / 8 | 0.0642 vs 0.0402 | +0.306 | 2.38e-01 |
| Onco/Integrating Viruses vs Acute RNA Controls | `cos_Wmp_cpgfree_vs_human` | 18 / 8 | 0.2492 vs 0.2264 | -0.056 | 8.49e-01 |
| Onco/Integrating Viruses vs Acute RNA Controls | `cos_Wmp_beyond_M1_vs_human` | 18 / 8 | 0.1415 vs 0.0344 | +0.597 | 1.58e-02 |
| Onco/Integrating Viruses vs Acute RNA Controls | `cos_Wmp_beyond_M2_vs_human` | 18 / 8 | 0.0796 vs -0.0133 | +0.708 | 3.35e-03 |
| Onco/Integrating Viruses vs Acute RNA Controls | `pure_nonenv_3600_unb` | 18 / 8 | 0.3478 vs 0.1994 | +0.653 | 7.63e-03 |
| Onco/Integrating Viruses vs Acute RNA Controls | `z_iid_W_mm` | 18 / 8 | 1.5780 vs 2.9773 | -0.403 | 1.15e-01 |
| Onco/Integrating Viruses vs Acute RNA Controls | `cpg_free_share` | 18 / 8 | 0.8006 vs 0.8648 | -0.222 | 3.97e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cos_Wmp_vs_human_grand` | 4 / 4 | 0.2530 vs 0.2580 | +0.125 | 8.86e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cos_Wmp_unb_vs_human_grand` | 4 / 4 | 0.4994 vs 0.5581 | +0.000 | 1.00e+00 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cos_Wmp_nonenv3600_vs_human` | 4 / 4 | 0.0884 vs 0.0692 | +0.250 | 6.86e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cos_Wmp_cpgfree_vs_human` | 4 / 4 | 0.2298 vs 0.2365 | +0.125 | 8.86e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cos_Wmp_beyond_M1_vs_human` | 4 / 4 | 0.1415 vs 0.1105 | +0.375 | 4.86e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cos_Wmp_beyond_M2_vs_human` | 4 / 4 | 0.0567 vs 0.0533 | +0.125 | 8.86e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `pure_nonenv_3600_unb` | 4 / 4 | 0.2941 vs 0.1431 | +0.750 | 1.14e-01 |
| High-Risk DNA Oncoviruses vs Benign Controls | `z_iid_W_mm` | 4 / 4 | -0.2327 vs -0.2185 | +0.000 | 1.00e+00 |
| High-Risk DNA Oncoviruses vs Benign Controls | `cpg_free_share` | 4 / 4 | 0.9288 vs 0.9112 | +0.250 | 6.86e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cos_Wmp_vs_human_grand` | 30 / 3 | 0.2568 vs 0.1534 | +0.822 | 1.50e-02 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cos_Wmp_unb_vs_human_grand` | 28 / 3 | 0.4954 vs 0.2082 | +0.857 | 1.02e-02 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cos_Wmp_nonenv3600_vs_human` | 30 / 3 | 0.0626 vs 0.0377 | +0.333 | 3.82e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cos_Wmp_cpgfree_vs_human` | 30 / 3 | 0.2424 vs 0.1771 | +0.533 | 1.49e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cos_Wmp_beyond_M1_vs_human` | 30 / 3 | 0.1068 vs 0.0725 | +0.067 | 8.83e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cos_Wmp_beyond_M2_vs_human` | 30 / 3 | 0.0547 vs 0.0180 | +0.222 | 5.71e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `pure_nonenv_3600_unb` | 30 / 3 | 0.2510 vs 0.3375 | -0.444 | 2.35e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `z_iid_W_mm` | 30 / 3 | 1.7237 vs 3.2296 | -0.556 | 1.31e-01 |
| Human Viruses (n=30) vs Bacteriophages (n=3) | `cpg_free_share` | 30 / 3 | 0.8648 vs 0.7530 | +0.222 | 5.71e-01 |
| Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4) | `APOBEC_SBS2_chi_mm_unb_x1e4` | 24 / 4 | 0.1194 vs 0.1463 | -0.354 | 2.91e-01 |
| Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4) | `APOBEC_SBS2_denovo_chi_mm_x1e4` | 24 / 4 | 0.1149 vs 0.1225 | -0.396 | 2.35e-01 |
| Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4) | `APOBEC_SBS2_unb_chiral_share` | 24 / 4 | 0.4871 vs 0.4843 | +0.104 | 7.76e-01 |
| Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4) | `pure_nonenv_3600_unb` | 24 / 4 | 0.4507 vs 0.4498 | -0.042 | 9.24e-01 |
| Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4) | `pal_d6_Wmp_share` | 24 / 4 | 0.0814 vs 0.0412 | +0.667 | 3.51e-02 |
| Length-Matched 100-kb T2T: Drivers vs Controls | `T2T_100kb_NonEnv_3600_unb` | 24 / 4 | 0.4654 vs 0.6504 | -0.479 | 1.44e-01 |
| Length-Matched 100-kb T2T: Drivers vs Controls | `T2T_100kb_Pal_d6_share` | 24 / 4 | 0.0808 vs 0.0991 | -0.417 | 2.09e-01 |
| Length-Matched 100-kb T2T: Drivers vs Controls | `T2T_100kb_cos_vs_chr_euch` | 24 / 4 | 0.8669 vs 0.8552 | +0.083 | 8.25e-01 |

## 10. Theorems V10 & V10.2 — GC- & Dinucleotide-Controlled Viral APOBEC3 TCW Evasion (33 Viruses)

| Virus | Class | GC % | TCW % | TCW O/E (Base) | TCW O/E (Dinuc) | CpG O/E | SBS2(+) Unb Chiral % | SBS2(+) De Novo χ_--^unb (×10⁴) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **HIV1_HXB2** | retrovirus_int | 42.1% | 2.42% | **1.053** | 0.981 | 0.206 | 45.74% | 0.085 |
| **HIV2_BEN** | retrovirus_int | 45.7% | 2.40% | **1.045** | 0.947 | 0.318 | 48.87% | 0.101 |
| **HTLV1** | retrovirus_int | 53.5% | 3.75% | **0.983** | 0.961 | 0.580 | 49.09% | 0.065 |
| **HTLV2** | retrovirus_int | 53.8% | 3.66% | **1.013** | 0.901 | 0.563 | 49.69% | 0.080 |
| **SIV** | retrovirus_int | 44.4% | 2.66% | **1.168** | 0.976 | 0.345 | 45.90% | 0.097 |
| **MoMLV** | retrovirus_int | 53.0% | 3.30% | **1.153** | 1.075 | 0.505 | 54.10% | 0.139 |
| **HPV16** | dna_onco_small | 36.5% | 1.83% | **0.542** | 0.862 | 0.426 | 45.77% | 0.043 |
| **HPV18** | dna_onco_small | 40.4% | 1.59% | **0.476** | 0.793 | 0.525 | 48.98% | 0.038 |
| **HPV6b** | dna_benign_ctrl | 40.9% | 1.67% | **0.516** | 0.846 | 0.486 | 51.56% | 0.039 |
| **HPV11** | dna_benign_ctrl | 41.0% | 1.69% | **0.522** | 0.881 | 0.466 | 64.67% | 0.036 |
| **MCPyV** | dna_onco_small | 40.7% | 4.64% | **1.326** | 0.986 | 0.157 | 49.16% | 0.137 |
| **SV40** | dna_onco_small | 40.8% | 3.76% | **1.001** | 0.890 | 0.124 | 49.64% | 0.085 |
| **JCV** | dna_benign_ctrl | 40.4% | 3.57% | **0.998** | 0.950 | 0.076 | 54.17% | 0.103 |
| **BKV** | dna_benign_ctrl | 39.5% | 3.50% | **0.956** | 0.932 | 0.060 | 47.50% | 0.109 |
| **HBV** | hepatitis_hcc | 48.5% | 4.76% | **1.209** | 0.984 | 0.536 | 50.59% | 0.109 |
| **HCV** | hepatitis_hcc | 58.2% | 3.19% | **1.149** | 1.042 | 0.731 | 46.02% | 0.058 |
| **EBV_HHV4** | herpes_large | 59.5% | 2.96% | **1.187** | 1.060 | 0.614 | 49.83% | 0.076 |
| **KSHV_HHV8** | herpes_large | 53.8% | 2.90% | **0.995** | 0.993 | 0.816 | 48.45% | 0.053 |
| **HSV1** | herpes_large | 68.3% | 1.56% | **0.920** | 0.990 | 1.009 | 50.43% | 0.029 |
| **HCMV_HHV5** | herpes_large | 57.5% | 2.68% | **1.049** | 1.071 | 1.193 | 49.33% | 0.053 |
| **AAV2** | vector_dna | 53.8% | 3.75% | **1.459** | 1.146 | 0.786 | 46.56% | 0.100 |
| **HAdV5** | vector_dna | 55.2% | 2.48% | **0.916** | 1.021 | 0.887 | 49.35% | 0.052 |
| **SARS_CoV_2** | rna_acute_ctrl | 38.0% | 3.65% | **0.998** | 1.018 | 0.408 | 49.04% | 0.105 |
| **SARS_CoV_1** | rna_acute_ctrl | 40.8% | 3.93% | **1.081** | 1.017 | 0.460 | 46.72% | 0.102 |
| **MERS_CoV** | rna_acute_ctrl | 41.2% | 3.90% | **1.003** | 0.998 | 0.555 | 48.22% | 0.081 |
| **Ebola_Zaire** | rna_acute_ctrl | 41.1% | 4.00% | **1.183** | 0.988 | 0.605 | 43.72% | 0.094 |
| **Measles** | rna_acute_ctrl | 47.4% | 3.86% | **1.314** | 1.006 | 0.493 | 45.81% | 0.101 |
| **Poliovirus** | rna_acute_ctrl | 46.3% | 3.18% | **1.055** | 0.967 | 0.538 | 51.91% | 0.095 |
| **Zika** | rna_acute_ctrl | 50.9% | 3.13% | **1.366** | 1.043 | 0.427 | 48.83% | 0.131 |
| **Dengue2** | rna_acute_ctrl | 45.8% | 2.99% | **1.277** | 0.980 | 0.411 | 42.32% | 0.146 |
| **Phage_Lambda** | phage_outgroup | 49.9% | 3.11% | **1.070** | 1.111 | 1.037 | 48.07% | 0.097 |
| **Phage_T4** | phage_outgroup | 35.3% | 4.48% | **1.120** | 1.036 | 0.942 | 46.87% | 0.121 |
| **Phage_PhiX174** | phage_outgroup | 44.8% | 3.63% | **0.977** | 1.068 | 0.991 | 45.63% | 0.112 |
> **Mann–Whitney TCW O/E Base (Integrating/Onco vs Acute RNA):** median `1.047` vs `1.132`, `r_rb = -0.361`, `p = 1.60e-01`.
> **Mann–Whitney TCW O/E Dinuc (Integrating/Onco vs Acute RNA):** median `0.985` vs `1.002`, `r_rb = -0.181`, `p = 4.95e-01`.
> **Mann–Whitney TCW O/E Base (Retroviruses+Papilloma/Polyoma vs Acute RNA):** median `0.999` vs `1.132`, `r_rb = -0.607`, `p = 1.97e-02`.

## 11. Theorems V11 & V11.2 — Beyond-Dinucleotide (Markov-1) & Beyond-Trinucleotide (Markov-2) Viral Mimicry Proof

| Virus | Class | GC % | cos Total (1024-D) | cos 3600-D | cos CpG-Free | cos 3600-D CpG-Free | cos Beyond-M1 (Dinuc) | cos Beyond-M2 (Trinuc) | cos(Self, M2) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **KSHV_HHV8** | herpes_large | 53.8% | +0.504 | +0.170 | +0.507 | +0.146 | +0.333 | **+0.232** | +0.700 |
| **EBV_HHV4** | herpes_large | 59.5% | +0.431 | +0.063 | +0.431 | +0.054 | +0.223 | **+0.165** | +0.696 |
| **HCMV_HHV5** | herpes_large | 57.5% | +0.081 | +0.058 | +0.151 | +0.029 | +0.217 | **+0.146** | +0.685 |
| **HAdV5** | vector_dna | 55.2% | +0.303 | +0.096 | +0.343 | +0.081 | +0.195 | **+0.140** | +0.559 |
| **HSV1** | herpes_large | 68.3% | +0.055 | +0.035 | +0.079 | +0.032 | +0.151 | **+0.127** | +0.540 |
| **Phage_Lambda** | phage_outgroup | 49.9% | +0.181 | +0.038 | +0.218 | +0.053 | +0.217 | **+0.105** | +0.697 |
| **HPV16** | dna_onco_small | 36.5% | +0.187 | +0.085 | +0.171 | +0.079 | +0.140 | **+0.092** | +0.443 |
| **AAV2** | vector_dna | 53.8% | +0.211 | +0.065 | +0.216 | +0.060 | +0.079 | **+0.089** | +0.289 |
| **HIV2_BEN** | retrovirus_int | 45.7% | +0.352 | +0.131 | +0.326 | +0.112 | +0.161 | **+0.088** | +0.409 |
| **HTLV2** | retrovirus_int | 53.8% | +0.238 | +0.046 | +0.215 | +0.035 | +0.134 | **+0.087** | +0.386 |
| **Poliovirus** | rna_acute_ctrl | 46.3% | +0.212 | +0.032 | +0.193 | +0.024 | +0.062 | **+0.077** | +0.359 |
| **BKV** | dna_benign_ctrl | 39.5% | +0.342 | +0.102 | +0.319 | +0.051 | +0.163 | **+0.072** | +0.475 |
| **HCV** | hepatitis_hcc | 58.2% | +0.253 | +0.042 | +0.245 | +0.050 | +0.094 | **+0.072** | +0.312 |
| **HTLV1** | retrovirus_int | 53.5% | +0.273 | +0.062 | +0.255 | +0.044 | +0.150 | **+0.069** | +0.442 |
| **SV40** | dna_onco_small | 40.8% | +0.347 | +0.130 | +0.335 | +0.100 | +0.160 | **+0.059** | +0.496 |
| **HPV11** | dna_benign_ctrl | 41.0% | +0.174 | +0.038 | +0.154 | +0.031 | +0.082 | **+0.055** | +0.491 |
| **HPV18** | dna_onco_small | 40.4% | +0.214 | +0.092 | +0.202 | +0.083 | +0.143 | **+0.055** | +0.525 |
| **JCV** | dna_benign_ctrl | 40.4% | +0.371 | +0.101 | +0.341 | +0.051 | +0.139 | **+0.052** | +0.522 |
| **HPV6b** | dna_benign_ctrl | 40.9% | +0.161 | +0.020 | +0.143 | +0.013 | +0.062 | **+0.046** | +0.487 |
| **Dengue2** | rna_acute_ctrl | 45.8% | +0.407 | +0.126 | +0.388 | +0.088 | +0.172 | **+0.036** | +0.590 |
| **HBV** | hepatitis_hcc | 48.5% | +0.188 | +0.023 | +0.180 | +0.010 | +0.049 | **+0.029** | +0.324 |
| **Zika** | rna_acute_ctrl | 50.9% | +0.406 | +0.096 | +0.377 | +0.051 | +0.116 | **+0.028** | +0.574 |
| **MoMLV** | retrovirus_int | 53.0% | +0.195 | -0.011 | +0.158 | -0.030 | -0.040 | **+0.025** | +0.437 |
| **MCPyV** | dna_onco_small | 40.7% | +0.292 | +0.067 | +0.258 | +0.023 | +0.098 | **+0.020** | +0.444 |
| **Phage_PhiX174** | phage_outgroup | 44.8% | +0.104 | +0.074 | +0.125 | +0.081 | +0.035 | **+0.018** | +0.342 |
| **Ebola_Zaire** | rna_acute_ctrl | 41.1% | +0.216 | +0.036 | +0.200 | +0.013 | -0.006 | **-0.007** | +0.422 |
| **SIV** | retrovirus_int | 44.4% | +0.272 | +0.059 | +0.254 | +0.048 | +0.086 | **-0.014** | +0.429 |
| **Measles** | rna_acute_ctrl | 47.4% | +0.254 | +0.017 | +0.213 | -0.013 | -0.039 | **-0.019** | +0.521 |
| **HIV1_HXB2** | retrovirus_int | 42.1% | +0.317 | +0.075 | +0.285 | +0.040 | +0.084 | **-0.020** | +0.461 |
| **Phage_T4** | phage_outgroup | 35.3% | +0.153 | -0.008 | +0.177 | -0.019 | +0.073 | **-0.029** | +0.596 |
| **SARS_CoV_1** | rna_acute_ctrl | 40.8% | +0.360 | +0.069 | +0.333 | +0.035 | +0.060 | **-0.029** | +0.623 |
| **MERS_CoV** | rna_acute_ctrl | 41.2% | +0.260 | +0.044 | +0.240 | +0.017 | +0.009 | **-0.048** | +0.505 |
| **SARS_CoV_2** | rna_acute_ctrl | 38.0% | +0.242 | -0.012 | +0.205 | -0.041 | -0.021 | **-0.079** | +0.592 |
> **Spearman Control (GC % vs Total 1024-D Alignment, n=33):** `rho_S = -0.005`, `p = 9.79e-01`.
> **Spearman Control (GC % vs Beyond-Markov-2 Alignment, n=33):** `rho_S = +0.559`, `p = 7.11e-04`.
> **Spearman Control (TCW O/E vs De Novo APOBEC3 χ_--^unb, n=33):** `rho_S = +0.626`, `p = 9.76e-05`.
> **Spearman Control (log10 Length vs Beyond-Markov-2 Alignment, n=33):** `rho_S = +0.103`, `p = 5.70e-01`.
> **Covariate-Adjusted Control (Beyond-M2 Excess Onco/Int vs Acute RNA (GC + log10L Residualized), n=18/8):** resid median `+0.0164` vs `-0.0434`, `r_rb = +0.597`, `p = 1.58e-02`.
> **Covariate-Adjusted Control (Beyond-M2 Excess Nuclear dsDNA vs Acute RNA (GC + log10L Residualized), n=14/8):** resid median `+0.0271` vs `-0.0570`, `r_rb = +0.893`, `p = 1.88e-04`.
> **Covariate-Adjusted Control (Beyond-M2 Gram-Schmidt Nuclear dsDNA vs Acute RNA (GC + log10L Residualized), n=14/8):** resid median `+0.0196` vs `-0.0275`, `r_rb = +0.679`, `p = 8.12e-03`.

## 12. Theorem V12 — Exact All-Orders Finite-Dose Mutagenesis Phase Transition & Double-Hit Kataegis Hessian

| Target | Signature | Hessian ‖H_--‖²×10⁴ | cos(J_--, H_--) | PR2-ok (ε=0) | PR2-ok (ε=0.10) | PR2-ok (ε=0.50) | PR2-ok (ε=1.00) | Nonlinear Dev (ε=0.25) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **TP53** | `SBS1_CpG_sym` | 0.424 | -0.749 | 80.07% | 80.60% | 81.53% | 80.41% | 1.0% |
| **TP53** | `SBS2_APOBEC_plus` | 0.388 | -0.461 | 80.07% | 80.46% | 78.97% | 71.97% | 0.5% |
| **TP53** | `SBS13_APOBEC_plus` | 0.386 | -0.481 | 80.07% | 80.26% | 77.30% | 67.80% | 0.3% |
| **PIK3CA** | `SBS1_CpG_sym` | 0.024 | -0.744 | 87.13% | 87.54% | 88.45% | 88.18% | 0.1% |
| **PIK3CA** | `SBS2_APOBEC_plus` | 0.102 | +0.103 | 87.13% | 87.40% | 79.46% | 65.99% | 1.6% |
| **PIK3CA** | `SBS13_APOBEC_plus` | 0.096 | -0.095 | 87.13% | 86.76% | 75.87% | 58.25% | 0.6% |
| **TERT** | `SBS1_CpG_sym` | 2.061 | -0.772 | 64.48% | 66.08% | 69.62% | 67.08% | 0.4% |
| **TERT** | `SBS2_APOBEC_plus` | 0.117 | -0.192 | 64.48% | 64.34% | 62.42% | 58.36% | 0.1% |
| **TERT** | `SBS13_APOBEC_plus` | 0.106 | +0.075 | 64.48% | 64.30% | 61.11% | 54.34% | 0.2% |
| **HRAS** | `SBS1_CpG_sym` | 4.141 | -0.744 | 58.80% | 60.66% | 64.91% | 62.74% | 0.8% |
| **HRAS** | `SBS2_APOBEC_plus` | 0.936 | -0.595 | 58.80% | 59.27% | 59.86% | 58.21% | 0.1% |
| **HRAS** | `SBS13_APOBEC_plus` | 0.791 | -0.425 | 58.80% | 59.15% | 59.32% | 57.05% | 0.3% |
| **MYC** | `SBS1_CpG_sym` | 1.273 | -0.753 | 64.54% | 65.90% | 69.06% | 67.35% | 0.2% |
| **MYC** | `SBS2_APOBEC_plus` | 0.178 | -0.157 | 64.54% | 65.57% | 67.20% | 63.79% | 0.1% |
| **MYC** | `SBS13_APOBEC_plus` | 0.175 | -0.298 | 64.54% | 64.83% | 63.05% | 56.46% | 0.2% |
| **BRCA1** | `SBS1_CpG_sym` | 0.037 | -0.768 | 96.60% | 96.77% | 97.05% | 96.79% | 0.1% |
| **BRCA1** | `SBS2_APOBEC_plus` | 0.173 | -0.178 | 96.60% | 96.44% | 90.28% | 76.66% | 2.9% |
| **BRCA1** | `SBS13_APOBEC_plus` | 0.159 | -0.232 | 96.60% | 96.36% | 88.64% | 70.58% | 2.6% |
| **FHIT** | `SBS1_CpG_sym` | 0.001 | -0.788 | 97.77% | 97.85% | 97.99% | 97.89% | 0.3% |
| **FHIT** | `SBS2_APOBEC_plus` | 0.141 | +0.073 | 97.77% | 96.82% | 78.32% | 60.08% | 2.4% |
| **FHIT** | `SBS13_APOBEC_plus` | 0.118 | -0.110 | 97.77% | 96.50% | 76.32% | 54.44% | 0.8% |
| **HPV16** | `SBS1_CpG_sym` | 0.469 | -0.728 | 59.00% | 59.36% | 60.78% | 61.76% | 0.4% |
| **HPV16** | `SBS2_APOBEC_plus` | 0.018 | -0.160 | 59.00% | 59.03% | 58.87% | 58.14% | 0.0% |
| **HPV16** | `SBS13_APOBEC_plus` | 0.019 | -0.070 | 59.00% | 59.19% | 59.24% | 57.89% | 0.0% |
| **HIV1_HXB2** | `SBS1_CpG_sym` | 0.272 | -0.760 | 52.67% | 52.96% | 53.62% | 53.27% | 0.1% |
| **HIV1_HXB2** | `SBS2_APOBEC_plus` | 0.069 | -0.245 | 52.67% | 52.92% | 53.79% | 54.16% | 0.1% |
| **HIV1_HXB2** | `SBS13_APOBEC_plus` | 0.061 | -0.283 | 52.67% | 53.01% | 53.95% | 53.58% | 0.2% |
| **MCPyV** | `SBS1_CpG_sym` | 0.305 | -0.754 | 57.11% | 57.30% | 57.85% | 58.04% | 0.3% |
| **MCPyV** | `SBS2_APOBEC_plus` | 0.611 | -0.251 | 57.11% | 57.94% | 60.20% | 59.18% | 0.1% |
| **MCPyV** | `SBS13_APOBEC_plus` | 0.453 | -0.363 | 57.11% | 57.53% | 57.58% | 54.53% | 0.1% |
| **HBV** | `SBS1_CpG_sym` | 5.341 | -0.760 | 50.83% | 51.34% | 52.99% | 52.93% | 0.5% |
| **HBV** | `SBS2_APOBEC_plus` | 0.492 | -0.221 | 50.83% | 50.64% | 49.67% | 48.97% | 0.0% |
| **HBV** | `SBS13_APOBEC_plus` | 0.405 | -0.318 | 50.83% | 50.87% | 51.08% | 50.93% | 0.1% |

## 13. Theorem V13 — Exact Klein-Orbit 6-mer Motif Attribution (Viral Mimicry & APOBEC3 Kataegis)

### 13A. KSHV_HHV8 vs. Human Euchromatin (Total 1024-D W_-+)

| Rank | 6-mer Motif | Mirror 6-mer | Orbit Size | Share of Total Inner Product % | CpG? |
|---:|---|---|---:|---:|:---:|
| 1 | `CCC-AGG` | `CCT-GGG` | 4 | +2.81% | No |
| 2 | `GAA-AAC` | `GTT-TTC` | 4 | +2.32% | No |
| 3 | `CCA-GGA` | `TCC-TGG` | 4 | +2.16% | No |
| 4 | `GCC-ACC` | `GGT-GGC` | 4 | +1.88% | No |
| 5 | `CAG-GAG` | `CTC-CTG` | 4 | +1.81% | No |
| 6 | `GGA-GGC` | `GCC-TCC` | 4 | +1.52% | No |
| 7 | `GCC-AGG` | `CCT-GGC` | 4 | +1.42% | No |
| 8 | `GGC-AGC` | `GCT-GCC` | 4 | +1.34% | No |

### 13B. KSHV_HHV8 vs. Human Euchromatin (Beyond-Markov-2 Residual W_-+^⊥M2)

| Rank | 6-mer Motif | Mirror 6-mer | Orbit Size | Share of Total Inner Product % | CpG? |
|---:|---|---|---:|---:|:---:|
| 1 | `CCC-AGG` | `CCT-GGG` | 4 | +4.69% | No |
| 2 | `AAA-CAA` | `TTG-TTT` | 4 | +3.41% | No |
| 3 | `AGC-TGG` | `CCA-GCT` | 4 | +3.31% | No |
| 4 | `TTT-ATT` | `AAT-AAA` | 4 | +2.97% | No |
| 5 | `CCT-CCC` | `GGG-AGG` | 4 | +2.82% | No |
| 6 | `CCA-GGA` | `TCC-TGG` | 4 | +2.69% | No |
| 7 | `GCC-AGG` | `CCT-GGC` | 4 | +2.47% | No |
| 8 | `CAG-GAG` | `CTC-CTG` | 4 | +2.31% | No |

### 13C. EBV_HHV4 vs. Human Euchromatin (Beyond-Markov-2 Residual W_-+^⊥M2)

| Rank | 6-mer Motif | Mirror 6-mer | Orbit Size | Share of Total Inner Product % | CpG? |
|---:|---|---|---:|---:|:---:|
| 1 | `GCC-TCC` | `GGA-GGC` | 4 | +6.13% | No |
| 2 | `CCA-GCC` | `GGC-TGG` | 4 | +5.43% | No |
| 3 | `GCC-ACC` | `GGT-GGC` | 4 | +5.31% | No |
| 4 | `CCA-GGC` | `GCC-TGG` | 4 | +5.00% | No |
| 5 | `AAA-GAA` | `TTC-TTT` | 4 | +4.74% | No |
| 6 | `CCC-GGG` | `CCC-GGG` | 2 | +4.62% | Yes |
| 7 | `GCC-GGG` | `CCC-GGC` | 4 | +3.61% | Yes |
| 8 | `CTC-CTG` | `CAG-GAG` | 4 | +3.56% | No |

### 13D. TP53 APOBEC3 SBS2(+) Chiral Kataegis Injection (J_-- Self-Energy)

| Rank | 6-mer Motif | Mirror 6-mer | Orbit Size | Share of Total Inner Product % | CpG? |
|---:|---|---|---:|---:|:---:|
| 1 | `CTA-TCT` | `AGA-TAG` | 4 | +5.29% | No |
| 2 | `TTG-AGA` | `TCT-CAA` | 4 | +5.15% | No |
| 3 | `CTC-CAT` | `ATG-GAG` | 4 | +4.48% | No |
| 4 | `TAC-CTC` | `GAG-GTA` | 4 | +3.39% | No |
| 5 | `AAA-CTC` | `GAG-TTT` | 4 | +3.17% | No |
| 6 | `CTC-CGT` | `ACG-GAG` | 4 | +2.26% | Yes |
| 7 | `TCA-ATC` | `GAT-TGA` | 4 | +2.08% | No |
| 8 | `TGA-GAC` | `GTC-TCA` | 4 | +1.99% | No |

### 13E. PIK3CA APOBEC3 SBS2(+) Chiral Kataegis Injection (J_-- Self-Energy)

| Rank | 6-mer Motif | Mirror 6-mer | Orbit Size | Share of Total Inner Product % | CpG? |
|---:|---|---|---:|---:|:---:|
| 1 | `CTT-TTT` | `AAA-AAG` | 4 | +7.73% | No |
| 2 | `TTT-TTG` | `CAA-AAA` | 4 | +4.17% | No |
| 3 | `CAT-TTT` | `AAA-ATG` | 4 | +3.44% | No |
| 4 | `AAA-TAA` | `TTA-TTT` | 4 | +3.07% | No |
| 5 | `ATT-TTA` | `TAA-AAT` | 4 | +2.63% | No |
| 6 | `TTT-TAG` | `CTA-AAA` | 4 | +2.54% | No |
| 7 | `TTT-TTC` | `GAA-AAA` | 4 | +2.27% | No |
| 8 | `ATT-TTT` | `AAA-AAT` | 4 | +1.74% | No |


## 14. Theorem V18 — Global BH-FDR Ledger, ZAP 3-Covariate Control & T2T Rank-1 Diagnostics

> **T2T Rank-1 Factorization (924 pairs):** `cos(v, l) ~ cos(v, H)*cos(l, H)` explains **`R² = 0.7743`** (slope `0.956`). Raw top winners: `{'KSHV_HHV8': 25, 'EBV_HHV4': 3}`; Human-residual `W^⊥H` top winners: `{'HCMV_HHV5': 7, 'EBV_HHV4': 5, 'HTLV2': 3, 'HPV18': 3, 'MCPyV': 2, 'HSV1': 2, 'Phage_Lambda': 1, 'HPV16': 1, 'SARS_CoV_1': 1, 'HBV': 1, 'JCV': 1, 'Dengue2': 1}`.

> **Global Benjamini–Hochberg FDR (`m = 50` tests; `2` primary, `48` exploratory; `5` survive `q_BH ≤ 0.05`):**

| Rank | Tier | Section | Test | Two-Sided p | p_min | q_BH |
|---:|:---:|:---:|---|---:|---:|---:|
| 1 | `primary` | V11 | `tcw_oe_vs_apobec_denovo` | 9.76e-05 | n/a | 0.0031 |
| 2 | `primary` | V11 | `beyond_m2_nucDNA_vs_cytoRNA_adj_gc_L` | 1.88e-04 | 6.25e-06 | 0.0031 |
| 3 | `exploratory` | V11 | `beyond_m2_nucDNA_vs_cytoRNA_adj_gc_L_cpg` | 1.88e-04 | 6.25e-06 | 0.0031 |
| 4 | `exploratory` | V11 | `gc_vs_cos_beyond_m2` | 7.11e-04 | n/a | 0.0089 |
| 5 | `exploratory` | V9 | `Onco/Integrating Viruses vs Acute RNA Controls | cos_Wmp_beyond_M2_vs_human` | 3.35e-03 | 1.28e-06 | 0.0335 |
| 6 | `exploratory` | V11 | `beyond_m2_gs_nucDNA_vs_cytoRNA_adj_gc_L_cpg` | 6.37e-03 | 6.25e-06 | 0.0508 |
| 7 | `exploratory` | V9 | `Onco/Integrating Viruses vs Acute RNA Controls | pure_nonenv_3600_unb` | 7.63e-03 | 1.28e-06 | 0.0508 |
| 8 | `exploratory` | V11 | `beyond_m2_gs_nucDNA_vs_cytoRNA_adj_gc_L` | 8.12e-03 | 6.25e-06 | 0.0508 |
| 9 | `exploratory` | V11 | `beyond_m2_int_vs_rna_adj_gc_L_cpg` | 9.23e-03 | 1.28e-06 | 0.0512 |
| 10 | `exploratory` | V9 | `Human Viruses (n=30) vs Bacteriophages (n=3) | cos_Wmp_unb_vs_human_grand` | 1.02e-02 | 4.45e-04 | 0.0512 |
| 11 | `exploratory` | V9 | `Human Viruses (n=30) vs Bacteriophages (n=3) | cos_Wmp_vs_human_grand` | 1.50e-02 | 3.67e-04 | 0.0608 |
| 12 | `exploratory` | V9 | `Onco/Integrating Viruses vs Acute RNA Controls | cos_Wmp_beyond_M1_vs_human` | 1.58e-02 | 1.28e-06 | 0.0608 |
| 13 | `exploratory` | V11 | `beyond_m2_int_vs_rna_adj_gc_L` | 1.58e-02 | 1.28e-06 | 0.0608 |
| 14 | `exploratory` | V10 | `mw_tcw_oe_base_retro_pap_vs_rna` | 1.97e-02 | 6.25e-06 | 0.0702 |
| 15 | `exploratory` | V9 | `Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4) | pal_d6_Wmp_share` | 3.51e-02 | 9.77e-05 | 0.1169 |
