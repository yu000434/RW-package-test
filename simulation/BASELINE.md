# Baseline

Frozen 2026-09-03 before simulation cleanup.

Method: `MICE_soft_PMM_donor_source_research`.

| Item | Value |
|---|---|
| Source files | 35 |
| Source tree SHA-256 | `81a89b5b850984751ef2bb48c0c539e42b04cd92a1770811128e52793a85e031` |
| Result files | 1,793 |
| Result bytes | 40,297,258 |
| Result tree SHA-256 | `876fbd53ef386597618b407b3ea83d42316dda2141cfefd622ee025e35569c44` |
| GS reference | `results/gs5_5/raw/gs_task0001.csv`, seed 8500000 |
| GS reference SHA-256 | `ae483042ca7be3a8c98d59b57b4d60643a186d99443a2da45a3bc7b94c9aa4a6` |
| RW reference | `results/rw150_20_5/raw/rw_task0001.csv`, S1 seed 12000000 |
| RW reference SHA-256 | `a1498799840efb12e6be434d3c7106958cac1d58ff6881e26838fd53d664c2a9` |

Run `Rscript tests/test_mice_pmm_parity.R` and `Rscript tests/smoke.R` after each cleanup change. The smoke comparison permits numerical differences below `1e-10`, which covers CSV round-trip differences without accepting a material change.
