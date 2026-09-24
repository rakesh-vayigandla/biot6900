# BIOT 6900 coursework
Rakesh Vayigandla

Environment: Python 3.11 (biot6900 conda env).

## Module 1
`module1_setup.ipynb`. All of Parts A to C ran cleanly. Part D left commented out.

## Module 2: Multi-omics target identification
- `BIOT6900_Module2_Starter.ipynb`: main notebook (Parts 1 to 3). Runs top to bottom with Restart & Run All.
- `BIOT6900_Module2_Part2_LoadRealCPTAC copy.ipynb`: turns the LinkedOmics CPTAC breast cancer downloads into the Part 1 input files.
- `Module2_Report.md`: written report for Part 3.
- `targets_t2d.csv`: ranked target list (hand-off to Week 3).
- `prepare_t2d_data.py`: downloads and builds the Part 3 tables.

Data, one folder per part:
- `data/part1_cptac_brca/`: processed CPTAC breast cancer matrices used by Part 1.
- `data/part2_cptac_brca_raw/`: raw LinkedOmics downloads used in Part 2.
- `data/part3_t2d/`: type 2 diabetes tables used by Part 3 (islet RNA from GEO GSE76894, UK Biobank plasma proteomics from Gadd et al. 2024, Open Targets GWAS evidence). The raw downloads (~100 MB) are not committed; run `python prepare_t2d_data.py` from the repo root to rebuild them (needs `pip install openpyxl`).

Notes / what didn't work:
- For Part 3 I used type 2 diabetes instead of Alzheimer's (disease choice was open). I first tried rheumatoid arthritis, but the only proteomics table I found had 168 proteins, which left too few genes after the join.
- The T2D join is limited to 1,024 genes by the plasma protein panel (~1,460 proteins), so some well-known T2D genes (PPARG, TCF7L2, SLC30A8, KCNJ11) drop out. This is discussed in the report.
