# BIOT 6900 coursework
Rakesh Vayigandla

Environment: Python 3.11 (biot6900 conda env).

Each module has its own folder. Notebooks use paths relative to their own folder, so open and run them from there.

```
module1/   setup notebook
module2/   multi-omics target identification (T2D)
module3/   target triage with protein language models (T2D)
module4/   structure prediction and pocket assessment (in progress)
```

## Module 1
`module1/module1_setup.ipynb`. All of Parts A to C ran cleanly. Part D left commented out.

## Module 2: Multi-omics target identification
- `module2/BIOT6900_Module2_Starter.ipynb`: main notebook (Parts 1 to 3). Runs top to bottom with Restart & Run All.
- `module2/BIOT6900_Module2_Part2_LoadRealCPTAC copy.ipynb`: turns the LinkedOmics CPTAC breast cancer downloads into the Part 1 input files.
- `module2/Module2_Report.md`: written report for Part 3.
- `module2/targets_t2d.csv`: ranked target list (hand-off to Week 3).
- `module2/prepare_t2d_data.py`: downloads and builds the Part 3 tables.

Data, one folder per part, in `module2/data/`:
- `part1_cptac_brca/`: processed CPTAC breast cancer matrices used by Part 1.
- `part2_cptac_brca_raw/`: raw LinkedOmics downloads used in Part 2.
- `part3_t2d/`: type 2 diabetes tables used by Part 3 (islet RNA from GEO GSE76894, UK Biobank plasma proteomics from Gadd et al. 2024, Open Targets GWAS evidence). The raw downloads (~100 MB) are not committed; run `python prepare_t2d_data.py` from `module2/` to rebuild them (needs `pip install openpyxl`).

Notes / what didn't work:
- For Part 3 I used type 2 diabetes instead of Alzheimer's (disease choice was open). I first tried rheumatoid arthritis, but the only proteomics table I found had 168 proteins, which left too few genes after the join.
- The T2D join is limited to 1,024 genes by the plasma protein panel (~1,460 proteins), so some well-known T2D genes (PPARG, TCF7L2, SLC30A8, KCNJ11) drop out. This is discussed in the report.

## Module 3: Target triage with protein language models
Disease: type 2 diabetes (continuing from Module 2).
- `module3/BIOT6900_Module3_Assignment.ipynb`: assignment notebook (Steps 1 to 7). Runs top to bottom with Restart & Run All.
- `module3/BIOT6900_Module3_Part2_RealData.ipynb`: lab Part 2 notebook.
- `module3/Module3_Report.pdf`: written report (`Module3_Report.md` is the editable source).
- `module3/targets_t2d_w3.csv`: re-ranked targets with small-molecule and antibody tractability priors (hand-off to Week 4).
- `module3/figures/module3_modality_map.png`: modality map used in the report.
- `module3/data/`: the Module 3 data pack (ESM C embeddings, training panel, anchor variant scores).

Notes:
- Input is my Module 2 file, read from `../module2/targets_t2d.csv`. It already uses `gene` and `score` columns, so no renaming was needed. The only change to the provided settings is that `TARGETS_FILE` path.
- All 15 top targets mapped to the proteome, so there were no unmapped genes. CRH was in the training panel and was removed before fitting (leakage guard).
- I added one extra cell after Step 4 that validates the small-molecule and antibody priors separately with the family split.
