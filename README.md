# BIOT 6900 coursework
Rakesh Vayigandla

Environment: Python 3.11 (biot6900 conda env).

Each module has its own folder. Notebooks use paths relative to their own folder, so open and run them from there.

```
module1/   setup notebook
module2/   multi-omics target identification (T2D)
module3/   target triage with protein language models (T2D)
module4/   structure-based druggability (T2D)
module5/   variant interpretation (CYP2C19)
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

## Module 4: Structure-based druggability
Disease: type 2 diabetes. Targets: **LPL** (P06858, Module 2 rank 1) and **CTSZ** (Q9UBR2, Module 2 rank 15), the two targets the Module 3 report sent forward to Week 4.
- `module4/BIOT6900_Module4_Assignment.ipynb`: assignment notebook (Steps 1 to 7). Runs top to bottom with Restart & Run All.
- `module4/BIOT6900_Module4_Part1_ToyModel.ipynb`, `BIOT6900_Module4_Part2_RealData_1.ipynb`, `BIOT6900_Module4_Demo_DrugPocket.ipynb`: lab notebooks (Step 8 p53 numbers come from Part 2).
- `module4/Module4_Report.pdf`: written report (`Module4_Report.md` is the editable source).
- `module4/structure_assessment_t2d.csv`: summary table with the druggability calls.
- `module4/figures/`: plots from the notebook used in the report.

Notes:
- Experimental structures: LPL 6OB0 chain A (LPL-GPIHBP1 complex with a bound small molecule); CTSZ 1DEU chain A (procathepsin X, no small-molecule ligand, so the ligand-site check applies to LPL only).
- `module4/structures/` (AlphaFold models, PDB files, fpocket output) is not committed; the notebook downloads everything again.
- fpocket was installed into the biot6900 env with `conda install -c conda-forge fpocket`. If the notebook says fpocket is not found, make sure the kernel is the biot6900 env (its `bin/` must be on PATH).

## Module 5: Variant interpretation
Gene: **CYP2C19**. I was not assigned a gene on the allocation list, so I chose CYP2C19 myself.
- `module5/BIOT6900_Module5_Part3_Assignment.ipynb`: assignment notebook (Q1 to Q6). Runs top to bottom with Restart & Run All.
- `module5/Module5_Report.pdf`: written report (`Module5_Report.md` is the editable source).
- `module5/pgx_classified_CYP2C19.csv`: the 33 VUS with their CADD band, evidence and provisional call.

Notes / what didn't work:
- **AlphaGenome Atlas: not run.** I did not have an API key, so `SKIP_ATLAS = True`; the Atlas part of the report is answered from reasoning.
- **Calibration:** CYP2C19 has no P/LP variants in ClinVar, so CADD cannot be calibrated on it alone. I pooled it with DPYD (fetched and cached the same way); 130 of the 135 labelled variants in the pooled calibration come from DPYD, which the report discusses.
- `module5/data/` (MyVariant.info cache) is not committed; the notebook re-downloads it.
- The template's text cells had cell type `md`, which Jupyter rejects as invalid, so I changed them to `markdown` to make the notebook run. No text was changed.
