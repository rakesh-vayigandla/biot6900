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
module6/   polygenic risk score portability (Alzheimer's PGS000334)
module7/   knowledge graph in Neo4j (T2D)
module8/   single-cell and spatial transcriptomics (lab notebooks)
module9/   clinical NLP: signal vs leakage (DKA)
portfolio1/ Portfolio Deliverable 1: T2D multi-omics target discovery package (Labs 2-4)
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

## Module 6: PRS portability across ancestries
Score: PGS Catalog **PGS000334** (22-SNP Alzheimer's score, provided with the course data), applied to 2,336 people from 1000 Genomes phase 3.
- `module6/BIOT6900_Module6_Part3_Graded_Starter.ipynb`: graded notebook. `USE_DEMO = False`; Section 2 shows `source: REAL` with the file checksums. Runs top to bottom with Restart & Run All.
- `module6/Module6_Report.pdf`: written report, including the proof-of-data block (`Module6_Report.md` is the editable source).
- `module6/figures/`: PRS distribution by ancestry and calibration by ancestry.
- `module6/BIOT6900_Module6_Toy_PRS.ipynb`, `BIOT6900_Module6_Part2_RealPRS.ipynb`: lab Parts 1 and 2.

Notes:
- The phenotype in this cohort is simulated (1000 Genomes has no disease labels), so AUC and calibration are illustrative.
- My earlier modules used type 2 diabetes, but the genotype file contains only the 22 Alzheimer's score SNPs and the phenotype was simulated from that score, so I used the provided Alzheimer's score rather than a T2D one.
- `module6/data/` (the Canvas data pack) is not committed.
- **Extra, not graded (`module6/t2d_extra/`):** my Part I disease is T2D, so I also ran the distribution analysis on two real PGS Catalog T2D scores (PGS003733 EUR-derived, PGS003735 AFR-derived), extracting the 287 genotypes they need straight from the 1000 Genomes VCFs for the same 2,336 people. There is no T2D phenotype for 1000 Genomes, so it reports the PRS distribution only, with no AUC. See the appendix of the report, `t2d_extra/build_t2d_prs.py` and `t2d_extra/t2d_prs_by_ancestry.csv`. The raw VCF extract (`t2d_extra/raw/`) is not committed.

## Portfolio Deliverable 1: Multi-omics target discovery package
Disease: type 2 diabetes (the same disease as Modules 2-4). A 6-page synthesis of Labs 2-4, ending in a ranked target shortlist (LPL, CSF1, IL17RB, CTSZ).
- `portfolio1/Portfolio1_T2D_Target_Package.pdf`: the package (`.md` is the editable source).
- `portfolio1/figures/`: Venn diagram of layer coverage, multi-omics integration plot, and AlphaFold structure figure with pockets. These were re-plotted from the committed Lab 2 tables and Lab 4 AlphaFold/fpocket outputs; no new analysis.
- Underlying work: `module2/` (ranked table `targets_t2d.csv`), `module3/`, `module4/`.

## Module 7: Knowledge graph (Neo4j)
Disease: type 2 diabetes. A **77-node, 84-relationship** graph in Neo4j AuraDB Free: 25 genes, 19 GWAS lead variants, 17 Reactome pathways, 9 approved drugs, 5 phenotypes, and 2 diseases (T2D plus non-alcoholic fatty liver disease as the related disease that shares phenotypes).
- `module7/BIOT6900_Module7_Starter.ipynb`: assignment notebook (Parts A to E). Runs top to bottom with Restart & Run All.
- `module7/Module7_Report.pdf`: written report (`Module7_Report.md` is the editable source).
- `module7/my_knowledge_graph.png`: the graph coloured by node kind, with a legend.
- `module7/nodes.csv`, `module7/edges.csv`: the node and edge tables.
- `module7/biomarker_ranking.csv`: PageRank fused with the Module 2/3 gene evidence.
- `module7/data/targets_t2d.csv`: the gene evidence table built from my Module 2/3 results plus Open Targets.

Notes / what didn't work:
- **Credentials are not committed.** The notebook reads them from `module7/neo4j_credentials.txt` (the Aura "Download credentials" file), which is git-ignored.
- The AuraDB Free instance had paused itself after a few days, so the first connection attempt failed to resolve the host; resuming it in the Aura console fixed it.
- Data sources: Open Targets (GWAS credible sets, HPO phenotypes, Reactome pathways, ChEMBL drug mechanisms), queried October 2026. Every edge carries a `source` and `confidence` property.
- I added one extra cell after the provided visualization cell that redraws the same graph with a node-kind legend and degree-scaled nodes, and one after the PageRank cell that fuses PageRank with the gene evidence into a single ranked table.

## Module 9: Clinical NLP, signal vs leakage
Disease: **diabetic ketoacidosis (DKA)**, an acute complication of the type 2 diabetes thread from Modules 2-7. All patients and notes are synthetic, generated inside the notebook from my configuration (`SEED = 2020300`, my student ID).
- `module9/BIOT6900_Module9_Notebook3_Graded.ipynb`: the graded notebook, with the D1-D4 answers written in its markdown cells. Runs top to bottom with Restart & Run All.

Configuration and findings:
- **Signals (causes, documented at admission):** severe acidosis, altered mental status, precipitating infection.
- **Consequences (interventions on deteriorating patients, my planted leakage):** pressors, intubation, comfort care.
- **The two injected wildcards:** `lactate_not_clearing` is genuine early signal (admission note only, 1.2x enriched in those who died); `renal_replacement` is leakage (final note only, 5.7x enriched, and the single most important feature in the model).
- **AUROC:** 0.929 with everything, 0.858 with my declared leakage removed, **0.739 with all leakage removed**. The one wildcard I had not declared cost more (-0.119) than all three declared leaks combined (-0.071).
- Assertion check: naive keyword counting found 206 patients with severe acidosis; assertion-aware extraction found 132, so negation handling removed 74 false positives.

Notes:
- I added one cell before Part D (clearly labelled) that reports which note each concept appears in and how enriched it is by outcome. That timing evidence is what the D2 answer argues from, rather than relying on the importance plot alone.
- Nothing failed to run.
