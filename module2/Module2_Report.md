# Multi-Omics Target Identification for Type 2 Diabetes

*BIOT 6900 · Module 2, Part 3 · Rakesh Vayigandla · September 2026*

## 1. Disease and data choice

### Why type 2 diabetes

I chose type 2 diabetes (T2D) for three reasons. First, it is one of the most common chronic diseases worldwide, and there is still a need for targets that act on its causes rather than only lowering blood glucose. Second, its genetics are unusually well mapped: large GWAS have identified hundreds of risk loci, so the genomic layer is rich. Third, T2D already has well-characterised drug targets and textbook genes (e.g. PPARG, KCNJ11, TCF7L2, SLC30A8), which makes it possible to check whether the pipeline recovers known biology. I first tried rheumatoid arthritis, but the only proteomics table available covered 168 proteins, which left too few genes after integration. For T2D, large open datasets exist for all three layers.

### Datasets

| Layer | Dataset | Source | Access | Size |
|:---|:---|:---|:---|:---|
| Transcriptomic | Pancreatic islets from 103 organ donors, 19 T2D vs 84 non-diabetic (Solimena et al. 2018) | NCBI GEO, GSE76894 | Open | 14,077 genes |
| Proteomic | 1,458 plasma proteins (Olink) vs later T2D diagnosis in 47,600 participants, 2,822 of whom developed T2D (Gadd et al. 2024) | UK Biobank Pharma Proteomics Project, published summary statistics (Supplementary Table 4) | Summary statistics open; individual-level data gated (UK Biobank application required) | 1,458 proteins |
| Genomic | GWAS credible-set evidence score (0 to 1) for T2D | Open Targets Platform (MONDO_0005148) | Open | 3,319 genes |

For the transcriptomic layer, each gene has a log2 fold change (T2D minus non-diabetic) and a Welch t-test p-value. For the proteomic layer, each protein has a hazard ratio (HR) per standard deviation, used as log2 HR. For the genomic layer, I used only the GWAS evidence and excluded literature and clinical-drug evidence, because approved drug targets would otherwise score themselves.

### Matched or unmatched?

The three layers are **unmatched**: they come from three unrelated groups of people (European organ donors, UK Biobank volunteers, and GWAS consortia). They also come from different tissues: islets, blood plasma, and germline DNA. No sample was measured at more than one layer, so I integrated at the **gene level**.

### Integration method

RNA and protein were inner-joined on gene symbol, because a gene needs both measurements, which left 1,024 genes. The GWAS table lists only genes with evidence, so it was left-joined, and genes missing from it were scored 0; 180 of the 1,024 genes have GWAS evidence. Concordance was defined as sign agreement: a gene is concordant if it is higher in T2D islets *and* higher plasma levels raise T2D risk, or both go the other way. Each layer's magnitude was converted to a rank percentile and averaged with equal weights.

## 2. Top targets and known biology

The join is limited by the protein panel: 1,458 plasma proteins shrink to 1,024 genes after intersecting with islet RNA. Only 557 of the 1,024 genes (54%) are sign-concordant, barely above the 50% expected by chance, and the RNA and protein effects are almost uncorrelated (Spearman ρ = 0.08). The layers therefore carry largely independent information. Table 1 compares the top 10 hits with the Open Targets overall T2D association, which combines all evidence types, including literature, across 9,907 T2D-linked genes.

**Table 1. Top 10 ranked T2D targets.** HR < 1 means a higher plasma protein level goes with lower T2D risk. *OT rank* is the gene's position among the 9,907 genes Open Targets associates with T2D (1 = strongest).

| Gene | Islet RNA log2FC | Plasma HR per SD | GWAS score | Concordant | Score | OT rank |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|
| *LPL* | 0.33 | 0.56 | 0.92 | no | 0.927 | 84 |
| *BAIAP2* | -0.41 | 1.78 | 0.35 | no | 0.925 | 1,446 |
| *CRH* | -0.88 | 0.70 | 0.30 | yes | 0.912 | 1,610 |
| *CSF1* | 0.34 | 1.58 | 0.59 | yes | 0.900 | 546 |
| *IL17RB* | -0.50 | 0.67 | 0.17 | yes | 0.899 | 2,105 |
| *IFNLR1* | 0.39 | 1.43 | 0.62 | yes | 0.895 | 491 |
| *ROBO1* | -0.33 | 1.60 | 0.42 | no | 0.891 | 1,189 |
| *PCDH17* | -0.42 | 1.39 | 0.51 | no | 0.888 | 830 |
| *NPY* | -0.40 | 1.47 | 0.36 | no | 0.885 | 1,260 |
| *SLIT2* | 0.47 | 1.32 | 0.66 | yes | 0.883 | 398 |

### Known genes recovered

**LPL** (lipoprotein lipase) ranks first. It is also the best-established T2D gene in the list: #84 in Open Targets, with strong literature support (literature score 0.87). LPL breaks down triglycerides carried in circulating lipoproteins, and higher plasma LPL is associated with lower T2D risk. **NPY** (Open Targets literature score 0.89) and **DLK1** (rank 18 in my list, literature score 0.66) also have substantial published links to islet biology. Both are lower in T2D islets, consistent with β-cells losing their secretory identity.

### Known genes missing, and why

Several textbook T2D genes do not appear at all, and in each case one layer is missing. PPARG and KCNJ11 appear only in the GWAS data. TCF7L2 and SLC30A8 have islet RNA and GWAS evidence but are not on the plasma protein panel. LEP and IGFBP1 are measured in plasma but absent from the islet RNA, which makes sense because they are made mainly by fat and liver, not islets. The strongest plasma marker overall, ACY1, is lost for the same reason. The absence of these genes reflects the data, not the biology.

### Non-obvious hits worth a second look

Most top hits are supported by GWAS evidence in Open Targets but have almost no published work linking them to T2D. The most interesting is the **SLIT-ROBO pathway**: SLIT2 (rank 10), ROBO1 (rank 7) and ROBO2 (rank 21) all rank highly. The same pathway appearing three times independently is more convincing than a single gene. SLIT2 has solid genetic support (GWAS score 0.66) but a literature score of only 0.05, so it is under-studied in T2D. **CSF1** and **IFNLR1** are concordant and GWAS-supported, pointing to an immune or inflammatory component.

### Effect of weighting

All 15 top hits have GWAS evidence. This is partly an artefact of the scoring: 844 of the 1,024 genes tie at a GWAS score of 0, so they all receive the same low percentile (≈0.41), while any gene with GWAS evidence jumps above ≈0.82. Up-weighting GWAS to 0.5 returns the same ten genes in a different order, so the ranking is robust to the weights. The on/off behaviour of GWAS would be better fixed with a continuous gene-level GWAS statistic.

## 3. A discordant gene: LPL

LPL is discordant: its RNA is slightly *higher* in T2D islets (log2FC 0.33, p = 0.04), yet a higher plasma LPL level predicts *lower* T2D risk (HR 0.56 per SD). This is expected rather than a data error, for three reasons.

1. **Different tissues.** Plasma LPL is released mainly from the capillaries of fat and muscle, not from islets, so islet RNA says little about how much LPL is in the blood.
2. **Different timing.** Plasma was sampled from healthy people years before diagnosis, whereas islet RNA comes from donors who already had T2D. In those islets, higher LPL may be a local response to excess lipids (lipotoxicity) rather than a cause of disease.
3. **Regulation after transcription.** How much LPL circulates depends on its release, binding and clearance, which are controlled at the protein level.

The RNA and protein layers here therefore answer different questions: what changes inside diseased islets, and what in the blood predicts future disease.

## 4. What these data let me claim, and not claim

**I can claim:**

- These genes show evidence from several independent sources. LPL, for example, is supported by genetics, a pre-diagnostic plasma signal, and islet expression, which makes it a stronger candidate than a gene seen in one layer.
- Plasma levels of proteins such as LPL are *associated* with future T2D risk, because the proteins were measured before diagnosis.
- Genes with GWAS evidence are more likely to be causal, because germline variants are present before disease onset.
- The ranked list is a reasonable prioritisation to take forward to Week 3.

**I cannot claim:**

- **Anything about individual patients.** No person was measured at more than one layer, so I cannot say that someone with high islet LPL RNA also has high plasma LPL. Concordance can only be judged by the *direction* of change, not measured as a correlation.
- **That RNA and protein are coupled for a given gene.** In Part 1, the matched CPTAC breast-cancer data measured the same 122 tumours at every layer. There, coupling could be measured directly as a per-gene correlation across patients (median ρ = 0.48, range -0.23 to 0.92, 6,306 genes). With unmatched data this is impossible.
- **Causality from the protein or RNA layers.** Plasma associations may reflect obesity or insulin resistance, and islet changes in donors with established T2D may be a *consequence* of high blood glucose.
- **That missing genes are unimportant.** PPARG, TCF7L2 and others are absent because of which proteins and genes each dataset measured, not because of their biology.

Further limitations are the small T2D islet group (19 donors), the use of organ donors, a plasma panel covering only a small fraction of the genome, GWAS scores that depend on how variants are assigned to genes, and a simple t-test with no adjustment for age, sex or BMI.

## 5. Conclusion

Integrating islet transcriptomics, plasma proteomics and GWAS nominated LPL as the top T2D target and highlighted the SLIT-ROBO pathway as an under-studied candidate. The weak RNA-protein agreement and the many known genes lost at the join show that, with unmatched data, the choice of cohorts and tissues shapes the result as much as the scoring method does.

## References

1. Solimena M, Schulte AM, Marselli L, et al. Systems biology of the IMIDIA biobank from organ donors and pancreatectomised patients defines a novel transcriptomic signature of islets from individuals with type 2 diabetes. *Diabetologia*. 2018;61(3):641-657. doi:10.1007/s00125-017-4500-3. Data: NCBI GEO accession GSE76894.
2. Gadd DA, Hillary RF, Kuncheva Z, et al. Blood protein assessment of leading incident diseases and mortality in the UK Biobank. *Nature Aging*. 2024;4(7):939-948. doi:10.1038/s43587-024-00655-7.
3. Buniello A, Suveges D, Cruz-Castillo C, et al. Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Research*. 2025;53(D1):D1467-D1475. doi:10.1093/nar/gkae1128. Data: type 2 diabetes mellitus (MONDO_0005148), accessed 24 September 2026.
