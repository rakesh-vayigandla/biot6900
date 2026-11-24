# Genomics and Precision Medicine in Type 2 Diabetes

## From variant burden to cell-type resolution on a multi-omics target set

Rakesh Vayigandla · BIOT 6900 · Portfolio Deliverable 2 · November 2026

---

## 1. Disease, targets, and the precision-medicine question

**Disease: type 2 diabetes (T2D).** This case study continues the target set I nominated in Part I, where integrating islet transcriptomics, plasma proteomics and GWAS evidence produced a ranked shortlist headed by **LPL, CSF1, IL17RB and CTSZ**. Part II asks a different question of those same targets. Part I asked *which genes matter*. Part II asks *which patients, which drug, and which cells*.

A precision-medicine decision in T2D is not a diagnosis: glucose and HbA1c diagnose it cheaply and reliably, so no genomic test competes with that. Genomics earns its place in three narrower decisions. **Who to test and when**, stratifying by inherited risk before onset so prevention can be targeted. **Which drug, at what dose**, since T2D patients carry a heavy co-medication burden running through pharmacogenes, and the actionable variant in a T2D clinic is frequently not in a diabetes gene at all. **Which cells to act on**, since T2D involves at least three tissues (beta cell, liver, adipose) and a target that looks attractive in bulk data may sit in a cell type no drug can usefully reach. This report takes those questions in turn across Labs 5 to 8, then asks whether the evidence converges.

**Two divergences, stated at the outset.** The evidence base assembled here is not uniformly on the T2D thread, because two labs were constrained by the data provided. **Lab 5 analysed CYP2C19**, since I was not assigned a gene from the allocation list. **Lab 8 analysed breast cancer tissue**, since the CELLxGENE Census query and the Visium section were fixed by the assignment. I have not disguised either. Section 2 argues that CYP2C19 is a legitimate pharmacogenomic axis *for T2D patients*, which is a claim about co-medication rather than about diabetes genetics. Section 5 presents Lab 8 as a method whose conclusion transfers while its tissue does not, and states plainly what it cannot support.

![**Figure 1.** Integrated evidence across the Part I multi-omics screen and the Lab 7 knowledge graph. **(A)** The fused evidence score combines network centrality, GWAS strength, druggability prior and plasma protein effect. The three highest-ranked genes are all existing drug targets, which is a circularity discussed in Section 4, not a discovery. **(B)** Restricting to genes with no approved T2D drug *and* a measured plasma protein effect isolates the genuine candidates. TNFRSF11A and LPL carry both the strongest genetics and the largest protein effects. Re-plotted from `module7/biomarker_ranking.csv`.](figures/fig1_integrated_evidence.png)

---

## 2. Variant interpretation and pharmacogenomics

### 2.1 Why CYP2C19 is the right pharmacogene for a T2D patient

The clinically actionable pharmacogenomics in T2D sits mostly in the comorbidities, not in the glucose pathway. Cardiovascular disease is the leading cause of death in T2D, and a large fraction of T2D patients are therefore prescribed **clopidogrel** after an event or a stent. Clopidogrel is a prodrug that requires CYP2C19 to be activated, and CPIC assigns a Level A recommendation: poor and intermediate metabolisers should receive an alternative agent such as prasugrel or ticagrelor rather than standard clopidogrel (Lee et al., 2022). CYP2C19 also clears proton-pump inhibitors, voriconazole and several antidepressants, all common in this population.

So the pharmacogenomic question I can actually answer for a T2D patient is: *before this patient starts clopidogrel, do they activate it?* That is a different question from the Part I target question, and I treat it as such throughout.

### 2.2 Allele frequency thresholds invert in a pharmacogene

The fetch returned 80 ClinVar variants in CYP2C19: 38 drug-response or other assertions, 33 variants of uncertain significance (VUS), 9 benign or likely benign, and **zero pathogenic or likely pathogenic**. Both anchors were recovered: CYP2C19\*2 (rs4244285) and \*17 (rs12248560).

The absence of pathogenic variants is biology, not a data failure. CYP2C19 is not a disease gene: losing its function is harmless until the patient takes a drug that depends on it, so no purifying selection keeps loss-of-function alleles rare. Here \*2 sits at **16.8%** allele frequency and \*17 at **20.0%**, which inverts the usual rare-disease filter. I set **MAX_AF = 0.5**, removing only 2 of 80 variants, both ClinVar benign and both above 90% where the reference base is actually the minor allele. The conventional rare-disease threshold of 0.001 would instead have **discarded 19 variants, including \*2, \*17, \*3, \*4 and \*8**: very nearly every allele CPIC uses to guide prescribing. A filter imported from Mendelian genetics deletes the clinical content of this gene.

I selected **CADD** (Rentzsch et al., 2019), which scores every variant class (85% of these variants), over REVEL and AlphaMissense, which score only missense variants (67.5%). Most CYP2C19 function alleles are not missense (\*2 synonymous, \*3 stop-gained, \*4 start-lost, \*17 upstream), so a missense-only tool is blind to exactly the variants that matter.

### 2.3 The calibration is borrowed, and it demonstrably fails

CYP2C19 has no pathogenic variants, so a positive likelihood ratio cannot be computed within the gene. I pooled it with DPYD, giving 135 labelled variants (39 pathogenic, 96 benign) across five CADD bands.

| CADD band | n | P/LP | B/LB | LR+ | Evidence reading |
|:--|:--:|:--:|:--:|:--:|:--|
| 0 to 10 | 65 | 0 | 65 | 0.00 | strong, benign |
| 10 to 20 | 24 | 0 | 24 | 0.00 | strong, benign |
| 20 to 25 | 11 | 5 | 6 | 2.05 | no evidence either way |
| 25 to 30 | 4 | 3 | 1 | 7.38 | n < 10, not used |
| 30 to 50 | 31 | 31 | 0 | infinite | strong, pathogenic |

**130 of those 135 variants come from DPYD**; CYP2C19 contributes five benign variants. This is a DPYD calibration applied to CYP2C19, and the data contain direct proof that the transfer is unsafe: **CYP2C19\*4 (M1V, start-lost) scores CADD 15.9 and therefore lands in the "strong, benign" band**, despite being an established no-function allele. The DPYD pathogenic set is largely truncating (mean CADD 35.5) while the CYP2C19 VUS are largely missense, so the two score distributions are not comparable.

Applied to the 33 VUS (Figure 2), **16 lean benign and 17 remain unchanged; none move toward pathogenic.** Given the \*4 failure I describe the 16 benign-leaning calls as weak provisional calls, not classifications.

![**Figure 2.** The 33 CYP2C19 VUS distributed across the transferred CADD bands. Nothing moves toward pathogenic, and the six variants with no CADD score are large structural deletions that no predictor in this course can evaluate. A blank score means the tool could not see the variant, not that the variant is harmless. Re-plotted from `module5/pgx_classified_CYP2C19.csv`.](figures/fig2_vus_bands.png)

### 2.4 A discordant variant read mechanistically: CYP2C19\*2

The single most instructive result in Lab 5 is a flat contradiction. CYP2C19\*2 (rs4244285) is:

- annotated **synonymous** (P227P) by snpEff,
- scored **CADD 5.7**, which falls in my "strong, benign" band,
- unscored by REVEL and AlphaMissense,
- and simultaneously **the defining no-function allele in CPIC's clopidogrel guideline**.

Both outputs are correct about different things. For a synonymous annotation and a CADD of 5.7 to coexist with complete loss of function, the variant must act **through RNA rather than protein**, and the gene must be **under weak selection**. Both hold. The c.681G>A change creates an aberrant splice site in exon 5; the transcript is mis-spliced, the reading frame shifts, and no functional enzyme is produced (de Morais et al., 1994). The amino acid sequence is untouched, which is exactly what the annotation reports and what a protein-level predictor cannot detect.

I prefer this to the alternative that rs4244285 merely tags a nearby causal variant, because the splice defect has been demonstrated directly, whereas linkage disequilibrium offers no mechanism.

**The actionable consequence for a T2D patient** is therefore: genotype CYP2C19 before starting clopidogrel, and act on the diplotype rather than on any single variant. What I could not determine from variant-level data is **zygosity and phase**: a \*2/\*2 patient is a poor metaboliser while \*2/\*17 is intermediate, and my analysis scores variants one at a time. Nor can a variant table capture **phenoconversion**, where a co-prescribed CYP2C19 inhibitor such as omeprazole makes a genotypic normal metaboliser behave like a \*2 carrier. In a diabetes clinic, where PPI co-prescription is routine, that is a serious practical limit.

---

## 3. Polygenic risk modelling

### 3.1 What the graded model shows, and what it cannot

I applied PGS000334 (Lambert et al., 2021), a 22-SNP Alzheimer's score built on the oligogenic architecture reported for late-onset disease (Zhang et al., 2020) and including APOE epsilon-4 (rs429358) and epsilon-2 (rs7412), to 2,336 individuals from 1000 Genomes phase 3 (EUR 500, AFR 500, EAS 500, SAS 489, AMR 347; The 1000 Genomes Project Consortium, 2015). All 22 variants were present; 13 required allele flipping and none were dropped for mismatch.

**The phenotype is simulated.** 1000 Genomes carries no disease labels, so case status was generated from hidden causal variants in linkage disequilibrium with the score SNPs. The genotypes and the resulting score distributions are real; every AUC and calibration figure below is illustrative of the *method*, not a measurement of anyone's risk. I report it in that spirit, and Section 3.3 shows why the conclusion survives anyway.

| | EUR | AFR | EAS | SAS | AMR | Overall |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| Mean PRS (SD) | +0.04 (0.94) | **+0.74** (1.11) | -0.47 (0.81) | -0.25 (0.78) | -0.10 (0.82) | |
| % above single top-10% cutoff | 7.8 | **28.0** | **2.2** | 4.3 | 6.3 | 10 |
| **AUC** | **0.672** | **0.493** | 0.575 | 0.611 | 0.555 | 0.599 |
| Predicted / observed risk | 0.31 / 0.31 | 0.42 / 0.37 | 0.25 / 0.28 | 0.27 / 0.26 | 0.29 / 0.27 | |
| Observed risk above cutoff | 0.67 | **0.39** | 0.55 | 0.57 | 0.50 | |

**Integration with clinical covariates.** A score is only worth ordering if it adds to what the clinic already knows for free. Added to a model containing age, the PRS lifts cross-validated discrimination by **0.056 AUC (0.562 for age alone, 0.618 combined)**. A gain of that size positions it as a cohort-enrichment tool, not an individual diagnostic.

### 3.2 Discrimination and calibration fail separately, and both must be reported

**Discrimination** asks whether the score ranks cases above controls within a group. It is moderate in EUR (0.672) and **at chance in AFR (0.493)**.

**Calibration** asks whether the predicted risks are the right size. The model sits on the diagonal in EUR (0.31 predicted, 0.31 observed), **over-predicts in AFR** (0.42 against 0.37), and under-predicts in EAS (0.25 against 0.28).

These two failures do not coincide, which is the reason the rubric asks for both. SAS has acceptable discrimination (0.611) and good average calibration. EAS ranks moderately (0.575) but its entire distribution is displaced downward (-0.47 SD) despite a base rate close to EUR's: an AUC alone would hide that shift, because AUC ignores the absolute level. Conversely, a calibration-in-the-large check alone would make AFR look only mildly miscalibrated, because it ignores the collapse in ranking.

The mechanism is linkage disequilibrium decay. The score's SNPs are mostly tags selected in European GWAS, and a tag carries risk information only while it stays correlated with the causal variant: mean tag-to-causal r-squared here is **0.82 in EUR but 0.44 in AFR**, matching the AUC ordering. Allele frequencies compound it, with APOE epsilon-4 at 0.26 in AFR against 0.15 in EUR, so European weights multiply frequencies at which they were never estimated (Martin et al., 2019).

![**Figure 3.** Left: PRS distribution by ancestry under PGS000334. The AFR distribution is displaced upward by 0.74 SD, which is a property of allele frequencies rather than of risk. Right: calibration by ancestry. The AFR curve is nearly flat, reaching only 39% observed risk where 67% is predicted, so the ranking has collapsed even where the average looks tolerable.](figures/fig3_prs_ancestry_calibration.png)

A single cross-ancestry threshold is therefore indefensible. It flags **28.0% of African-ancestry individuals and 2.2% of East Asian-ancestry individuals** against a nominal 10%, and the people it flags do not carry equal risk: observed risk above the cutoff is 0.67 in EUR but 0.39 in AFR, barely above the AFR base rate of 0.37. Deployed unchanged, the cutoff would act on ancestry rather than on risk, and it would do so most aggressively for the group least represented in the source GWAS.

### 3.3 The same pattern on real T2D scores

Because the graded cohort is simulated, I repeated the distribution analysis on **real published T2D scores** outside the graded run, extracting the 287 genotypes they require directly from the 1000 Genomes VCFs for the same 2,336 people. I scored two T2D scores trained in different populations: **PGS003733 (European-derived)** and **PGS003735 (African-derived)**. Harmonisation was clean, with zero allele mismatches and zero missing genotypes. There is no T2D phenotype for this cohort, so this reports distribution only and no AUC should be inferred.

| | EUR | AFR | EAS | SAS | AMR |
|:--|:--:|:--:|:--:|:--:|:--:|
| Mean PRS, EUR-derived score | -0.39 | -0.46 | **+0.59** | +0.46 | -0.27 |
| % above whole-cohort top 10% | 3.0 | 2.8 | **21.0** | 18.0 | 3.5 |
| Mean PRS, AFR-derived score | -0.37 | -0.02 | +0.31 | +0.19 | -0.16 |

Two findings make the Section 3.2 conclusion robust rather than an artefact of simulation. First, **a real T2D score shifts just as sharply**, flagging 21% of East Asian-ancestry and 18% of South Asian-ancestry individuals against 3% of European-ancestry individuals. The affected groups are *different* from the Alzheimer's case, which confirms that the direction of the shift depends on each score's specific allele frequencies and not on any general ordering of populations. Second, the two T2D scores **disagree about the same people**: across all 2,336 individuals they correlate at only **r = 0.65**, despite targeting the same disease at the same 287 positions. The ancestry of the source GWAS, not the biology of diabetes, accounts for roughly a third of the variance in an individual's result.

![**Figure 4.** Real T2D polygenic scores across the same 2,336 individuals. The ancestry displacement reproduces on real scores, and the affected groups differ from the Alzheimer's case, which is the evidence that the phenomenon is a property of score construction rather than of the simulated cohort.](figures/fig4_t2d_prs_by_ancestry.png)

---

## 4. Knowledge graph integration

I built a T2D knowledge graph in Neo4j AuraDB containing **77 nodes and 84 relationships**: 25 genes, 19 GWAS lead variants, 17 Reactome pathways (Milacic et al., 2024), 9 approved drugs (Zdrazil et al., 2024), 5 Human Phenotype Ontology terms (Gargano et al., 2024) and 2 diseases. The 25 genes comprise my 15 Part I targets, 7 genes targeted by approved T2D drugs (so that the graph contains real pharmacology rather than only my candidates), and 3 genes from a deliberately included second disease, non-alcoholic fatty liver disease. Every edge carries a `source` and a `confidence` property, because the final question in this section is precisely how much weight each edge can bear.

![**Figure 5.** The T2D knowledge graph, coloured by node kind with degree-scaled nodes. Genes form the hub because every other layer attaches to a gene; pathway nodes are what turn a gene list into a graph, connecting genes such as LPL and PPARG that share no direct edge.](figures/fig5_knowledge_graph.png)

### 4.1 Centrality recovers known pharmacology, which is the wrong kind of success

PageRank ranks GLP1R highest (0.0217), then PPARG, INSR and SLC5A2 (0.0165), then KCNJ11 (0.0139). **Every one is an existing approved drug target.**

This is circularity, and it is the most important thing to say about the result. PageRank rewards highly connected nodes, and the highly connected nodes are the ones I attached drugs and pathways to *because* they are already drugged. A graph built partly from pharmacology will rank pharmacology highly. Read honestly, PageRank here is a wiring check, not a discovery.

To extract something other than a popularity contest, I fused PageRank with the Part I evidence already carried on the gene nodes (GWAS significance, druggability prior, absolute plasma protein effect), each as a percentile rank. Figure 1A shows the result, with approved targets shaded grey. **The informative hit is TNFRSF11A at rank 4**, the only gene in the top five that no approved T2D drug touches. It earns the position on evidence rather than connectivity: a strong credible set (p = 2 x 10^-21) and the **largest plasma protein effect of any gene in my Part I list** (absolute log2 hazard ratio 0.91). TNFRSF11A encodes RANK, the RANKL receptor in the non-canonical NF-kappa-B pathway, an inflammatory axis consistent with the CSF1 and IFNLR1 immune signal I flagged in Part I.

### 4.2 Community structure recovers a known mechanism unprompted

Modularity clustering split the graph into 17 communities that are mechanistically coherent. **Community 6 is the clearest evidence that the graph is behaving**: it contains `{KCNJ11, ABCC8, Glimepiride, Regulation of insulin secretion, rs5215}`. That is a complete self-contained mechanism, namely the two subunits of the K-ATP channel, the pathway they act in, the sulfonylurea that blocks them, and the lead variant shared by both genes. This is the sulfonylurea mechanism recovered from the data without being told, and it identifies the T2D subtype where the defect is insulin secretion rather than insulin resistance. Community 0 is its modern counterpart (GLP1R with semaglutide and liraglutide), and community 2 is a lipid and adipocyte module containing **LPL, PPARG and pioglitazone**, which becomes important in Section 6.

### 4.3 The graph question: can an approved drug for another disease reach T2D?

Querying across shared phenotypes (5 hops) returned one candidate:

```
(Resmetirom) -[TARGETS]-> (THRB)                    ChEMBL mechanism | high
(THRB) -[ASSOCIATED_WITH]-> (Fatty liver disease)   clinical precedence | LOW
(Fatty liver disease) -[HAS_PHENOTYPE]-> (Insulin resistance)   curated | medium
(Insulin resistance) <-[HAS_PHENOTYPE]- (Type 2 diabetes)       HPO | high
```

**Hypothesis: resmetirom, a thyroid hormone receptor-beta agonist approved for MASH, may improve insulin resistance in T2D.** The mechanism is plausible, since THR-beta activation reduces hepatic fat and hepatic fat is tightly coupled to hepatic insulin resistance, which is the shared phenotype the path turned on (Harrison et al., 2024).

**An edge is evidence, not truth, and this path proves it.** Open Targets (Buniello et al., 2025) lists only two HPO terms for T2D, so I added hyperglycemia, obesity and hypertriglyceridemia from the clinical definition at medium confidence. **All three shared phenotypes in this path are edges I curated myself.** The single most interesting result in Lab 7 rests on its least defensible evidence. The `THRB to fatty liver disease` edge is explicitly labelled low confidence because THRB is in the graph only because a drug targets it, not because genetics implicates it: that edge means no more than "someone made a drug for this". With ontology-derived phenotype annotations in place of my judgement, the path would either become trustworthy or vanish, and either outcome would be informative.

---

## 5. Single-cell and spatial context

**This section is a method transfer, and I state its limit before its result.** The Lab 8 analysis was performed on breast cancer, because the CELLxGENE Census query and the Visium section were fixed by the assignment. It therefore provides **no cell-type evidence about any T2D target**. What transfers is the question and the demonstrated capacity to answer it; what does not transfer is the tissue.

The question is the one a target list cannot answer on its own: *does the target actually sit in the cells we assume, and can a drug physically reach it?* On 18,388 cells from 10 donors, processed in Scanpy (Wolf et al., 2018) and integrated with Harmony (Korsunsky et al., 2019; donor silhouette moved from +0.095 to -0.076 while cell-type structure was retained at 0.043), the answer for ERBB2 was clear: expression is **10.7 times higher in malignant cells than in genuine T cells** (0.64 against 0.06), and effectively absent from myeloid and endothelial compartments. Paired Visium data, analysed with Squidpy (Palla et al., 2022), placed that signal in the epithelial-rich domains and showed the tissue to be immune-excluded, with 35 of 36 off-diagonal neighbourhood pairs negative. The target is abundant and tumour-intrinsic, but concentrated exactly where immune effector cells are least present: favourable for an antibody-drug conjugate, unfavourable for a T-cell engager.

![**Figure 6.** Cell-type resolution of a target. Left: integrated embedding coloured by corrected cell-type label. Right: target expression, concentrated in the malignant compartment. The analytical capability demonstrated here is what a T2D target set requires; the tissue is not T2D.](figures/fig6_singlecell_target.png)

### 5.1 The discordant finding: a label error that inverted the biology

Using the provided annotation helper, ERBB2 appeared **highest in T cells** (mean 0.52 across 10,880 cells) and near the bottom in epithelial cells. That is not a believable result, since ERBB2 is an epithelial receptor tyrosine kinase.

The cause is a substring match. The helper tests whether `"t cell"` appears in the cell-type string, and the Census label `"malignant cell"` ends in those characters. **Of the 10,880 cells labelled T, 8,526 (78%) are actually malignant cells.** The same trap catches `"mast cell"` and `"myofibroblast cell"`. I left the provided function unmodified and re-derived the labels in a separate marker-validated cell.

This is the most transferable finding in the entire report, and it is methodological rather than biological. The error was detectable **only because ERBB2's correct answer was known in advance.** For a novel target such as TNFRSF11A or CSF1 there would have been no alarm, and the analysis would have confidently reported the wrong cell type. Any cell-type claim about a T2D target must therefore be validated against canonical markers before it is believed, not after.

### 5.2 What a T2D version would require

For the Part I targets the relevant compartments are known in advance, which makes the design specifiable. LPL is an adipose, macrophage and capillary endothelial gene; CSF1 acts on the monocyte and macrophage lineage; the secretion module recovered in Community 6 (KCNJ11, ABCC8) is pancreatic beta cell. The analysis would require an islet and adipose atlas rather than a tumour atlas, and the question would be whether the inflammatory targets sit in resident macrophages (druggable from circulation) or in a stromal niche. Notably, the Lab 8 cohort did contain 2,473 hepatocytes, and liver is central to T2D, but ERBB2 expression in hepatocytes speaks to antibody-drug conjugate toxicity rather than to diabetes.

---

## 6. Integrated precision-medicine synthesis

### 6.1 What I would tell a precision-medicine team

**Who to test, and how to report it.** Do not deploy a polygenic score with a single cross-ancestry threshold. The evidence against it is quantitative and reproduces on real T2D scores: one cutoff flags 2.2% to 28.0% of a group depending on ancestry, and two scores for the same disease at the same 287 positions agree on individuals at only r = 0.65. The defensible minimum is **ancestry-specific percentiles with ancestry-specific calibration**, and a result report that states measured performance in the patient's own ancestry group. Where performance is at chance, as in AFR here, the honest option is to return no score rather than a number carrying unearned authority. Positioned correctly, a T2D PRS enriches prevention cohorts (it adds 0.056 AUC over age) and does not diagnose individuals.

**Which drug.** Genotype **CYP2C19 before initiating clopidogrel**, which in a T2D population is a frequent prescription. Act on the diplotype, not on single variants, and record concurrent PPI use because phenoconversion can override the genotype. This recommendation rests on CPIC Level A evidence rather than on my calibration, which Section 2.3 shows cannot be trusted at variant level.

**Which target to pursue.** Two candidates survive the integration, and they differ in character:

- **LPL** is the Part I rank-1 target, sits at fused rank 6, and clusters in the graph's lipid and adipocyte module alongside **PPARG and pioglitazone**. That adjacency is the useful finding: approved pharmacology already acts on this module, so the mechanism is drug-reachable. **The critical caveat is one the graph structurally cannot express.** My Part I data show LPL transcript up (log fold change +0.33) but plasma LPL protein *inversely* associated with disease (log hazard ratio -0.84, p = 8 x 10^-168), and the integration flagged the two layers as discordant. Higher circulating LPL is protective. The therapeutic direction is therefore **activation, not inhibition**, and since `ASSOCIATED_WITH` records no effect direction, a graph-driven pipeline would cheerfully propose an inhibitor. This single omission is the most consequential limitation in the whole analysis.
- **TNFRSF11A** is the non-obvious candidate: fused rank 4, the largest plasma protein effect in the Part I set (0.91), a strong credible set (p = 2 x 10^-21), and no approved T2D drug. It is not an obvious diabetes gene, which is precisely the argument for examining it. RANK/RANKL inhibition is already a clinically validated axis in bone disease, so tractability is plausible. It also carries the most risk, since nothing independent corroborates it.

**Which cells to look at.** For either candidate the answer is the adipose, macrophage and islet compartments, and Section 5 establishes that the required analysis must be marker-validated rather than trusted from automated labels. The convergence here is honest but incomplete: variant, risk and network evidence can be placed on the same targets, while the cell-type layer supplies method rather than measurement.

### 6.2 Where the four lines agree, and where they do not

The four labs do not uniformly reinforce one another, and claiming otherwise would misrepresent the evidence. Two convergences are genuine: an **inflammatory axis** appearing independently in Part I proteomics (CSF1, IFNLR1), in the graph's non-obvious hit (TNFRSF11A, RANK/NF-kappa-B) and in the fused ranking; and a **lipid and adipocyte axis** appearing as Part I rank 1 (LPL) and recovered unprompted as graph community 2 with its pharmacology attached. Two divergences are equally real: the pharmacogenomic layer addresses comorbidity rather than diabetes biology, and the cell-type layer is on the wrong tissue. The target nomination is reasonably supported; the cell-type claim is not yet made.

---

## 7. Limitations and validation

**Cross-cohort integration is not real integration.** No patient appears in more than one of these analyses. The targets come from islet and plasma cohorts, the PRS from 1000 Genomes, the graph from aggregated databases, and the single-cell data from unrelated donors. The layers are joined by gene identity, not by patient identity, so nothing here is a multi-modal measurement on one person. That is the single largest structural limitation, and it is the gap a capstone-scale project would close.

**The PRS results rest on a simulated phenotype.** Every AUC and calibration number in Section 3.1 is illustrative, and the Section 3.3 T2D analysis has no phenotype at all, so it reports distribution only. The ancestry conclusion is robust because it reproduces across two independent score families; the performance numbers are not.

**The variant calibration is borrowed and verifiably fails**, so the 16 benign-leaning VUS calls are hypotheses rather than classifications. The AlphaGenome Atlas, which could have scored the 25 non-missense variants every other tool missed, was not run because I had no API key, so the regulatory and splicing layer of this gene is unassessed.

**Graph edges are evidence of varying quality, and the weakest carry the most weight.** The three shared phenotypes driving the resmetirom path are my own curation; the THRB association is pharmacological precedent with no genetic support; `OCCURS_IN` assumes the GWAS lead variant is causal, which it usually is not; and `ASSOCIATED_WITH` encodes no effect direction, which as Section 6.1 shows is the difference between activating and inhibiting LPL.

**The cell-type layer is absent for T2D.** Section 5 supplies method, not measurement, and the labelling bug shows that even a correct pipeline silently produces inverted biology when an annotation helper is trusted.

### 7.1 How I would validate the top finding

Taking **TNFRSF11A** as the top novel finding, in ascending order of cost:

1. **Check the literature and trial registries first.** The graph knows only what I loaded. Query ChEMBL and ClinicalTrials.gov for RANK/RANKL agents with glycaemic endpoints, and specifically for denosumab cohorts, where glucose outcomes may already have been recorded in osteoporosis trials. This is free and could settle the question immediately.
2. **Mendelian randomisation as a free experiment.** Use TNFRSF11A variants as instruments to ask whether genetically altered RANK signalling changes T2D risk, and critically in which direction. Given the LPL lesson, direction is the first thing to establish, not the last. If the genetics point the wrong way the hypothesis dies before any wet work.
3. **Resolve the cell type properly.** Measure TNFRSF11A in an islet and adipose single-cell atlas with marker-based validation of every label, asking whether expression localises to resident macrophages or to stroma. This is the Section 5 method applied to the right tissue.
4. **Replicate the proteomic association in a multi-ancestry cohort.** The 0.91 effect comes from one plasma panel limited to roughly 1,460 proteins. Independent replication in All of Us or the Million Veteran Program, with ancestry-stratified estimates, would test both the effect and whether it is as European-specific as the PRS results suggest.

**What would change my mind.** If Mendelian randomisation shows RANK signalling raising T2D risk rather than lowering it, or shows no effect, then TNFRSF11A is a biomarker of inflammation accompanying diabetes rather than a target, and the inflammatory axis should be pursued through CSF1 instead.

---

## References

Buniello, A., Suveges, D., Cruz-Castillo, C., Llinares, M. B., Cornu, H., Lopez, I., ... Dunham, I. (2025). Open Targets Platform: Facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Research, 53*(D1), D1467-D1475. https://doi.org/10.1093/nar/gkae1128

de Morais, S. M., Wilkinson, G. R., Blaisdell, J., Nakamura, K., Meyer, U. A., & Goldstein, J. A. (1994). The major genetic defect responsible for the polymorphism of S-mephenytoin metabolism in humans. *Journal of Biological Chemistry, 269*(22), 15419-15422.

Gargano, M. A., Matentzoglu, N., Coleman, B., Addo-Lartey, E. B., Anagnostopoulos, A. V., Anderton, J., ... Robinson, P. N. (2024). The Human Phenotype Ontology in 2024. *Nucleic Acids Research, 52*(D1), D1333-D1346. https://doi.org/10.1093/nar/gkad1005

Harrison, S. A., Bedossa, P., Guy, C. D., Schattenberg, J. M., Loomba, R., Taub, R., ... Ratziu, V. (2024). A phase 3, randomized, controlled trial of resmetirom in NASH with liver fibrosis. *New England Journal of Medicine, 390*(6), 497-509. https://doi.org/10.1056/NEJMoa2309000

Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., ... Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods, 16*(12), 1289-1296. https://doi.org/10.1038/s41592-019-0619-0

Lambert, S. A., Gil, L., Jupp, S., Ritchie, S. C., Xu, Y., Buniello, A., ... Inouye, M. (2021). The Polygenic Score Catalog as an open database for reproducibility and systematic evaluation. *Nature Genetics, 53*(4), 420-425. https://doi.org/10.1038/s41588-021-00783-5

Lee, C. R., Luzum, J. A., Sangkuhl, K., Gammal, R. S., Sabatine, M. S., Stein, C. M., ... Scott, S. A. (2022). Clinical Pharmacogenetics Implementation Consortium guideline for CYP2C19 genotype and clopidogrel therapy: 2022 update. *Clinical Pharmacology & Therapeutics, 112*(5), 959-967. https://doi.org/10.1002/cpt.2526

Martin, A. R., Kanai, M., Kamatani, Y., Okada, Y., Neale, B. M., & Daly, M. J. (2019). Clinical use of current polygenic risk scores may exacerbate health disparities. *Nature Genetics, 51*(4), 584-591. https://doi.org/10.1038/s41588-019-0379-x

Milacic, M., Beavers, D., Conley, P., Gong, C., Gillespie, M., Griss, J., ... D'Eustachio, P. (2024). The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Research, 52*(D1), D672-D678. https://doi.org/10.1093/nar/gkad1025

Palla, G., Spitzer, H., Klein, M., Fischer, D., Schaar, A. C., Kuemmerle, L. B., ... Theis, F. J. (2022). Squidpy: A scalable framework for spatial omics analysis. *Nature Methods, 19*(2), 171-178. https://doi.org/10.1038/s41592-021-01358-2

Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. *Nucleic Acids Research, 47*(D1), D886-D894. https://doi.org/10.1093/nar/gky1016

The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature, 526*(7571), 68-74. https://doi.org/10.1038/nature15393

Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: Large-scale single-cell gene expression data analysis. *Genome Biology, 19*, 15. https://doi.org/10.1186/s13059-017-1382-0

Zdrazil, B., Felix, E., Hunter, F., Manners, E. J., Blackshaw, J., Corbett, S., ... Leach, A. R. (2024). The ChEMBL Database in 2023. *Nucleic Acids Research, 52*(D1), D1180-D1192. https://doi.org/10.1093/nar/gkad1004

Zhang, Q., Sidorenko, J., Couvy-Duchesne, B., Marioni, R. E., Wright, M. J., Goate, A. M., ... Visscher, P. M. (2020). Risk prediction of late-onset Alzheimer's disease implies an oligogenic architecture. *Nature Communications, 11*, 4799. https://doi.org/10.1038/s41467-020-18534-1

---

*Underlying work: `module5/` (variant interpretation), `module6/` and `module6/t2d_extra/` (polygenic risk), `module7/` (knowledge graph), `module8/` (single-cell and spatial), `module2/` to `module4/` and `portfolio1/` (Part I target nomination). Figures 1 and 2 were re-plotted for this report from committed lab tables; Figures 3 to 6 are the lab figures as generated. No new analysis was performed.*
