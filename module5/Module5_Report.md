# Variant Interpretation in CYP2C19

*BIOT 6900 · Module 5 Assignment · Rakesh Vayigandla · October 2026*

**Gene.** CYP2C19. I was not assigned a gene on the allocation list, so I chose it myself. CYP2C19 activates clopidogrel and clears proton-pump inhibitors, voriconazole and several antidepressants. **Fetch and quality checks.** The fetch returned 80 ClinVar variants (38 "other", mostly drug response; 33 VUS; 9 B/LB; **0 P/LP**). Both anchors were present (CYP2C19\*2, rs4244285; \*17, rs12248560), and the score columns were populated (CADD 68/80, REVEL 54, AlphaMissense 54, gnomAD 55). **Fallbacks, disclosed.** (1) CYP2C19 has no P/LP variants, so I pooled it with DPYD to calibrate. (2) I had no AlphaGenome key, so the Atlas was skipped.

## 1. Thresholds and calibration

**Threshold: MAX_AF = 0.5.** CYP2C19 is not a disease gene. Its important alleles are common because losing CYP2C19 is harmless until a patient takes a drug that depends on it, so there is no purifying selection to keep them rare. In my data, \*2 has an allele frequency of 16.8% and \*17 (increased function) 20.0%. A credible maximum therefore has to sit above any functional allele, and 0.5 does so. It removed only 2 of 80 variants, rs3758581 (I331V, 95.2%) and rs17885098 (synonymous, 91.9%). Both are ClinVar Benign, and at that frequency the reference base is the minor allele. **No P/LP or drug-response variant was removed.**

**The 0.001 contrast would have discarded 19 variants**, including \*2 (rs4244285), \*17 (rs12248560), \*3 (W212X, rs4986893, 0.25%), \*4 (M1V, rs28399504, 0.24%), \*8 (W120R, rs41291556, 0.16%) and nine more drug-response variants. In other words, a rare-disease filter would have deleted almost every allele that CPIC uses to guide clopidogrel prescribing.

**Predictor: CADD.** It scores every variant class (85% of my variants), whereas REVEL and AlphaMissense score only missense variants (67.5%). Most CYP2C19 function alleles are not missense (\*2 is labelled synonymous, \*3 stop-gained, \*4 start-lost, \*17 upstream), so a missense-only tool would miss the variants that matter.

**Calibration was not possible on CYP2C19 alone.** It has 0 P/LP and 7 B/LB variants after filtering (5 with a CADD score), so LR+ is undefined. I took the documented route and pooled with DPYD, the only one of the five genes with enough P/LP variants.

**Table 1. Pooled CADD calibration** (CYP2C19 + DPYD; bands of CADD phred).

| Band | n | P/LP | B/LB | LR+ | Evidence |
|:--|:--:|:--:|:--:|:--:|:--|
| 0–10 | 65 | 0 | 65 | 0.00 | strong, benign |
| 10–20 | 24 | 0 | 24 | 0.00 | strong, benign |
| 20–25 | 11 | 5 | 6 | 2.05 | no evidence either way |
| 25–30 | **4** | 3 | 1 | 7.38 | (moderate, pathogenic) **n < 10: not used** |
| 30–50 | 31 | 31 | 0 | ∞ | strong, pathogenic |

**What this rests on.** 135 labelled variants (39 P/LP, 96 B/LB), but **130 of them are from DPYD**; CYP2C19 contributes only 5 benign variants. So this is really a DPYD calibration applied to CYP2C19. The DPYD P/LP variants all have CADD ≥ 23.4 (mean 35.5), which suggests they are mostly truncating, whereas CYP2C19's VUS are missense. The LRs of 0 and ∞ are not exact either; they mean the band contained no variants of one label, which, with only 39 P/LP variants, could easily be chance. The 25–30 band rests on four variants, and I treated it as no evidence.

There is direct proof that the transfer is unsafe: **CYP2C19\*4 (M1V, start-lost, CADD 15.9) falls in the "strong, benign" band**, yet it is a well-established no-function allele. My confidence in these bands is therefore low. They can say a variant "looks like DPYD's benign variants", not that it is benign in CYP2C19.

**VUS classification.** Of 33 VUS, **16 lean benign** (CADD ≤ 20) and **17 stay put**: 8 fall in the 20–25 band (LR 2.05, no evidence), 3 in the small 25–30 band (L87Q, G135V, R144S), and 6 are unscored structural deletions. **None move towards pathogenic.** Given the caveats above, I would describe the 16 "leaning benign" calls as weak.

## 2. Coverage gaps and the regulatory genome

**No missense predictor could score:** synonymous (7), upstream (5), intron (1), intergenic (1), splice donor (1), stop-gained (1), stop-lost (1), or the 8 large deletions, which have no consequence assigned at all. Ten variants were scored by **no** tool. A blank score here means the tool *could not see* the variant, not that it found the variant harmless.

**Specific unscored variant: rs72558186** (chr10:g.94781999T>A, a splice-donor change at the intron-5 donor site; CYP2C19\*7). CADD, REVEL and AlphaMissense all returned nothing, and it has no gnomAD frequency, yet ClinVar records a drug-response assertion. To classify it I would need **splicing evidence**: a SpliceAI delta score for donor loss as the in-silico line, then an **RNA assay** (a minigene, or RT-PCR from heterozygous liver or lymphoblastoid cells) to show exon skipping or intron retention and the resulting frameshift. If the assay shows a broken transcript, that is functional evidence (PS3) and would justify a no-function call far more strongly than any predictor.

**The Atlas (skipped, no key).** It could have scored the variants every other tool in the course misses. It would have prioritised the 25 non-missense variants. It could ask whether \*17 (rs12248560, upstream) and rs7902257 (upstream, 2.2%) change **promoter activity in liver** (via DNase, ChIP and RNA-seq tracks), and whether rs4244285 and rs72558186 act through **splicing** (via MERGED_SPLICING). Nothing else in the course can do this: CADD gives non-coding variants a number but no mechanism or tissue, and REVEL and AlphaMissense cannot see them at all. I would also have to read the attribution: if conservation (CACTUS or phastCons) carried the score, the mechanism would still be unresolved.

## 3. A discordant variant, read mechanistically: CYP2C19\*2 (rs4244285)

**Rule tripped:** `flag_common_but_plp`. It is common (gnomAD 16.8%) yet carries a functional ClinVar assertion (drug response).

**What each tool said:**

- **snpEff:** synonymous (P227P).
- **CADD:** 5.7, which falls in my "strong, benign" band.
- **REVEL and AlphaMissense:** no score.
- **ClinVar and CPIC:** the defining no-function allele, and the main reason patients get poor clopidogrel activation.

**My explanation: both outputs are right about different things.** For a synonymous call and a CADD of 5.7 to coexist with loss of function, the variant must act **through RNA, not protein**, and the gene must be **under weak selection**. Both are true. c.681G>A creates an **aberrant splice site in exon 5**, so the mRNA is mis-spliced, the reading frame shifts, and no functional enzyme is made. The DNA change is "silent" at the amino-acid level, which is exactly what snpEff reports and what a missense tool cannot see. CADD's low score largely reflects that the site is a synonymous position of modest conservation. That is a fair summary for a gene where losing function costs nothing without a drug, so evolution has not removed the allele.

I prefer this explanation to the alternative, that rs4244285 is just a marker in linkage disequilibrium with a nearby causal variant. The splice defect has been shown directly, and that alternative does not explain how the allele works.

## 4. What I could not determine about \*2

1. **Zygosity and phase.** My analysis scores single variants, but the phenotype depends on the diplotype. A \*2/\*2 patient is a poor metaboliser, and \*2/\*17 is intermediate. rs12769205 (intronic, 17.3%), which also flagged, sits on the same haplotype as \*2. Without phased genotypes I cannot say which alleles are on the same chromosome.
2. **Modifiers, including phenoconversion.** A patient's real CYP2C19 activity also depends on co-medication. Omeprazole and fluvoxamine inhibit CYP2C19 and can make a genotypic normal metaboliser behave like a \*2 carrier. Nothing in a variant table can capture this.

**The single experiment.** Allele-specific RT-PCR (or long-read RNA sequencing) of CYP2C19 transcripts from **human liver samples heterozygous for rs4244285**. If the A allele produces the aberrant exon-5 transcript and the G allele does not, the splicing mechanism is confirmed. **What would change my mind:** if both alleles spliced normally, the effect must come from a linked variant (rs12769205 is the obvious candidate), and I would re-assign the causal call to it.

## References

1. Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M. CADD: predicting the deleteriousness of variants throughout the human genome. *Nucleic Acids Res*. 2019;47(D1):D886-D894. doi:10.1093/nar/gky1016.
2. Ioannidis NM, Rothstein JH, Pejaver V, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *Am J Hum Genet*. 2016;99(4):877-885. doi:10.1016/j.ajhg.2016.08.016.
3. Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381(6664):eadg7492. doi:10.1126/science.adg7492.
4. Landrum MJ, Lee JM, Benson M, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Res*. 2018;46(D1):D1062-D1067. doi:10.1093/nar/gkx1153.
5. Xin J, Mark A, Afrasiabi C, et al. High-performance web services for querying gene and variant annotation. *Genome Biol*. 2016;17:91. doi:10.1186/s13059-016-0953-9 (MyVariant.info).
6. de Morais SM, Wilkinson GR, Blaisdell J, et al. The major genetic defect responsible for the polymorphism of S-mephenytoin metabolism in humans. *J Biol Chem*. 1994;269(22):15419-15422 (CYP2C19\*2 splice defect).
7. Sim SC, Risinger C, Dahl ML, et al. A common novel CYP2C19 gene variant causes ultrarapid drug metabolism relevant for the drug response to proton pump inhibitors and antidepressants. *Clin Pharmacol Ther*. 2006;79(1):103-113. doi:10.1016/j.clpt.2005.10.002 (CYP2C19\*17).
8. Lee CR, Luzum JA, Sangkuhl K, et al. Clinical Pharmacogenetics Implementation Consortium guideline for CYP2C19 genotype and clopidogrel therapy: 2022 update. *Clin Pharmacol Ther*. 2022;112(5):959-967. doi:10.1002/cpt.2526.
