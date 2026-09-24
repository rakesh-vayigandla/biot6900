"""Build gene-level type 2 diabetes (T2D) tables for Module 2 Part 3.

Three independent cohorts (nothing sample-matched -> integrate at the gene level):

  RNA      GEO GSE76894 - human pancreatic islets from organ donors, T2D vs non-diabetic
           (Solimena et al., Diabetologia 2018; Affymetrix U133 Plus 2.0, RMA-normalized)
  Protein  UK Biobank Pharma Proteomics Project - plasma Olink proteins vs incident T2D
           (Gadd et al., Nature Aging 2024, Supplementary Table 4; Cox HR per SD, age-adjusted)
  Genetics Open Targets Platform - GWAS credible-set evidence score for T2D (MONDO_0005148)

Writes to data/part3_t2d/:
  t2d_rna.tsv      gene, log2fc, pval    (Welch t-test, T2D vs non-diabetic islets)
  t2d_protein.tsv  gene, log2hr, pval    (log2 hazard ratio per SD of plasma protein)
  t2d_gwas.tsv     gene, gwas_score      (0-1; genes absent here have no GWAS evidence)

Run from the repo root:  python prepare_t2d_data.py
(needs openpyxl for the protein table:  pip install openpyxl)
"""
import gzip
import io
import json
import os
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

RAW_DIR = "data/part3_t2d/raw"
OUT_DIR = "data/part3_t2d"

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76894/matrix/GSE76894_series_matrix.txt.gz"
GPL_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=data"
PROT_URL = ("https://media.springernature.com/original/springer-static/esm/"
            "art%3A10.1038%2Fs43587-024-00655-7/MediaObjects/43587_2024_655_MOESM3_ESM.xlsx")
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
T2D_ID = "MONDO_0005148"


def fetch(url, name):
    """Download once into RAW_DIR, then reuse the local copy."""
    path = os.path.join(RAW_DIR, name)
    if not os.path.exists(path):
        os.makedirs(RAW_DIR, exist_ok=True)
        print(f"downloading {name} ...")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req) as r, open(path, "wb") as f:
            f.write(r.read())
    return path


def build_rna():
    with gzip.open(fetch(GEO_URL, "GSE76894_series_matrix.txt.gz"), "rt") as f:
        lines = f.read().splitlines()
    status = next(l for l in lines if l.startswith("!Sample_characteristics_ch1") and "diabetes status" in l)
    status = [s.strip('"').split(": ")[-1] for s in status.split("\t")[1:]]
    start = lines.index("!series_matrix_table_begin") + 1
    end = lines.index("!series_matrix_table_end")
    expr = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t", index_col=0)
    expr.columns = status  # "T2D" or "ND"

    # probe -> gene symbol; drop probes that map to no gene or to several genes
    gpl_lines = open(fetch(GPL_URL, "GPL570.txt")).read().splitlines()
    header = next(i for i, l in enumerate(gpl_lines) if l.startswith("ID\t"))
    gpl = pd.read_csv(io.StringIO("\n".join(gpl_lines[header:])), sep="\t",
                      usecols=["ID", "Gene Symbol"], dtype=str).dropna()
    gpl = gpl[~gpl["Gene Symbol"].str.contains("///")].set_index("ID")["Gene Symbol"]
    expr = expr.loc[expr.index.intersection(gpl.index)]

    # one probe per gene: keep the most highly expressed probe
    order = expr.mean(axis=1).sort_values(ascending=False).index
    expr = expr.loc[order]
    expr.index = gpl.loc[order].values
    expr = expr[~expr.index.duplicated(keep="first")]

    t2d, nd = expr.loc[:, expr.columns == "T2D"], expr.loc[:, expr.columns == "ND"]
    out = pd.DataFrame({"gene": expr.index,
                        "log2fc": (t2d.mean(axis=1) - nd.mean(axis=1)).values,  # data already log2 (RMA)
                        "pval": stats.ttest_ind(t2d, nd, axis=1, equal_var=False).pvalue})
    return out.sort_values("pval").reset_index(drop=True), t2d.shape[1], nd.shape[1]


def build_protein():
    tab = pd.read_excel(fetch(PROT_URL, "Gadd2024_NatAging_SuppTables.xlsx"),
                        sheet_name="Supplementary Table 4", header=7)
    tab = tab[tab["Outcome"] == "Type 2 diabetes"].copy()
    tab["gene"] = tab["Predictor"].str.split(".").str[0]  # "GDF15.Q99988.OID20251.v1" -> GDF15
    tab["log2hr"] = np.log2(tab["HR"])
    tab = tab.rename(columns={"P.Value": "pval"})
    # CXCL8, IL6 and TNF were measured on 3 panels each: keep the strongest measurement
    tab = tab.sort_values("pval").drop_duplicates("gene")
    return tab[["gene", "log2hr", "pval"]].reset_index(drop=True), int(tab["N Cases"].iloc[0])


def build_gwas():
    query = """query($id:String!, $i:Int!) { disease(efoId:$id) {
      associatedTargets(page:{index:$i, size:3000}) { count
        rows { target { approvedSymbol } datasourceScores { id score } } } } }"""
    path = os.path.join(RAW_DIR, "opentargets_T2D.json")
    if not os.path.exists(path):
        os.makedirs(RAW_DIR, exist_ok=True)
        print("querying Open Targets ...")
        rows, i = [], 0
        while True:
            body = json.dumps({"query": query, "variables": {"id": T2D_ID, "i": i}}).encode()
            req = urllib.request.Request(OT_URL, data=body, headers={"Content-Type": "application/json"})
            page = json.load(urllib.request.urlopen(req))["data"]["disease"]["associatedTargets"]
            rows += page["rows"]
            if not page["rows"] or len(rows) >= page["count"]:
                break
            i += 1
        json.dump(rows, open(path, "w"))
    rows = json.load(open(path))
    recs = [{"gene": r["target"]["approvedSymbol"], "gwas_score": s["score"]}
            for r in rows for s in r["datasourceScores"] if s["id"] == "gwas_credible_sets"]
    return pd.DataFrame(recs).sort_values("gwas_score", ascending=False).reset_index(drop=True), len(rows)


if __name__ == "__main__":
    rna, n_t2d, n_nd = build_rna()
    prot, n_cases = build_protein()
    gwas, n_ot = build_gwas()

    rna.to_csv(f"{OUT_DIR}/t2d_rna.tsv", sep="\t", index=False)
    prot.to_csv(f"{OUT_DIR}/t2d_protein.tsv", sep="\t", index=False)
    gwas.to_csv(f"{OUT_DIR}/t2d_gwas.tsv", sep="\t", index=False)

    print(f"RNA:      {len(rna):>6} genes     (islets: {n_t2d} T2D vs {n_nd} non-diabetic)")
    print(f"Protein:  {len(prot):>6} proteins  (plasma: {n_cases} incident T2D cases, UK Biobank)")
    print(f"Genetics: {len(gwas):>6} genes     with GWAS evidence (of {n_ot} Open Targets T2D targets)")
