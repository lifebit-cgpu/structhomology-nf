# Analysis run: structure vs sequence search for uncharacterized human proteins

This document covers the **structhomology-nf** pipeline and the analysis run made with it on 2026-11-03. General pipeline usage is in the top-level `README.md`.

* **Repository:** https://github.com/lifebit-cgpu/structhomology-nf
* **Pipeline commit used for every run:** `3470193c595946129c1ad2ac75d97f275db191e2`. The checkout used for the full run is byte-identical to GitHub at that commit (36 of 36 pipeline files). This `analysis/` folder was added afterwards and does not change the pipeline.
* **CloudOS job IDs: none.** CloudOS refused to import the repository ("could not read the GitHub repository ... with the user's GitHub account"). All runs below ran in an interactive Linux session with Apptainer. To run on CloudOS, link a GitHub account that can see `lifebit-cgpu/structhomology-nf`, then launch with Nextflow >= 24.04.2 (not the platform default, 22.10.8) and `-profile docker`.

## 1. Result in one paragraph
All three positive controls (TPI1, MB, LDHA) recovered their known PDB families at rank 1 (TM > 0.98). Myoglobin's distant relatives, soybean leghemoglobin (1BIN) and *Vitreoscilla* haemoglobin (1VHB), were found **only** by structure. For the 20 uncharacterized proteins, MMseqs2 found **no** significant hit (E < 1e-3) in PDB100, while Foldseek found hits for 7. Structure gave credible **new** hypotheses for 2: C21orf140/FAM243A (a Tc1/mariner transposase RNase H-like domain) and C17orf113 (a hAT transposase fold). Both are medium confidence, and both have lost most of their catalytic residues. For 3 proteins (C8orf74, C5orf34, C8orf76) structure confirmed families that Pfam/PANTHER profiles already assign. For 13 it found nothing coherent.

| Rank | Accession | Gene | Structural suggestion | MMseqs2 same family? | Confidence |
|---|---|---|---|---|---|
| 1 | B9A014 | C21orf140 | Tc1/mariner DDE-transposase catalytic (RNase H-like) domain, res 23-249; catalytic triad not conserved (1/3 acidic: Mos1 D284<->D183) -> likely domesticated, non-catalytic | no | medium |
| 2 | A0A1B0GUU1 | C17orf113 | hAT (Hermes-like) transposase-derived fold over res 109-644 plus N-terminal zinc finger (domesticated ZBED-like architecture); catalytic DDE 1/3 conserved (Hermes D180<->D207) | no | medium |
| 3 | Q6ZUG5 | - | Weak local similarity to BLUF photoreceptor / dynein light-chain fragments (res ~190-420) | no | low (likely spurious) |
| 4 | A6NEN9 | CXorf65 | N-terminal ubiquitin-like / PB1-like beta-grasp domain (res 1-82) | no | low |
| 5 | Q6P047 | C8orf74 | Ciliary central-apparatus C1a-32 (CLAMP/CCDC189) family; structure confirms the profile annotation | no | high (fold) / no new function |
| 6 | Q96MH7 | C5orf34 | Polo-box-like domain (res 338-541) resembling the PLK1 polo-box duplicated region | no | medium (fold) / no new function |
| 7 | Q96K31 | C8orf76 | TPR-like helical repeat region | no | low |
| 8-20 | Q9NWW7, Q8IXR9, Q5U649, A6NNL5, Q96LI9, Q9NUD7, Q9GZU0, Q5VUE5, Q5TEZ5, Q6UXQ4, Q96MD7, A6NGY3, Q96KX1 | C2orf42, C12orf56, C12orf60, C15orf61, CXorf58, C20orf96, C6orf62, C1orf53, C6orf163, C2orf66, C9orf85, C5orf52, C4orf36 | No structural neighbour with E<0.05 in a coherent family | - | none |

## 2. Data sources (all accessed 2026-11-03)
| Source | URL | Used for |
|---|---|---|
| UniProtKB REST | https://rest.uniprot.org/uniprotkb/search | protein selection; hit annotation |
| AlphaFold DB API (models v6) | https://alphafold.ebi.ac.uk/api/prediction/{accession} | model version (`latestVersion`), model and confidence files |
| Foldseek PDB100 | `foldseek databases PDB` (archive https://foldseek.steineggerlab.workers.dev/pdb100.tar.gz, 2.33 GB) | target set for both arms: PDB snapshot 2025-01-01 (PDB_DATE 250101); 324,204 chains from 151,880 PDB IDs; 1,533,071,455 bytes (1.43 GiB) on disk |
| RCSB PDB Data API | https://data.rcsb.org/graphql | hit titles, organisms, EC/Pfam/InterPro/CATH/GO, UniProt mapping |
| InterPro / Pfam / PANTHER | https://www.ebi.ac.uk/interpro | sequence-profile context used in interpretation |

**UniProt query:** `reviewed:true AND organism_id:9606 AND protein_name:uncharacterized`. It returned 483 entries; 479 are named "Uncharacterized"; all 479 have an AFDB v6 model; 54 pass the pLDDT gate; 33 are eligible; 20 were selected. Re-checked when this README was written: x-total-results = 483.

**Selected accessions (20):** Q96KX1, Q6UXQ4, Q5VUE5, A6NNL5, A6NGY3, Q96MD7, A6NEN9, Q9GZU0, Q5U649, B9A014, Q6P047, Q5TEZ5, Q96LI9, Q9NUD7, Q96K31, Q6ZUG5, Q9NWW7, Q8IXR9, Q96MH7, A0A1B0GUU1

**Positive controls:** P60174 (TPI1), P02144 (MB), P00338 (LDHA). Expected relatives are in `assets/positive_controls.csv`.

## 3. Software, containers and versions
| Component | Version / container |
|---|---|
| Nextflow | 26.10.4 for the full run; 26.04.6 for the lint, test, resume and controls checks (the pipeline requires >= 24.04.2) |
| nf-schema plugin | 2.8.0 |
| Foldseek | 10.941cd33, `quay.io/biocontainers/foldseek:10.941cd33--h5021889_1` |
| MMseqs2 | 18.8cc5c, `quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0` |
| Python (pipeline scripts, standard library only) | 3.12.12, `quay.io/biocontainers/python:3.12.12` |
| Apptainer (session runtime) | 1.5.4 |
| Report build | Python 3.12, pandas, matplotlib, plotly 5.24.1, reportlab 4.2.5 |

Software versions table emitted by the full run:

| Process | Tool | Version |
|---|---|---|
| AFDB_FETCH | python | 3.12.12 |
| ANNOTATE_HITS | python | 3.12.12 |
| CHECK_CONTROLS | python | 3.12.12 |
| DB_TO_FASTA | foldseek | 10.941cd33 |
| FOLDSEEK_EASYSEARCH | foldseek | 10.941cd33 |
| MMSEQS_CREATEDB | mmseqs | 18.8cc5c |
| MMSEQS_EASYSEARCH | mmseqs | 18.8cc5c |
| PLDDT_QC | python | 3.12.12 |
| SOFTWARE_VERSIONS | python | 3.12.12 |
| SUMMARISE | python | 3.12.12 |
| Workflow | Nextflow | 26.10.4 |
| Workflow | structhomology-nf | 1.0.0 |

## 4. Parameters
The full parameter table is in the top-level `README.md` (generated from `nextflow_schema.json`). **Settings for this run:** all defaults except `--foldseek_db`, which points at the pre-downloaded PDB100. Key defaults are `--min_mean_plddt 70`, `--min_frac_confident 0.5`, `--foldseek_db_name PDB`, `--search_evalue 10`, `--max_seqs 1000`, `--mmseqs_sensitivity 7.5`, `--sig_evalue 0.001`, `--min_tmscore 0.5`, `--min_aln_len 50`, `--top_n 10` and `--fail_on_control_failure true`.

pLDDT gate: mean >= 70 and fraction of residues >= 70 at least 0.5. 70 is AlphaFold's boundary for a reliable backbone, which Foldseek's 3Di alphabet needs; the fraction rule keeps a folded domain that has disordered tails. Significance: E <= 1e-3 in both arms (Foldseek hits can also qualify with TM >= 0.5 over >= 50 residues). The report's hit counts use the E-value criterion only. Control rule: a known relative (by PDB ID or family name) among the top 10 Foldseek hits with E < 1e-3 or TM > 0.5.

## 5. Exact commands, in order
```bash
# 0. Requirements: Java 17+, Nextflow >= 24.04.2, Docker or Apptainer/Singularity, ~8 GB free disk for PDB100
git clone https://github.com/lifebit-cgpu/structhomology-nf.git
cd structhomology-nf && git checkout 3470193c595946129c1ad2ac75d97f275db191e2

# 1. (optional) re-derive the protein selection; the samplesheets it produced are already in assets/
#    (select_proteins.py lives in analysis/ on main, not at 3470193: fetch it from main)
pip install requests pandas
python analysis/select_proteins.py --count-only
python analysis/select_proteins.py --outdir selection

# 2. test profile (3 controls vs a 24-entry PDB subset, a few minutes), then show that -resume caches everything
nextflow run . -profile test,docker --outdir results_test
nextflow run . -profile test,docker --outdir results_test -resume

# 3. positive controls against the full PDB100 (FOLDSEEK_DATABASES lists the databases, then downloads 'PDB')
nextflow run . -profile docker --input assets/samplesheet_controls.csv \
    --controls assets/positive_controls.csv --outdir results_controls

# 4. full set: 20 uncharacterized + 3 controls (same launch dir + -resume reuses the database download)
nextflow run . -profile docker --input assets/samplesheet_full.csv \
    --controls assets/positive_controls.csv --outdir results_full -resume
#    or reuse an existing copy: add  --foldseek_db /path/to/pdb100_dir

# 5. report, interactive page and tables: build_report.py is delivered with the analysis outputs
#    (workspace), not committed here. RUN_DIR must contain results_full/ and deliver/ (published tables
#    plus interpretation.csv and catalytic_residue_check.csv); in_selection.csv, in_control_recovery.csv
#    and in_candidates.csv go in the working directory
RUN_DIR=/path/to/run python build_report.py
```
**Command used for the full run here.** The paths are specific to the session. `session_sandbox.config` only bind-mounts `/proc` and turns off trace/report/timeline for that sandbox; it is not part of the pipeline.
```bash
nextflow run <checkout of 3470193> -profile apptainer -c session_sandbox.config \
  --input assets/samplesheet_full.csv --controls assets/positive_controls.csv \
  --foldseek_db <PDB100 downloaded with 'foldseek databases PDB'> --outdir results_full -resume
```
The first attempt stopped when the session's execution time ran out. The rerun with `-resume` reused 52 completed tasks and ran the remaining 68 (about 9 min in total).

**Interpretation step (not scripted in the repository).** The InterPro look-up, the mapping of catalytic residues through the Foldseek alignments, the structure-gain score and the confidence calls were made in interactive session code. Their results are recorded in `interpretation.csv` and `catalytic_residue_check.csv`. The structure-gain score is strength of the best consistent structural hit x novelty against sequence profiles (1 / 0.5 / 0) x (0 if MMseqs2 finds the same family). The confidence levels are expert judgement, not a statistic.

## 6. Outputs of this analysis (delivered in the workspace)
| File | Content |
|---|---|
| `report.pdf` | written report: introduction, build-vs-reuse, methods, results with figures, hypotheses and experiments, caveats, conclusions |
| `results.html` | self-contained offline page (4.8 MB): combined hit table that you can sort and filter (group, found-by, E-value, TM-score, annotation status, text); per-protein structure-vs-sequence table; Plotly figures |
| `combined_hits.csv` | pipeline combined hit table (204 rows) |
| `per_protein_summary.csv`, `plddt_qc.csv` | pipeline per-protein summary and pLDDT QC |
| `interpretation.csv`, `catalytic_residue_check.csv` | ranking, hypotheses and confidence; catalytic-site mapping |
| `structure_vs_sequence_counts.csv` | significant hits: both / structure-only / sequence-only |
| `control_check.csv`, `control_recovery.csv` | positive-control verdicts and recovery of each expected relative |
| `software_versions.csv`, `pipeline_parameters.csv`, `run_checkout_vs_github.csv` | provenance |
| `uncharacterized_selection.csv` | the 20 selected proteins with pLDDT metrics |
| `build_report.py`, `select_proteins.py`, `figures.zip` | analysis scripts and figures |

## 7. Caveats
* **Structural similarity is not function.** Both leading candidates share a fold with transposases but have lost most of their catalytic residues.
* **A predicted model is not an experimental structure.** Low-pLDDT regions may be genuinely disordered rather than badly predicted. Only 54 of 479 uncharacterized proteins passed the gate.
* **A Foldseek hit to a protein that is itself poorly annotated moves nothing forward.** Such hits are flagged (`poorly_annotated`).
* MMseqs2 here is a single-sequence search. Profile methods (Pfam, PANTHER, HHpred) already explain 3 of the structure-only signals.
* PDB100 holds only experimental structures. AFDB-based Foldseek databases may find more, but mostly against targets with no annotation.
* These are session runs, not CloudOS jobs. The full run's PDB100 version field is blank because the database was supplied as a folder; the version above comes from the database survey.
