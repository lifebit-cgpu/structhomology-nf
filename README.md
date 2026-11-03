# structhomology-nf

**Does a predicted structure find relatives that sequence search misses?**

`structhomology-nf` is a Nextflow DSL2 pipeline that runs two searches for each UniProt accession against **the same set of target proteins**:

- a **structure** search with [Foldseek](https://github.com/steineggerlab/foldseek) using the protein's AlphaFold DB model;
- a **sequence** search with [MMseqs2](https://github.com/soedinglab/MMseqs2) using the same protein's sequence.

Only the search method differs, so for every hit you can see whether structure, sequence, or both found it. You also get the hit's functional annotation and whether that annotation actually tells you anything.

```
samplesheet (UniProt accessions)
  -> AFDB_FETCH            AlphaFold DB record + model; version read from the API (latestVersion)
  -> PLDDT_QC              mean pLDDT, confidence-band fractions, PASS / FAIL gate
  -> target database       foldseek databases <name> | --foldseek_db <dir> | --target_pdb_ids <file>
  -> DB_TO_FASTA           the same targets as amino-acid sequences (+ database version / size / entries)
  -> FOLDSEEK_EASYSEARCH   structure arm          MMSEQS_EASYSEARCH  sequence arm (-s 7.5)
  -> ANNOTATE_HITS         shared / structure-only / sequence-only; RCSB + UniProt annotation
  -> SUMMARISE             per-protein summary + one combined hit table
  -> CHECK_CONTROLS        positive controls must recover known relatives, or the run fails
  -> SOFTWARE_VERSIONS     versions.yml from every process -> software_versions.{yml,tsv}
```

Why build this instead of reusing an existing pipeline? See [docs/build_vs_reuse.md](docs/build_vs_reuse.md). In short:

- nf-core/proteinfold, ProtPen and the Orengo domain pipeline each cover part of this job.
- None of them takes UniProt accessions, resolves the AlphaFold DB model version from the API, gates models on pLDDT, and runs a matched sequence search over the same targets.

---

## Quick start

Requirements:

- Nextflow >= 24.04.2 (tested with 26.10.4);
- Docker, Singularity or Apptainer;
- internet access to alphafold.ebi.ac.uk, files.rcsb.org, data.rcsb.org and rest.uniprot.org.

```bash
# Three positive controls against a 24-entry PDB target set; a few minutes.
nextflow run lifebit-cgpu/structhomology-nf -r main -profile test,docker --outdir results_test
```

The test profile runs three proteins:

| Protein | Fold |
|---|---|
| TPI1 (P60174) | TIM barrel |
| myoglobin (P02144) | globin |
| LDHA (P00338) | Rossmann + LDH/MDH |

They are searched against a small Foldseek database built on the fly from 24 RCSB entries ([assets/test_target_pdb_ids.txt](assets/test_target_pdb_ids.txt)):

- 14 known relatives of the controls;
- 10 unrelated decoys.

The run passes only if `CHECK_CONTROLS` finds that every control recovers a known relative.

## Worked example

**1. The samplesheet.** It has one row per UniProt accession. `group` is optional and must be `query` or `control`.

```csv
accession,group
A6NNL5,query
Q9GZU0,query
P02144,control
```

**2. Run against PDB100.** The run first lists `foldseek databases` and stops if the name you give is not in that list. It then downloads the database once into `--db_cache_dir`:

```bash
nextflow run lifebit-cgpu/structhomology-nf -r main -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --foldseek_db_name PDB \
    --db_cache_dir /data/foldseek_dbs \
    --controls https://raw.githubusercontent.com/lifebit-cgpu/structhomology-nf/main/assets/positive_controls.csv
```

**3. Reuse a database you already have.** Point at the directory that holds it. The pipeline finds the database prefix from its `*_ss.dbtype` file:

```bash
nextflow run lifebit-cgpu/structhomology-nf -r main -profile docker \
    --input samplesheet.csv --outdir results --foldseek_db /data/foldseek_dbs/foldseek/pdb
```

**4. Add proteins and resume.** Finished tasks are taken from the cache, so only the new accessions are fetched and searched:

```bash
nextflow run lifebit-cgpu/structhomology-nf -r main -profile docker --input samplesheet_more.csv --outdir results -resume
```

**5. Read the results:**

- `summary/per_protein_summary.tsv` has one row per protein;
- `summary/combined_hits.tsv` has every top hit with both arms side by side;
- `controls/control_check.tsv` has the control verdicts.

### Samplesheet validation

The samplesheet is validated against [assets/schema_input.json](assets/schema_input.json) with nf-schema before anything runs. A malformed sheet stops the run with one message per bad row, for example:

```
* --input (bad.csv): Validation of file failed:
    -> Entry 4: Detected duplicate entries: [accession:P00338]
    -> Entry 1: Error for field 'accession' (TPIS_HUMAN): "TPIS_HUMAN" does not match regular expression [...] (Column 'accession' must be a canonical UniProtKB accession in upper case, e.g. P60174 or A0A1B0GUU1 (no isoform suffix ...
    -> Entry 2: Error for field 'group' (controls): Expected any of [[query, control]] ...
```

## Method details

### AlphaFold DB model version

`AFDB_FETCH` calls `https://alphafold.ebi.ac.uk/api/prediction/<accession>`.

- The API can return several models, including isoforms such as `Q8IXR9-2`. The pipeline keeps the record whose `uniprotAccession` exactly matches the query, which is the canonical `AF-<acc>-F1` model.
- It downloads the model from that record's `pdbUrl` (or `cifUrl`). The version comes from the record's `latestVersion` and is never hard-coded.
- The pLDDT file is taken from `plddtDocUrl`, or built from the entry ID and `latestVersion` when that field is missing.
- Everything the pipeline used is saved in `afdb/<acc>/<acc>.record.json`: the entry ID, model version, creation date, URLs and fetch time.
- The AlphaFold file server rate-limits requests. Downloads therefore retry with exponential backoff, and at most `--afdb_max_forks` (default 3) run at once.

### pLDDT gate (threshold and rationale)

A model is searched only if **mean pLDDT >= 70** *and* **at least 50% of residues have pLDDT >= 70**.

- **Why 70.** It is AlphaFold's own boundary between "confident" and "low" confidence, above which the backbone is generally reliable. Foldseek's 3Di alphabet encodes backbone geometry, so a model with an unreliable backbone gives an unreliable 3Di string.
- **Why also a fraction.** The fraction rule keeps a well-folded domain with a long disordered tail, which mean pLDDT alone can pull below 70. It still rejects models that are mostly low-confidence.

The bands follow AlphaFold DB: very high > 90, confident 70–90, low 50–70, very low <= 50. `PLDDT_QC` also reports the longest run of residues with pLDDT >= 70.

A low pLDDT does not necessarily mean a bad prediction. Many low-pLDDT regions are genuinely disordered. The gate only answers one question: is there a confident fold to search with?

### Target database and the matched sequence search

By default the pipeline uses the **`PDB`** entry from `foldseek databases` (Foldseek's PDB100). This is what our development run recorded:

| Property | Value |
|---|---|
| Snapshot | PDB 2025-01-01 (`PDB_DATE 250101`) |
| Size on disk | 4.3 GB (4,563,898,557 bytes) |
| Entries | 324,204 chains from 151,880 PDB IDs |
| Download | 2.33 GB archive |

`DB_TO_FASTA` runs `foldseek convert2fasta` on the same database. **MMseqs2 therefore searches exactly the same chains as Foldseek**, under the same IDs. Every run records the database name, version string, size in bytes and entry count in `database/<db>.db_info.json`, so a different snapshot is visible.

Both arms use the same reporting limits (`--search_evalue`, `--max_seqs`). MMseqs2 runs at `-s 7.5`, its most sensitive setting, so sequence search gets its best chance.

### Hit classification

Only the best alignment per target is kept in each arm.

- **Structure-significant:** Foldseek E-value <= `--sig_evalue` (1e-3), *or* max(qTM, tTM) >= `--min_tmscore` (0.5) over at least `--min_aln_len` (50) aligned residues. TM >= 0.5 usually means the same fold.
- **Sequence-significant:** MMseqs2 E-value <= `--sig_evalue`.
- **Category** of each hit: `both`, `structure_only`, `sequence_only`, or `weak` (found but significant in neither arm).
- **Poorly annotated:** the hit is described as uncharacterised, unknown function, hypothetical, a DUF, or a designed protein, or its UniProt annotation score is 1 with no EC number and no function text. Such hits, and hits to the query protein itself, are **not informative**. A match to another protein of unknown function moves nothing forward.
- **Annotation sources:** RCSB PDB GraphQL gives the entry title, entity description, organism, EC, Pfam, InterPro, CATH, SCOP, ECOD, GO and UniProt mapping. UniProt REST gives the protein name, annotation score and function text. Use `--skip_remote_annotation` to rely on database headers only.

For each protein, `structure_added_class` takes one of these values:

| Class | Meaning |
|---|---|
| `structure_only` | Only structure found informative relatives |
| `structure_extends_sequence` | Both found informative relatives; structure found additional PDB entries |
| `structure_agrees_with_sequence` | Structure found nothing that sequence missed |
| `sequence_only` | Only sequence found informative relatives |
| `only_self_or_poorly_annotated_hits` | Significant hits exist, but none is informative |
| `no_significant_hits` | Neither arm found a significant hit |

`structure_gain` is the number of informative PDB entries found by structure only.

### Positive controls

The controls are listed in [assets/positive_controls.csv](assets/positive_controls.csv):

| Control | Expected relatives |
|---|---|
| TPI1 | yeast, chicken and *T. brucei* TIM |
| MB | sperm-whale Mb, human Hb, soybean leghemoglobin, *Vitreoscilla* Hb |
| LDHA | dogfish, *B. stearothermophilus* and *P. falciparum* LDH, *E. coli* MDH |

A control **passes** if at least one of its top 10 Foldseek hits meets both conditions:

1. it is a known relative, meaning its PDB ID is listed for that control or its annotation matches the control's family regex;
2. it is significant: E < 1e-3 or TM >= 0.5.

The family-regex route is needed because PDB100 keeps one representative chain per identical-sequence cluster, so the representative may not be the listed PDB ID. With `--fail_on_control_failure true` (the default), a failing control stops the run. `control_check.tsv` also records which expected relatives were found by structure only.

## Parameters

This table was generated from `nextflow_schema.json` with `python3 bin/schema_to_markdown.py nextflow_schema.json`. `nextflow run lifebit-cgpu/structhomology-nf --help` prints the same information.

**Input/output options**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--input` | string (file-path) | `required` | CSV samplesheet: column 'accession' (canonical UniProtKB accession, required) and optional column 'group' ('query' or 'control'). |
| `--outdir` | string (directory-path) | `required` | Output directory (local path or cloud bucket URI). |
| `--publish_dir_mode` | string | `copy` | How result files are published to the output directory. |

**AlphaFold DB options**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--afdb_api_url` | string | `https://alphafold.ebi.ac.uk/api/prediction` | AlphaFold DB prediction API endpoint. The model version is read from the canonical record's 'latestVersion' field, never hard-coded. |
| `--afdb_max_forks` | integer | `3` | Maximum number of AlphaFold DB downloads in parallel (the file server rate-limits with HTTP 429). |

**pLDDT gate**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--min_mean_plddt` | number | `70` | Minimum mean pLDDT for a model to be searched. |
| `--min_frac_confident` | number | `0.5` | Minimum fraction of residues with pLDDT >= 70 for a model to be searched. |

**Target database**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--foldseek_db_name` | string | `PDB` | Name of a database offered by `foldseek databases` (e.g. PDB, Alphafold/Swiss-Prot, CATH50). The pipeline lists the available databases and stops if the name is not among them. |
| `--foldseek_db` | string (directory-path) | | Directory holding an existing Foldseek database (the prefix is detected from its *_ss.dbtype file). Skips the download. |
| `--target_pdb_ids` | string (file-path) | | Text file of PDB IDs (one per line, '#' comments allowed). Builds a small Foldseek database from RCSB mmCIF files; used by -profile test. |
| `--db_cache_dir` | string (directory-path) | | Persistent directory where downloaded Foldseek databases are stored and reused by later runs. Defaults to <outdir>/databases. |

**Search options**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--search_evalue` | number | `10` | Maximum E-value of hits reported by both Foldseek and MMseqs2. |
| `--max_seqs` | integer | `1000` | Maximum prefilter hits passed to alignment per query, in both arms. |
| `--mmseqs_sensitivity` | number | `7.5` | MMseqs2 -s sensitivity. 7.5 is the most sensitive setting, so sequence search is given its best chance. |

**Hit classification and annotation**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--sig_evalue` | number | `0.001` | E-value at or below which a hit counts as significant, in both arms. |
| `--min_tmscore` | number | `0.5` | A Foldseek hit also counts as significant when max(query-, target-normalised TM-score) reaches this value over at least --min_aln_len residues (TM >= 0.5 indicates the same fold). |
| `--min_aln_len` | integer | `50` | Minimum alignment length for the TM-score significance route. |
| `--top_n` | integer | `10` | Top structural and top sequence hits per protein carried into the combined table. |
| `--skip_remote_annotation` | boolean | `false` | Do not query RCSB PDB or UniProt; annotate hits from database headers only (offline use). |
| `--rcsb_graphql_url` | string | `https://data.rcsb.org/graphql` | RCSB PDB Data API (GraphQL) endpoint for entry titles, entity descriptions, EC numbers, Pfam/InterPro/CATH/SCOP/ECOD/GO and UniProt mappings. |
| `--uniprot_rest_url` | string | `https://rest.uniprot.org` | UniProt REST base URL for protein names, annotation scores and function text of hits. |

**Positive controls**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `--controls` | string (file-path) | | CSV of positive controls with columns accession, gene, fold, expected_pdb_ids (';'-separated), family_regex, rationale. Controls present in the samplesheet are checked. |
| `--fail_on_control_failure` | boolean | `true` | Stop the run with an error when any positive control fails to recover a known relative. |

**Resources.** These are set in [conf/base.config](conf/base.config) through process labels and capped by `process.resourceLimits` (16 CPUs / 64 GB / 48 h by default; 4 CPUs / 6 GB / 1 h in `-profile test`). To change them, pass your own config, for example `-c my.config` containing `process.resourceLimits = [cpus: 8, memory: 32.GB, time: 24.h]`.

## Profiles

| Profile | Use |
|---|---|
| `docker`, `singularity`, `apptainer` | container engine (one is required) |
| `test` | 3 positive controls vs a 24-entry PDB target set; a few minutes |
| `test_full` | the 20 uncharacterised proteins + 3 controls ([assets/samplesheet_full.csv](assets/samplesheet_full.csv)) against PDB100 |
| `cloudos` | Lifebit CloudOS defaults (Docker, retries); the executor is set by the platform |
| `debug` | keeps work directories, prints host names |

## Outputs

```
<outdir>/
├── afdb/<acc>/                 <acc>.model.pdb, .confidence.json, .fasta, .record.json (entry, version, URLs, fetch time)
├── plddt_qc/                   <acc>.plddt_qc.tsv / .json   (mean, median, band fractions, longest confident run, flag)
├── database/                   <db>.db_info.json (name, version, size_bytes, n_entries), foldseek_databases_list.txt,
│                               targets_manifest.tsv (test / --target_pdb_ids only)
├── foldseek/                   <acc>.m8   raw Foldseek hits (BLAST-tab with header)
├── mmseqs/                     <acc>.m8   raw MMseqs2 hits
├── per_protein/                <acc>.all_hits.tsv, <acc>.top_hits.tsv (annotated), <acc>.summary.json
├── summary/
│   ├── combined_hits.tsv       top-N structural + top-N sequence hits per protein, both arms side by side, annotation
│   ├── per_protein_summary.tsv one row per protein: QC, counts, structure_added_class, structure_gain, best hits
│   ├── plddt_qc_summary.tsv    QC for every protein, including those that failed the gate
│   └── run_summary.json        database info and run-level counts
├── controls/control_check.tsv  PASS / FAIL / NOT_RUN per positive control, with evidence
├── pipeline_info/              software_versions.yml / .tsv, execution_report.html, execution_timeline.html, execution_trace.txt
└── databases/foldseek/         downloaded Foldseek database (storeDir; or under --db_cache_dir)
```

Key columns of `combined_hits.tsv`:

| Columns | Content |
|---|---|
| `query`, `query_group`, `query_gene`, `query_mean_plddt` | the query protein |
| `target`, `pdb_id`, `chain` | the hit |
| `category`, `structure_significant`, `sequence_significant` | how the hit was found |
| `fs_rank`, `fs_evalue`, `fs_bits`, `fs_prob`, `fs_fident`, `fs_alnlen`, `fs_alntmscore`, `fs_qtmscore`, `fs_ttmscore`, `fs_lddt`, `fs_qstart`…`fs_tcov` | Foldseek alignment statistics |
| `seq_rank`, `seq_evalue`, `seq_bits`, `seq_fident`, `seq_alnlen`, `seq_qcov`, `seq_tcov` | MMseqs2 alignment statistics |
| `hit_title`, `hit_description`, `hit_organism`, `hit_uniprot`, `hit_uniprot_name`, `hit_annotation_score`, `hit_ec`, `hit_pfam`, `hit_interpro`, `hit_cath`, `hit_scop`, `hit_ecod`, `hit_go`, `hit_function` | functional annotation |
| `same_protein`, `poorly_annotated`, `informative` | annotation flags |

## Software and containers

Every tool runs in a container pinned to an exact tag; none uses `latest`. Image digests are recorded in the analysis-run files.

| Tool | Version | Container | Processes |
|---|---|---|---|
| Foldseek | 10.941cd33 | `quay.io/biocontainers/foldseek:10.941cd33--h5021889_1` | FOLDSEEK_DATABASES, FOLDSEEK_CREATEDB, FOLDSEEK_EASYSEARCH, DB_TO_FASTA |
| MMseqs2 | 18.8cc5c | `quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0` | MMSEQS_CREATEDB, MMSEQS_EASYSEARCH |
| Python (standard library only) | 3.12.12 | `quay.io/biocontainers/python:3.12.12` | AFDB_FETCH, PLDDT_QC, FETCH_PDB_TARGETS, ANNOTATE_HITS, SUMMARISE, CHECK_CONTROLS, SOFTWARE_VERSIONS |
| Nextflow | >= 24.04.2 (tested 26.10.4) | — | — |
| nf-schema plugin | 2.8.0 | — | parameter and samplesheet validation |

The exact versions used in a run are written to `pipeline_info/software_versions.tsv`, collected from each process's `versions.yml`.

**How the code is organised:**

- Each process is in its own file under `modules/local/` or `modules/nf-core/`.
- The four Foldseek and MMseqs2 modules are adapted from [nf-core/modules](https://github.com/nf-core/modules). The adaptation finds the database prefix inside a staged directory, so one module works with downloaded, user-supplied or freshly built databases.
- The Python helpers in `bin/` are passed to each process as staged input files and run with `python3`. This works on every executor without depending on git file permissions.

## Resume and caching

- Task inputs are deterministic, so `-resume` reuses every finished task. A stub run of `-profile test` followed by `-resume` re-submitted 0 of 22 tasks.
- Downloaded Foldseek databases are kept with `storeDir` (in `--db_cache_dir`, default `<outdir>/databases`) and reused even without `-resume`.
- On resume, the AlphaFold DB model version stays as it was on the first run. Start a fresh run, without `-resume`, to pick up a newer model version.

## This analysis run

The protein set was chosen on 2026-11-03.

**UniProt query** (`https://rest.uniprot.org/uniprotkb/search`, TSV, paginated): `reviewed:true AND organism_id:9606 AND protein_name:uncharacterized`

**How the 20 proteins were chosen:**

- The query returned 483 entries; 479 have a name starting with "Uncharacterized".
- All 479 have an AlphaFold DB model, version v6 according to `latestVersion`.
- 54 pass the pLDDT gate.
- 33 remain after also requiring:
  - length 100–1000 aa;
  - a model covering the full sequence;
  - no family-implying words in the name ("domain", "family", "kinase", "homolog", "-like", "containing").
- 20 were picked to span lengths: up to 4 per length band, ordering non-"Putative" entries first, then by lowest annotation score, then by highest pLDDT.

The 20 accessions:

`Q96KX1 Q6UXQ4 Q5VUE5 A6NNL5 A6NGY3 Q96MD7 A6NEN9 Q9GZU0 Q5U649 B9A014 Q6P047 Q5TEZ5 Q96LI9 Q9NUD7 Q96K31 Q6ZUG5 Q9NWW7 Q8IXR9 Q96MH7 A0A1B0GUU1`

Positive controls: `P60174` (TPI1), `P02144` (MB), `P00338` (LDHA). Both lists are in [assets/samplesheet_full.csv](assets/samplesheet_full.csv).

**Commands to rerun, in order:**

```bash
# 0. Get the code (or let `nextflow run lifebit-cgpu/structhomology-nf -r <commit>` pull it)
git clone https://github.com/lifebit-cgpu/structhomology-nf.git && cd structhomology-nf

# 1. Pipeline self-test: controls vs 24 PDB entries. CHECK_CONTROLS must pass.
nextflow run . -profile test,docker --outdir results_test

# 2. Positive controls against PDB100. Must pass before step 3 is interpreted.
nextflow run . -profile docker --input assets/samplesheet_controls.csv --outdir results_controls \
    --foldseek_db_name PDB --db_cache_dir dbcache --controls assets/positive_controls.csv

# 3. Full set: 20 uncharacterised proteins + 3 controls. Reuses the cached PDB100.
nextflow run . -profile test_full,docker --outdir results_full --db_cache_dir dbcache -resume
```

The analysis report, interactive results page and the interpretation of the results are delivered with the analysis run. They are not part of this repository's code.

## Caveats

- **Structural similarity is not function.** A shared fold can come from homology or from convergence. A hit suggests a hypothesis to test experimentally; it does not assign a function.
- **A predicted model is not an experimental structure.** Regions with low pLDDT may be genuinely disordered rather than badly predicted.
- **A hit to a poorly annotated protein moves nothing forward.** Such hits are flagged `poorly_annotated` and do not count toward `structure_gain`.
- **PDB100 covers only proteins with experimental structures.** Families with no solved structure cannot be found. `--foldseek_db_name Alphafold/Swiss-Prot` gives wider coverage, at about 10x the disk space.

## Citation and licence

See [CITATIONS.md](CITATIONS.md). The code is released under the [MIT licence](LICENSE).
