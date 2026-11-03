# Build or reuse?

Searched on 2026-11-03: the nf-core pipeline registry, GitHub repositories and code (Nextflow files), and the wider web.

**Decision: build a new pipeline, reusing nf-core modules and conventions.**

## What it has to do

1. Take UniProt accessions as input.
2. Resolve the AlphaFold DB model version from the API.
3. Gate models on pLDDT.
4. Search with Foldseek.
5. Run a sequence search over the same targets.

## What exists

| Candidate | URL | What it does | Requirements met |
|---|---|---|---|
| nf-core/proteinfold | https://github.com/nf-core/proteinfold | Predicts structures from FASTA (AF2, AF3, ColabFold, ESMFold, Boltz, …); optional Foldseek easy-search on the top model | 4 only |
| nf-core/proteinannotator | https://github.com/nf-core/proteinannotator | Sequence annotation (HMMER Pfam/FunFam, InterProScan, s4pred) | none |
| nf-core/proteinfamilies | https://github.com/nf-core/proteinfamilies | Builds protein families (MMseqs2 clustering, MSA, HMMs) | none |
| ProtPen | https://github.com/ProtPen/ProtPen | AlphaFold model download → Foldseek → UniProt enrichment, plus EggNOG-mapper. Python/SLURM, not Nextflow | 1, partly 2, 4; partly 5 (EggNOG orthology, not a matched search) |
| UCLOrengoGroup/domain-annotation-pipeline | https://github.com/UCLOrengoGroup/domain-annotation-pipeline | Splits AlphaFold models into domains, then Foldseek against CATH | 1, 4 (CATH only) |
| seanome/2024-kmerseek-analysis | https://github.com/seanome/2024-kmerseek-analysis | One-off benchmark; hard-codes `model_v4` | 4 |
| phagomics, dbcan4-advanced, ovo, UHVDB/toolkit, Arcadia 2023-polyphosphate | various | Project-specific pipelines that call Foldseek | 4 |

Three targeted searches found almost nothing:

| Search | Results |
|---|---|
| GitHub code search for Nextflow files calling `alphafold.ebi.ac.uk/api/prediction` | 0 |
| Repository search for "foldseek alphafold nextflow" | 0 |
| Code search for `"foldseek easy-search" extension:nf` | 10 files, all in the repositories above |

## Why build rather than reuse

- Each existing pipeline covers only part of the job.
  - nf-core/proteinfold has Foldseek, but it predicts its own structures from FASTA.
  - ProtPen and the Orengo pipeline fetch AlphaFold models, but they are not nf-core-style Nextflow, and the Orengo images are tagged `main-latest`.
- None of them runs a matched sequence search over the same targets. That matched comparison is the scientific claim here.
- Adapting any of them would mean replacing its entry point and most of its analysis. A small DSL2 pipeline built from nf-core modules is less work and easier to maintain.

## Reused

- nf-core modules `foldseek/createdb`, `foldseek/easysearch`, `mmseqs/createdb` and `mmseqs/easysearch`, adapted to find the database prefix inside a staged directory;
- the nf-schema plugin for parameter and samplesheet validation;
- nf-core conventions: one process per module, `versions.yml` from each process, and `test` profiles.

## Written new

- `FOLDSEEK_DATABASES`: no nf-core module exists for it;
- `AFDB_FETCH`, `PLDDT_QC`, `FETCH_PDB_TARGETS`, `DB_TO_FASTA`, `ANNOTATE_HITS`, `SUMMARISE`, `CHECK_CONTROLS` and `SOFTWARE_VERSIONS`.
