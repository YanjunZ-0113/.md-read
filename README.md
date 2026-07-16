# AmpSeq — amplicon indel analysis CLI

Command-line pipeline for paired-end amplicon data against an **expected integrated reference** (not WT editing efficiency).

## Install

```bash
pip install -r ampseq/requirements.txt
```

Needs Python 3.9+.

## Quick start

```bash
# Show examples
python ampseq.py example

# Re-run the AP1 demo sample
python ampseq.py run ^
  -r1 AP1_S28_L001_R1_001.fastq.gz ^
  -r2 AP1_S28_L001_R2_001.fastq.gz ^
  -o ampseq_out ^
  -n AP1 ^
  --use-ap1-defaults
```

On Linux/macOS use `\` instead of `^`. Equivalent:

```bash
python -m ampseq run -r1 R1.fq.gz -r2 R2.fq.gz -o out --use-ap1-defaults
```

## New sample (explicit sequences)

```bash
python ampseq.py run \
  -r1 sample_R1.fastq.gz \
  -r2 sample_R2.fastq.gz \
  -o sample_out \
  --fwd FORWARD_PRIMER \
  --rev REVERSE_PRIMER \
  --insert-ref INSERT_146BP_OR_FASTA \
  --window WINDOW_60BP_OR_FASTA \
  --junction 68
```

Or pass a full amplicon (primers + insert); AmpSeq strips primers to get the insert:

```bash
python ampseq.py run -r1 R1.fq.gz -r2 R2.fq.gz -o out \
  --amplicon full_amplicon.fa \
  --fwd FWD --rev REV \
  --window window.fa --junction 68
```

## Options

| Flag | Meaning |
|---|---|
| `-r1` / `-r2` | Paired FASTQ (`.gz` ok) |
| `-o` | Output directory |
| `--insert-ref` | Primer-trimmed expected insert |
| `--amplicon` | Full amplicon; insert derived via primers |
| `--fwd` / `--rev` | PCR primers |
| `--window` | Quantification window sequence inside insert |
| `--junction` | Dashed line position on insert (0-based) |
| `--no-trim-adapters` | Skip cutadapt |
| `--max-small` | Max bp for small-change allele view (default 2) |
| `--use-ap1-defaults` | Built-in AP1 insert / window / junction=68 |

## Outputs (`-o`)

```
out/
  qc_stats.json
  editing_summary.json
  small_indel_summary.json
  run_config.json
  cutadapt_report.txt          # if trimming on
  trimmed/                     # trimmed FASTQs
  tables/
    unique_inserts.tsv
    top_alleles.tsv
    small_indel_alleles.tsv
  figures/
    fig5_allele_composition.png|pdf
    fig6_indel_spectrum.png|pdf
    fig7_indel_location.png|pdf
    fig8_allele_plot_small_indels.png|pdf
    fig9_allele_plot_all_mapped.png|pdf
```

## Pipeline steps

1. Optional Illumina adapter trim (`cutadapt`)
2. Overlap-merge R1/R2 + PCR primer trim + collapse unique inserts
3. Needleman–Wunsch align (`parasail`) to insert reference; classify indels in window
4. CRISPResso-style allele plots (small changes + all mapped)

## Note on terminology

Reference = **expected integrated sequence**.  
`indel_fraction` = fraction of on-target reads with indels in the quantification window — not Cas9 editing efficiency vs WT.
