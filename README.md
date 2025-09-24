# nextflow-htsget

Minimal Nextflow workflow that fetches GA4GH HTSGET reads with the bundled Python client and runs a handful of QC utilities, all inside `docker.io/qclayssen/nextflow-htsget:latest`.

## Quick start

```bash
nextflow run . --outdir results --samplesheet samplesheet.csv
```

Provide either a samplesheet (for example the bundled `samplesheet.csv`) or one or more manual HTSGET URLs—never both at the same time.

## Samplesheet format

Provide a CSV or TSV with a header row. Each entry must include:
- `id` (or `name`): used for file naming and reporting.
- `filetype`: one of `bam`, `fastq`, or `vcf`; determines downstream QC steps.
- `uri`: the HTSGET endpoint to request.

Optional columns:
- `reference_name`, `start`, `end` to request a genomic slice (requires `reference_name` when `start`/`end` are set).

Example (`samplesheet.csv` ships with these rows; TSV uses the same columns with tab separators):

```csv
id,name,filetype,uri,reference_name,start,end
ga4gh_demo,NA12878,bam,https://htsget.ga4gh-demo.org/reads/htsnexus_test_NA12878,,,
ga4gh_demo_small,NA12878,bam,https://htsget.ga4gh-demo.org/reads/htsnexus_test_NA12878,1,100,2000
```

QC outputs are published under `<outdir>/qc` alongside the downloaded data in `<outdir>/downloads`.

## Direct URL input (demo)

For quick demos you can skip preparing a sheet and supply one or more HTSGET discovery URLs directly (clear the samplesheet field in the UI, or pass `--samplesheet ''` on the CLI):

```bash
nextflow run . --htsget_urls https://htsget.ga4gh-demo.org/reads/htsnexus_test_NA12878 \
                 --htsget_urls https://htsget.example.org/reads/sample2
```

Entries are treated as BAM fetches by default. Use `--htsget_filetype fastq` (or `vcf`) if the manual URLs should trigger a different QC path.

## Seqera reports

`tower.yml` exposes the MultiQC summary, FastQC HTML reports, samtools/bcftools stats, and the built-in Nextflow execution reports in Seqera Platform. Launch with `-with-tower` to view them on the run page.

## Container image

Rebuild the image after local changes with `scripts/build_docker.sh` (uses `docker buildx` for `linux/amd64`).
