# nextflow-htsget

Minimal Docker-based Nextflow workflow that fetches GA4GH HTSGET data three ways:

| Method | Tool | What it does |
| --- | --- | --- |
| `cli` | `htsget get` | Streams the BAM demo sample with the official client. |
| `python` | `bin/fetch_htsget.py` | Resolves the discovery document and writes the first data block. |
| `curl` | `curl` | Saves the raw discovery JSON returned by the endpoint. |

## Samplesheet

Inputs come from `samplesheet.csv` (override with `--samplesheet`). Each row
needs:

- `id`: identifier used to derive defaults
- `name`: optional label for logging (defaults to `id`)
- `filetype`: one of `cli`, `python`, or `curl`
- `uri`: the HTSGET discovery endpoint to fetch
- `filename`: file name to emit for that method

The included example points at the GA4GH demo service.

## Run

All processes run inside `docker.io/qclayssen/htsget-demo:latest`. Build and
push that image if you change the code.

Execute locally or on Seqera Platform:

```bash
nextflow run main.nf --outdir results
```

Outputs land in `<outdir>/cli`, `<outdir>/python`, and `<outdir>/curl`.
