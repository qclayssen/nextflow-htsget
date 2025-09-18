# qclayssen/nextflow-htsget

**Demonstration pipeline that downloads GA4GH HTSGET data via three methods**

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A524.04.2-brightgreen.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

This pipeline demonstrates three different methods for fetching genomic data using the GA4GH HTSGET protocol:

| Method | Tool | Description |
| --- | --- | --- |
| `cli` | `htsget get` | Uses the official HTSGET command-line client to stream BAM data |
| `python` | `htsget` Python library | Uses the Python HTSGET library to fetch and write data |
| `curl` | `curl` | Retrieves the raw discovery JSON returned by the HTSGET endpoint |

## Usage

> [!NOTE]  
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

The typical command for running the pipeline is as follows:

```bash
nextflow run qclayssen/nextflow-htsget --samplesheet samplesheet.csv --outdir results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>           # Finished results in specified location (defined with --outdir)
.nextflow_log      # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Samplesheet

You will need to create a samplesheet with information about the HTSGET endpoints you want to fetch before running the pipeline. Use the `--samplesheet` parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row:

| Column     | Description                                  |
| --------   | -------------------------------------------- |
| `id`       | Unique sample identifier                     |
| `name`     | Sample name for display (optional)          |
| `filetype` | File type (e.g., `bam`)                     |
| `uri`      | HTSGET endpoint URI to fetch                |

An [example samplesheet](samplesheet.csv) has been provided with the pipeline that points to the GA4GH demo service.

### Running the pipeline

Each process runs inside the `docker.io/qclayssen/nextflow-htsget:latest` container. The pipeline can be run on:

- Your local compute environment
- Seqera Platform (formerly Nextflow Tower)
- Any HPC or cloud compute environment

Execute the pipeline as follows:

```bash
nextflow run qclayssen/nextflow-htsget \
   --samplesheet samplesheet.csv \
   --outdir results \
   -profile docker
```

The pipeline outputs will land in:
- `<outdir>/cli` - BAM files fetched using the HTSGET CLI
- `<outdir>/python` - BAM files fetched using the Python library  
- `<outdir>/curl` - JSON discovery documents retrieved with curl

### Profiles

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at runtime. For more information and to see if your system is available in these configs, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!

### Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

#### `-c`

Specify the path to a specific config file (this is a core Nextflow feature). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Credits

This pipeline was originally written by Quentin Clayssen.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use this pipeline, please cite:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Holger Homann, Johannes Köster, Marvin Langen, Harshil Patel, Pavel Tomecek, Denis Moreno, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
