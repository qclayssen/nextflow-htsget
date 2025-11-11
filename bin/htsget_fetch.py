#!/usr/bin/env python3
"""Download records from an HTSGET endpoint using the python htsget client."""

from __future__ import annotations

import argparse
import sys

import htsget


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True, help="HTSGET endpoint URL")
    parser.add_argument("--output", required=True, help="Path to write the fetched data")
    parser.add_argument(
        "--format",
        dest="fmt",
        default=None,
        help="HTSGET format hint (e.g. BAM, VCF, FASTQ)"
    )
    parser.add_argument("--reference-name", dest="reference_name", default=None, help="Reference name to subset")
    parser.add_argument("--start", type=int, default=None, help="Start position (0-based, inclusive)")
    parser.add_argument("--end", type=int, default=None, help="End position (0-based, exclusive)")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    params = {}
    if args.fmt:
        params["data_format"] = args.fmt
    if args.reference_name:
        params["reference_name"] = args.reference_name
    if args.start is not None:
        params["start"] = args.start
    if args.end is not None:
        params["end"] = args.end

    url = args.url

    # support htsget URI format by converting it to URI format expected by downloader
    if url.starts_with("htsget://"):
        url = url.replace("htsget://", "https://", 1)

    with open(args.output, "wb") as handle:
        htsget.get(url, handle, **params)

    return 0


if __name__ == "__main__":
    sys.exit(main())
