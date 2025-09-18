#!/usr/bin/env python3
"""Resolve an HTSGET discovery URL and download the first data block."""

import argparse
from pathlib import Path

import requests


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("uri", help="HTSGET discovery endpoint")
    parser.add_argument("output", help="Destination file")
    args = parser.parse_args()

    discovery = args.uri.replace("htsget://", "https://", 1)
    response = requests.get(discovery, timeout=60)
    response.raise_for_status()
    payload = response.json()

    body = payload.get("htsget", payload)
    urls = body.get("urls") or body.get("accessMethods") or []
    if not urls:
        raise SystemExit("No data URLs available")

    first = urls[0]
    if isinstance(first, dict):
        data_url = first.get("url") or (first.get("accessUrl") or {}).get("url")
    else:
        data_url = first
    if not data_url:
        raise SystemExit("Missing download URL")

    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(data_url, stream=True, timeout=120) as handle:
        handle.raise_for_status()
        with target.open("wb") as out_handle:
            for chunk in handle.iter_content(1024 * 1024):
                if chunk:
                    out_handle.write(chunk)


if __name__ == "__main__":
    main()
