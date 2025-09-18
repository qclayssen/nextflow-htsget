#!/usr/bin/env python3
"""Simple HTSGET downloader script."""

import sys
import requests

def main():
    if len(sys.argv) != 3:
        print("Usage: fetch_htsget.py <uri> <output_file>")
        sys.exit(1)

    uri = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Downloading from {uri} to {output_file}")

    try:
        # For simple HTTP URIs, just download directly
        if uri.startswith('http'):
            response = requests.get(uri, stream=True)
            response.raise_for_status()

            with open(output_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            print(f"Successfully downloaded to {output_file}")
        else:
            print(f"Unsupported URI scheme: {uri}")
            sys.exit(1)

    except Exception as e:
        print(f"Error downloading: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()