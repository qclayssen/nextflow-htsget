#!/usr/bin/env python3

"""
Fetch data from HTSGET endpoint using the htsget Python library
"""

import sys
import argparse
import requests
import htsget


def main():
    parser = argparse.ArgumentParser(description='Fetch data from HTSGET endpoint')
    parser.add_argument('url', help='HTSGET endpoint URL')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('--reference-name', help='Reference sequence name')
    parser.add_argument('--start', type=int, help='Start position')
    parser.add_argument('--end', type=int, help='End position')
    
    args = parser.parse_args()
    
    try:
        # Use htsget library to fetch the data
        with open(args.output, 'wb') as output_file:
            kwargs = {}
            if args.reference_name:
                kwargs['reference_name'] = args.reference_name
            if args.start:
                kwargs['start'] = args.start
            if args.end:
                kwargs['end'] = args.end
                
            htsget.get(args.url, output_file, **kwargs)
            
    except Exception as e:
        print(f"Error fetching data: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()