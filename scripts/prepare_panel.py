#!/usr/bin/env python3
"""Fetch only locked public sources. Installed equivalent: asconcordance prepare-panel."""
import argparse
from pathlib import Path
import sys

# For direct use from a source checkout, import that checkout rather than relying
# on the caller's PYTHONPATH or the executable's scripts/ import path.
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ascompare.versions import prepare_panel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--versions', nargs='+')
    parser.add_argument('--lock', type=Path)
    args = parser.parse_args()
    try:
        print('Saved ' + str(prepare_panel(args.out, args.versions, args.lock)))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
