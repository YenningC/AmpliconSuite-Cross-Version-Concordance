#!/usr/bin/env python3
"""Run only the synthetic CNVkit integration test and preserve all evidence."""
import argparse
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-python', required=True, help='Python with program test dependencies; also supplies Rscript PATH')
    parser.add_argument('--program', default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--out', required=True, help='Directory containing prepared CNVkit environments and smoke evidence')
    parser.add_argument('--version', required=True, choices=['0.9.10', '0.9.13', '0.9.14'])
    parser.add_argument('--cnvkit-python', help='Optional explicit CNVkit interpreter')
    parser.add_argument('--reference-mode', default='matched_normal', choices=['matched_normal', 'flat'])
    args = parser.parse_args()
    out, base, program = Path(args.out).absolute(), Path(args.base_python).absolute(), Path(args.program).absolute()
    if not base.is_file() or not (program / 'tests/test_concordance_cnv.py').is_file():
        parser.error('Base interpreter or the program synthetic test is missing')
    out.mkdir(parents=True, exist_ok=True)
    python = args.cnvkit_python or (str(base) if args.version == '0.9.14'
                                  else str(out / ('cnvkit-' + args.version + '-env/bin/python')))
    if not Path(python).is_file():
        parser.error('CNVkit interpreter is missing; prepare the environment first')
    name = 'smoke-' + args.version + ('-flat' if args.reference_mode == 'flat' else '') + '-attempt-'
    attempt = 1
    while list(out.glob(name + str(attempt) + '*')):
        attempt += 1
    stem = name + str(attempt)
    extra = dict(CNVKIT_SMOKE_PYTHON=str(Path(python).absolute()), CNVKIT_SMOKE_VERSION=args.version,
                 CNVKIT_SMOKE_REFERENCE_MODE=args.reference_mode,
                 CNVKIT_SMOKE_EXTRA_PATH=str(base.parent),
                 CNVKIT_SMOKE_ARTIFACTS=str(out / (stem + '-artifacts')))
    env = dict(os.environ)
    env.update(extra)
    command = [str(base), '-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_concordance_cnv.py',
               '-k', 'test_actual_cnvkit', '-v']
    record = dict(command=command, environment=extra, label='SYNTHETIC_ONLY', program=str(program),
                  runtime_variant='isolated_exact_CNVkit_reusing_recorded_scientific_dependencies')
    with (out / (stem + '.log')).open('x') as stream:
        result = subprocess.run(command, cwd=program, env=env, stdout=stream, stderr=subprocess.STDOUT)
    record['returncode'] = result.returncode
    (out / (stem + '.json')).write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record))
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
