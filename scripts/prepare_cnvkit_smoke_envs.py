#!/usr/bin/env python3
"""Install exact public CNVkit code in a new environment for CBS-only smoke tests.

Scientific dependencies are inherited from the explicit base Python. This is
not a reconstruction of an historical container. No patient data is accessed.
"""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-python', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--version', required=True, choices=['0.9.10', '0.9.13'])
    parser.add_argument('--legacy-cbs-dependencies', action='store_true',
                        help='0.9.10 only: install the tested NumPy/pandas/SciPy subset; HMM remains unsupported')
    args = parser.parse_args()
    base = Path(args.base_python).absolute()
    if not base.is_file():
        parser.error('The base Python executable is missing')
    if args.legacy_cbs_dependencies and args.version != '0.9.10':
        parser.error('--legacy-cbs-dependencies applies only to CNVkit0.9.10')
    out = Path(args.out).absolute()
    env = out / ('cnvkit-' + args.version + '-env')
    if env.exists():
        parser.error('Refusing to overwrite an existing environment: ' + str(env))
    if any((out / name).exists() for name in ('install-' + args.version + '.json',
                                            'install-' + args.version + '.log',
                                            'commands-' + args.version + '.json')):
        parser.error('Existing install evidence found; use a new output directory')
    out.mkdir(parents=True, exist_ok=True)
    pip = [str(env / 'bin/python'), '-m', 'pip', 'install', '--no-deps', '--disable-pip-version-check',
           '--cache-dir', str(out / 'pip-cache')]
    commands = [[str(base), '-m', 'venv', '--system-site-packages', str(env)],
                [*pip, '--report', str(out / ('install-' + args.version + '.json')), 'CNVkit==' + args.version]]
    if args.legacy_cbs_dependencies:
        commands.append([*pip, '--report', str(out / 'install-legacy-cbs-dependencies.json'),
                         'numpy==1.26.4', 'pandas==2.2.3', 'scipy==1.14.1', 'pytz==2025.2', 'pomegranate==1.1.2'])
    record = dict(commands=commands, base_python=str(base), patient_data_access=False,
                  scope='CNVkit CBS local smoke only; historical HMM API is not supported by this environment',
                  runtime_variant='isolated_exact_CNVkit_reusing_recorded_scientific_dependencies')
    metadata = out / ('commands-' + args.version + '.json')
    metadata.write_text(json.dumps(record, indent=2) + '\n')
    with (out / ('install-' + args.version + '.log')).open('x') as stream:
        for command in commands:
            print('Running', command, flush=True)
            try:
                subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
            except subprocess.CalledProcessError as exc:
                record.update(status='failed', returncode=exc.returncode)
                metadata.write_text(json.dumps(record, indent=2) + '\n')
                raise
    record['status'] = 'installed_not_yet_tested'
    metadata.write_text(json.dumps(record, indent=2) + '\n')
    print(env / 'bin/python')


if __name__ == '__main__':
    main()
