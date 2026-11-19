"""Identify the selected sample's original CNVkit version from local declarations."""
import hashlib
import json
from pathlib import Path
import re


def identify_source_cnvkit(sample):
    """Populate the effective standard on a run-local copy before hashing its config."""
    declared = sample.get('source_cnvkit_version')
    evidence = {'declared': declared, 'origin': 'unavailable'}
    sample['source_cnvkit_evidence'] = evidence
    sample.pop('_source_cnvkit_issue', None)
    version = None
    try:
        def checked(value):
            if value is None or value == '':
                return None
            if not isinstance(value, str) or not re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+[a-zA-Z0-9.+_-]*', value):
                raise ValueError('Exact original CNVkit version is missing or invalid')
            return value
        version = checked(declared)
        spec = sample['source_result']
        root = spec.get('result_dir') or spec.get('run_dir')
        files = sorted(Path(root).glob('*_run_metadata.json')) if root else []
        if len(files) > 1:
            raise ValueError('Multiple Suite run metadata files in the supplied source result directory')
        if files:
            path = files[0]
            if not path.is_file() or path.is_symlink():
                raise ValueError('Expected a regular Suite run metadata JSON')
            with path.open('rb') as handle:
                data = handle.read(4 * 1024 * 1024 + 1)
            if len(data) > 4 * 1024 * 1024:
                raise ValueError('Suite run metadata exceeds 4 MiB')
            evidence.update(path=str(path), sha256=hashlib.sha256(data).hexdigest())
            raw = json.loads(data)
            if not isinstance(raw, dict):
                raise ValueError('Suite run metadata must be a JSON object')
            release = raw.get('AmpliconSuite-pipeline_version')
            if release and str(release).removeprefix('v') != sample['source_version'].removeprefix('suite-'):
                raise ValueError('Original metadata Suite version differs from sample source_version')
            build = raw.get('ref_genome')
            if build and build != sample.get('reference_build'):
                raise ValueError('Original metadata reference build differs from sample reference_build')
            command = str(raw.get('cnvkit_cmd') or '').strip().lower()
            called = command not in ('', 'na', 'n/a', 'none', 'not run', 'skipped')
            evidence['cnvkit_command_recorded'] = called
            if called:
                recorded = checked(raw.get('cnvkit_version'))
                evidence['recorded'] = recorded
                if recorded and version and recorded != version:
                    raise ValueError('Original metadata CNVkit version conflicts with sample declaration')
                if recorded:
                    version = recorded
                    evidence['origin'] = 'source_result_run_metadata'
        if not version:
            raise ValueError('Provide the exact source_cnvkit_version or original Suite metadata recording the CNVkit run')
        if evidence['origin'] == 'unavailable':
            evidence['origin'] = 'sample_declaration'
    except (ValueError, OSError, UnicodeError) as exc:
        sample['_source_cnvkit_issue'] = str(exc)
        version = None
    sample['source_cnvkit_version'] = version
