"""Strict preflight; source alignments, indexes and reference files are read only."""
import csv
import hashlib
import re
from pathlib import Path
from .common import fingerprint

def alignment_fingerprint(path, method):
    if method == "full":
        return fingerprint(path)
    if method != "sampled":
        raise ValueError("alignment_fingerprint must be full or sampled")
    p = Path(path).resolve(strict=True); st = p.stat(); h = hashlib.sha256()
    with p.open("rb") as f:
        for offset in sorted(set([0, max(0,st.st_size//2-524288),max(0,st.st_size-1048576)])):
            f.seek(offset); h.update(str(offset).encode()); h.update(f.read(1048576))
    return {"path":str(p),"bytes":st.st_size,"mtime_ns":st.st_mtime_ns,"sampled_sha256":h.hexdigest(),"method":"sampled SHA256 first/middle/last 1MiB + size + mtime; NOT full-file checksum"}


ID = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]{0,99}$')


def safe_id(value):
    if not ID.fullmatch(value) or value in {'.', '..'}:
        raise ValueError('Unsafe/empty ID: use 1–100 letters, numbers, dot, dash or underscore')
    return value


def read_samples(manifest, full_mode=False):
    manifest = Path(manifest).resolve()
    with manifest.open() as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'sample_id','alignment','index','reference_build'} | ({'normal_alignment','normal_index'} if full_mode else {'seeds_bed'})
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Manifest needs: ' + ', '.join(sorted(required)))
        rows = list(reader)
    if not rows:
        raise ValueError('Empty manifest')
    seen = set()
    for row in rows:
        safe_id(row['sample_id'])
        if row['sample_id'] in seen:
            raise ValueError('Duplicate sample_id')
        seen.add(row['sample_id'])
        for key in (('alignment','index','normal_alignment','normal_index') if full_mode else ('alignment','index','seeds_bed')):
            p = Path(row[key]).expanduser()
            row[key] = str((manifest.parent / p).resolve()) if not p.is_absolute() else str(p.resolve())
    return rows


def read_seeds(path, lengths, max_seeds=20, max_total_bp=10000000, allow_empty=False):
    seeds = []
    seen = set()
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.startswith(('#','track ','browser ')):
                continue
            fields = line.split()
            if len(fields) != 4:
                raise ValueError(f'BED line {n}: exactly chrom start end seed_id required')
            chrom,start,end,name = fields
            start,end = int(start),int(end)
            safe_id(name)
            if name in seen:
                raise ValueError('Duplicate seed_id')
            if chrom not in lengths or not 0 <= start < end <= lengths[chrom]:
                raise ValueError(f'BED line {n}: invalid BED0 coordinate or unknown contig')
            seen.add(name)
            seeds.append({'seed_id':name,'chrom':chrom,'start':start,'end':end})
    if (not seeds and not allow_empty) or len(seeds)>max_seeds or sum(x['end']-x['start'] for x in seeds)>max_total_bp:
        raise ValueError('Seed count/total length outside configured small-region limits')
    return seeds


def reference_files(config):
    ref = config['reference']
    root = Path(ref['aa_data_repo']).resolve()
    build_dir = root / ref['build']
    listing = build_dir / 'file_list.txt'
    mapping = {}
    for line in listing.read_text().splitlines():
        fields = line.split('#')[0].split()
        if len(fields)>=2:
            mapping[fields[0]] = (build_dir / fields[1]).resolve(strict=True)
    if 'fa_file' not in mapping:
        raise ValueError('AA file_list.txt has no fa_file')
    fasta = Path(ref['fasta']).resolve(strict=True)
    if mapping['fa_file'] != fasta and fingerprint(mapping['fa_file'])['sha256'] != fingerprint(fasta)['sha256']:
        raise ValueError('Selected FASTA differs from AA reference FASTA')
    return sorted(set([listing.resolve(), fasta, Path(str(fasta)+'.fai')] + list(mapping.values())))


def validate(config, samples):
    import pysam
    ref = config['reference']
    if ref['build'] not in ('hg19','GRCh37','GRCh38'):
        raise ValueError('Simplified adapter supports explicit hg19, GRCh37 or GRCh38 only')
    files = reference_files(config)
    fasta = pysam.FastaFile(ref['fasta'])
    lengths = dict(zip(fasta.references, fasta.lengths))
    results, sources, md5_cache = [], {}, {}
    for file in files:
        sources[str(file)] = fingerprint(file)
    for row in samples:
        seeds = []
        try:
            if row['reference_build'] != ref['build']:
                raise ValueError('Declared sample reference build differs from configured build')
            for key in (('alignment','index') if config.get('analysis_scope')=='full_suite' else ('alignment','index','seeds_bed')):
                p = Path(row[key])
                if not p.is_file() or (not p.stat().st_size and not (key=='seeds_bed' and config.get('allow_empty_seeds',False))):
                    raise ValueError(f'Missing/empty {key}')
                if str(p) not in sources:
                    sources[str(p)] = alignment_fingerprint(p, config.get('alignment_fingerprint','full')) if key=='alignment' else fingerprint(p)
            if config.get('analysis_scope') != 'full_suite':
                seeds = read_seeds(row['seeds_bed'], lengths, config.get('max_seeds',20), config.get('max_seed_bp',10000000),allow_empty=config.get('allow_empty_seeds',False))
            pysam.quickcheck(row['alignment'])
            with pysam.AlignmentFile(row['alignment'], reference_filename=ref['fasta'], index_filename=row['index']) as bam:
                if bam.header.to_dict().get('HD',{}).get('SO') != 'coordinate':
                    raise ValueError('Alignment is not declared coordinate sorted')
                if not bam.has_index() or not bam.check_index():
                    raise ValueError('Index cannot be loaded')
                sq = bam.header.to_dict().get('SQ',[])
                if not sq:
                    raise ValueError('Alignment has no sequence dictionary')
                missing_md5 = 0
                missing_contigs = []
                allowed_missing = set(ref.get('allowed_missing_reference_contigs',[]))
                for seq in sq:
                    name = seq['SN']
                    if name not in lengths:
                        if name in allowed_missing:
                            missing_contigs.append(name)
                            continue
                        raise ValueError('Alignment contig absent from reference; no explicit limited-scope exception')
                    if lengths[name] != seq['LN']:
                        raise ValueError('Alignment sequence length differs from reference')
                    if seq.get('M5'):
                        if name not in md5_cache:
                            h = hashlib.md5()
                            for pos in range(0, lengths[name], 4*1024*1024):
                                h.update(fasta.fetch(name,pos,min(pos+4*1024*1024,lengths[name])).upper().encode())
                            md5_cache[name] = h.hexdigest()
                        if seq['M5'].lower() != md5_cache[name]:
                            raise ValueError('Alignment M5 mismatch with reference sequence')
                    else:
                        missing_md5 += 1
                for seed in seeds:
                    if seed['chrom'] not in bam.references:
                        raise ValueError('Seed contig absent from alignment')
                    next(bam.fetch(seed['chrom'],seed['start'],seed['end']), None)
                warnings = []
                if missing_contigs:
                    warnings.append('LIMITED REFERENCE: extra alignment contigs missing from FASTA; human-only exploratory run; connected structures entering these contigs are not evaluable')
                if missing_md5:
                    warnings.append('M5 absent on some contigs; build identity relies on declared build, lengths and supplied FASTA provenance')
                results.append({**row,'status':'ok' if seeds or config.get('analysis_scope')=='full_suite' else 'no_seeds','reason':'' if seeds or config.get('analysis_scope')=='full_suite' else 'Upstream seed selection produced no regions; AA/AC not evaluated','seeds':seeds,'warnings':warnings,'reference_lengths':lengths,'missing_reference_contigs':missing_contigs})
        except Exception as exc:
            results.append({**row,'status':'failed','stage':'input_check','reason':str(exc),'seeds':seeds})
    fasta.close()
    return results, list(sources.values())


def validate_full(config, samples):
    expanded=[]
    for n,row in enumerate(samples):
        expanded.append(row)
        expanded.append({**row,'sample_id':f'NORMAL_CHECK_{n}','alignment':row['normal_alignment'],'index':row['normal_index']})
    checked,sources=validate({**config,'analysis_scope':'full_suite'},expanded)
    results=[]
    for n,row in enumerate(samples):
        tumor,normal=checked[2*n:2*n+2]
        tumor['normal_check']=normal
        try:
            if Path(row['alignment']).samefile(row['normal_alignment']):
                normal.update(status='failed',reason='Tumor and normal point to the same alignment file')
        except OSError:
            pass
        if normal['status']!='ok':
            tumor.update(status='failed',stage='input_check',reason='Matched normal failed: '+normal.get('reason','unknown'))
        tumor['warnings']=tumor.get('warnings',[])+['Normal: '+x for x in normal.get('warnings',[])]
        results.append(tumor)
    return results,sources
