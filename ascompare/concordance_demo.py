"""Explicitly synthetic, BAM-free end-to-end AmpliconSuite Cross-Version Concordance demo."""
from pathlib import Path

from .common import atomic_json
from .versions import load_lock


SOURCE_VERSION='suite-1.6.0'
REFERENCE_VERSION='suite-1.6.2'


def _write_result(directory, sample, amplicons, entry):
    directory=Path(directory);aa=directory/'aa';ac=directory/'ac'
    aa.mkdir(parents=True);ac.mkdir()
    beds=ac/'classification_classification_bed_files';beds.mkdir()
    support=entry.get('supported_labels',entry.get('supported_classes',[]))
    flags=[name for name in ('ecDNA','BFB','FAN') if name in support]
    header=['sample_name','amplicon_number','amplicon_decomposition_class']+[name+'+' for name in flags]
    profiles=['\t'.join(header)]
    (aa/f'{sample}_summary.txt').write_text(f'#Amplicons = {len(amplicons)}\n')
    for number,amplicon in enumerate(amplicons,1):
        intervals=amplicon['intervals'];labels=amplicon['labels'];name=f'{sample}_amplicon{number}'
        if not set(labels)<=set(support):raise ValueError('Demo fixture requests unsupported locked category')
        graph=[];segments=[]
        for sid,(chrom,start,end) in enumerate(intervals,1):
            graph.append(f'sequence\t{chrom}:{start}-\t{chrom}:{end-1}+\t8\t20\t{end-start}\t10\n')
            segments.append(f'Segment\t{sid}\t{chrom}\t{start}\t{end-1}\n')
        kind='Cyclic' if any(x in labels for x in ('ecDNA','BFB','FAN')) else labels[0]
        traversals=[f'{sid}+' for sid in range(1,len(intervals)+1)]
        if kind=='Linear':traversals=['0+',*traversals,'0-']
        (aa/(name+'_graph.txt')).write_text(''.join(graph))
        (aa/(name+'_cycles.txt')).write_text(''.join(segments)+'Cycle=1;Copy_count=8;Segments='+','.join(traversals)+'\n')
        profiles.append('\t'.join([sample,f'amplicon{number}',kind]+['Positive' if label in labels else 'None detected' for label in flags]))
        for label in labels:
            (beds/f'{name}_{label}_1_intervals.bed').write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in intervals))
    (ac/'classification_amplicon_classification_profiles.tsv').write_text('\n'.join(profiles)+'\n')
    atomic_json(directory/'stages.json',{'aa':{'status':'complete','evidence':'SYNTHETIC schema fixture; AA was not executed'},
        'ac':{'status':'complete','evidence':'SYNTHETIC schema fixture; AC was not executed'}})
    atomic_json(directory/'SYNTHETIC_DATA.json',{'synthetic':True,'sample_id':sample,
        'declared_fixture_components':{name:entry[name]['version'] for name in ('aa','ac')},
        'warning':'Invented result artifacts for software validation. No patient, BAM, AA or AC execution.'})


def write_demo_manifest(directory):
    """Create new synthetic inputs and return a portable relative-path manifest."""
    root=Path(directory).absolute()
    if root.exists():raise ValueError('Synthetic input directory must be new; existing inputs are never overwritten')
    root.mkdir(parents=True)
    entries={entry['id']:entry for entry in load_lock()['versions']}
    cases=[
        ('SYNTHETIC_01',
         [dict(labels=['ecDNA'],intervals=[['chr1',100000,180000]])],
         [dict(labels=['ecDNA'],intervals=[['chr1',100000,180000]])]),
        ('SYNTHETIC_02',
         [dict(labels=['ecDNA'],intervals=[['chr2',200000,290000]])],
         [dict(labels=['ecDNA','BFB'],intervals=[['chr2',200000,300000],['chr7',1000000,1020000]])]),
        ('SYNTHETIC_03',
         [dict(labels=['Linear'],intervals=[['chr3',400000,470000]]),
          dict(labels=['ecDNA'],intervals=[['chr8',2000000,2080000],['chr8',2090000,2130000]])],
         [dict(labels=['Linear'],intervals=[['chr3',400000,470000]]),
          dict(labels=['ecDNA','FAN'],intervals=[['chr8',2000000,2080000],['chr8',2090000,2130000]])])]
    samples=[]
    for sample,source,reference in cases:
        for version,amps in ((SOURCE_VERSION,source),(REFERENCE_VERSION,reference)):
            _write_result(root/'results'/sample/version,sample,amps,entries[version])
        seed_path=root/'seeds'/(sample+'.bed');seed_path.parent.mkdir(exist_ok=True)
        seed_path.write_text(''.join(f'{amp["intervals"][0][0]}\t{amp["intervals"][0][1]}\t{amp["intervals"][0][1]+1000}\tSYNTHETIC_CANDIDATE_{n}\n'
                                      for n,amp in enumerate(source,1)))
        def spec(version):
            return dict(run_dir=str(Path('results')/sample/version),
                        component_versions={name:entries[version][name]['version'] for name in ('aa','ac')})
        samples.append(dict(sample_id=sample,source_version=SOURCE_VERSION,reference_build='hg19',
            seeds_bed=str(seed_path.relative_to(root)),source_result=spec(SOURCE_VERSION),
            comparison_results=[dict(version=REFERENCE_VERSION,**spec(REFERENCE_VERSION))]))
    manifest=root/'manifest.json'
    atomic_json(manifest,dict(program='AmpliconSuite Cross-Version Concordance',mode='results_only',synthetic=True,
        reference_version=REFERENCE_VERSION,target_versions=[REFERENCE_VERSION,SOURCE_VERSION],
        reference={'build':'hg19'},samples=samples,runtimes=[],threads=1,
        statistics={'min_samples':20,'alpha':0.05},
        demo_scope='SYNTHETIC existing AA+AC artifacts only; no BAM and no executable runtimes. AC-only target remains missing.'))
    (root/'README.md').write_text('# SYNTHETIC TEST DATA\n\n'
        'These invented AA/AC artifacts exercise strict import, paired comparison and reporting. '
        'No patient input, BAM, AA reconstruction or AC classification was executed. '
        'SYNTHETIC_03 contains two complete amplicons; SYNTHETIC_02 retains a remote connected region.\n\n'
        'The manifest declares fixture Suite1.6.0 and Suite1.6.2 component combinations from the lock; '
        'these declarations are not evidence that those tools generated the invented outputs. '
        'The fixed report reference is suite-1.6.2. AC-only target results are unavailable without runtimes. '
        'Seeds are candidate regions, not confirmed amplifications.\n')
    return manifest


def demo(out):
    """Run the same public loader/orchestrator used by `asconcordance run`."""
    from .concordance import load_manifest,run_concordance
    out=Path(out).absolute()
    if out.exists():raise ValueError('Demo output must be new; use a different output path')
    inputs=out.parent/(out.name+'_synthetic_inputs')
    manifest=write_demo_manifest(inputs)
    config,lock,entries=load_manifest(manifest)
    result=run_concordance(config,lock,entries,out)
    atomic_json(out/'demo_provenance.json',{'synthetic':True,'manifest':str(manifest.resolve()),
        'input_directory':str(inputs.resolve()),'reference_version':REFERENCE_VERSION,
        'patient_data_used':False,'bam_used':False,'AA_executed':False,'AC_executed':False})
    return result
