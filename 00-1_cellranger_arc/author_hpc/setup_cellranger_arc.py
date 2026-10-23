#!/usr/bin/env python3
"""Generate non-destructive raw-data symlinks and Cell Ranger ARC library CSVs."""
import argparse
import csv
import os
import re
from pathlib import Path

DEFAULT_PROJECT = Path('/gpfs/projects/DamaghiGroup/projects/mcf7_aa_cacr_lipid')
DEFAULT_RAW = Path('/gpfs/projects/DamaghiGroup/raw/sequencing/2023-10-04_novogene_mcf7_aa_cacr_snmult')


def ensure_symlink(link: Path, target: Path):
    if not target.is_dir():
        raise FileNotFoundError(f'Missing raw directory: {target}')
    if link.is_symlink():
        if link.resolve() != target.resolve():
            raise RuntimeError(f'Conflicting symlink: {link} -> {os.readlink(link)}')
        return
    if link.exists():
        raise RuntimeError(f'Path exists; refusing to replace: {link}')
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(target, target_is_directory=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--project', type=Path, default=DEFAULT_PROJECT)
    p.add_argument('--raw', type=Path, default=DEFAULT_RAW)
    args = p.parse_args()
    project = args.project.resolve()
    raw = args.raw.resolve()
    if not raw.is_dir():
        raise SystemExit(f'Raw directory not found: {raw}')
    samples_path = project / 'config' / 'samples.csv'
    with samples_path.open(newline='') as f:
        samples = list(csv.DictReader(f))
    if not samples:
        raise SystemExit('Empty sample manifest')
    for row in samples:
        sample = row['sample_id']
        if not sample or '/' in sample or sample in {'.', '..'}:
            raise SystemExit(f'Invalid sample ID: {sample}')
        for modality, key in [('scRNA', 'gex_id'), ('scATAC', 'atac_id')]:
            source = raw / modality / row[key]
            link = project / 'data' / 'raw' / modality / row[key]
            ensure_symlink(link, source)
            for read in (('R1', 'R2') if modality == 'scRNA' else ('R1', 'R2', 'R3')):
                if not list(source.glob(f'{row[key]}*_{read}_*.fastq.gz')):
                    raise RuntimeError(f'Missing {read} FASTQ under {source}')
        output = project / 'config' / 'libraries' / f'{sample}.csv'


        def fastq_sample_name(modality, sample_id):
            source = raw / modality / sample_id
            fastqs = sorted(source.glob(f'{sample_id}*_R1_*.fastq.gz'))
            if not fastqs:
                raise RuntimeError(f'No R1 FASTQ found: {source}')

            names = {
                re.sub(r'_S\d+_L\d+_R1_\d+\.fastq\.gz$', '', f.name)
                for f in fastqs
            }
            if len(names) != 1:
                raise RuntimeError(f'Ambiguous FASTQ sample names: {names}')
            return names.pop()


        entries = [
            [
                str((project / 'data/raw/scRNA' / row['gex_id']).absolute()),
                fastq_sample_name('scRNA', row['gex_id']),
                'Gene Expression',
            ],
            [
                str((project / 'data/raw/scATAC' / row['atac_id']).absolute()),
                fastq_sample_name('scATAC', row['atac_id']),
                'Chromatin Accessibility',
            ],
        ]


        output.parent.mkdir(parents=True, exist_ok=True)
        import io
        buf = io.StringIO()
        writer = csv.writer(buf, lineterminator='\n')
        writer.writerow(['fastqs', 'sample', 'library_type'])
        writer.writerows(entries)
        contents = buf.getvalue()
        if output.exists() and output.read_text() != contents:
            raise RuntimeError(f'Libraries CSV exists with different content: {output}; review before replacing')
        if not output.exists():
            output.write_text(contents)
        print(f'Validated {sample}: {row["gex_id"]} + {row["atac_id"]} -> {output}')
    print('Setup complete. Raw contents were not modified.')


if __name__ == '__main__':
    main()
