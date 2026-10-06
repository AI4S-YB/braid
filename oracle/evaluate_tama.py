from pathlib import Path
import subprocess, random, json, difflib, sys
import argparse, hashlib, os, shutil
parser = argparse.ArgumentParser(description='TAMA differential matrix; ignores console logs, compares all successful output files')
parser.add_argument('--out', required=True)
parser.add_argument('--upstream', type=Path, default=Path(__file__).resolve().parents[2] / 'tama')
parser.add_argument('--python2', help='CPython 2.7 executable; otherwise use the documented compatibility adapter')
parser.add_argument('--merge-cases', type=int, default=120)
parser.add_argument('--collapse-cases', type=int, default=80)
options = parser.parse_args()
repo = Path(__file__).resolve().parents[1]
root = Path(options.out).resolve()
root.mkdir(parents=True, exist_ok=False)
rng = random.Random(20261004)
results = []
metadata = {'oracle': options.python2 or 'python2_compat.py (fallback)', 'sha256': {kind: hashlib.sha256((options.upstream / ('tama_' + kind + '.py')).read_bytes()).hexdigest() for kind in ['merge', 'collapse']}, 'binary_sha256': hashlib.sha256((repo / 'target/debug/plena').read_bytes()).hexdigest(), 'seed': 20261004}
(root / 'metadata.json').write_text(json.dumps(metadata, indent=2))

def compare(name, kind, args):
    d = root / name
    d.mkdir(exist_ok=True)
    success = {}
    for side in ['py', 'rs']:
        out = d / side
        out.mkdir(exist_ok=True)
        oracle = [options.python2] if options.python2 else [sys.executable, str(repo / 'oracle/python2_compat.py')]
        cmd = oracle + [str(options.upstream / ('tama_' + kind + '.py'))] if side == 'py' else [str(repo / 'target/debug/plena'), 'tama', kind]
        cmd += args + ['-p', str(out / 'out')]
        (out / 'command.json').write_text(json.dumps(cmd))
        try:
            p = subprocess.run(cmd, capture_output=True, timeout=60, env={**os.environ, 'PYTHONHASHSEED': '0'})
        except subprocess.TimeoutExpired:
            success[side] = False
            continue
        (out / 'stdout').write_bytes(p.stdout)
        (out / 'stderr').write_bytes(p.stderr)
        success[side] = p.returncode == 0
        if side == 'py':
            success[side] = success[side] and (b'TAMA Merge has completed successfully!' in p.stdout if kind == 'merge' else (out / 'out_report.txt').exists() and 'has run successfully' in (out / 'out_report.txt').read_text())
    files = {side: {f.name: f.read_text().replace(str(d / side), '<prefix>') for f in (d / side).glob('out*') if f.name != 'stdout'} for side in ['py', 'rs']}
    diff = [k for k in files['py'].keys() | files['rs'].keys() if files['py'].get(k) != files['rs'].get(k)]
    status = 'match' if all(success.values()) and (not diff) else 'both_failed' if not any(success.values()) else 'mismatch'
    results.append(dict(name=name, kind=kind, status=status, success=success, diff=sorted(diff)))
    if status == 'mismatch':
        print(name, success, diff, flush=True)
        if all(success.values()):
            for key in diff[:1]:
                print(''.join(list(difflib.unified_diff(files['py'].get(key, '').splitlines(True), files['rs'].get(key, '').splitlines(True)))[:14]), flush=True)
    if len(results) % 10 == 0:
        print('processed', len(results), flush=True)
    (root / 'results.json').write_text(json.dumps(results, indent=2))

def bed(name, starts, ends, strand):
    return '\t'.join(map(str, ['chr1', starts[0], ends[-1], name, 40, strand, starts[0] + 5, ends[-1] - 5, 0, len(starts), ','.join((str(b - a) for a, b in zip(starts, ends))), ','.join((str(a - starts[0]) for a in starts))])) + '\n'
for seed in range(options.merge_cases):
    d = root / f'merge-{seed:03}'
    d.mkdir(exist_ok=True)
    lines = []
    caps = ['capped'] * 3 if seed % 3 == 0 else ['no_cap'] * 3 if seed % 3 == 1 else ['capped', 'no_cap', 'no_cap']
    for sample in range(3):
        beds = []
        for t in range(16):
            locus = rng.randrange(3)
            n = rng.choice([1, 2, 3])
            strand = rng.choice(['+', '-'])
            base = locus * 3000 + 100
            coords = [base + x * 200 for x in range(3)]
            starts = [x + rng.choice([0, 0, 2, 10, 20]) for x in coords]
            ends = [x + 100 + rng.choice([0, 0, 2, 10, 20]) for x in coords]
            if strand == '+':
                starts = starts[-n:]
                ends = ends[-n:]
            else:
                starts = starts[:n]
                ends = ends[:n]
            beds.append(bed(f'G{locus};G{locus}.{sample}_{t}', starts, ends, strand))
        rng.shuffle(beds)
        path = d / f's{sample}.bed'
        path.write_text(''.join(beds))
        lines.append(f'{path}\t{caps[sample]}\t{rng.choice(['1,1,1', '2,1,3', '1,3,2'])}\tS{sample}')
    manifest = d / 'files.tsv'
    manifest.write_text('\n'.join(lines) + '\n')
    flags = ['-f', str(manifest), '-a', str(rng.choice([0, 10, 20, 50])), '-m', str(rng.choice([0, 10, 20])), '-z', str(rng.choice([0, 10, 20])), '-d', 'merge_dup', '-e', rng.choice(['common_ends', 'longest_ends'])]
    if seed % 2 == 0:
        flags += ['-s', 'S0,S2', '-cds', 'S1,S2']
    compare(d.name, 'merge', flags)
gmap = repo / 'tests/parity/gmap_collapse'
bam = root / 'gmap.bam'
import shutil
ref = root / 'genome.fa'
shutil.copyfile(gmap / 'test_genome.fa', ref)
p = subprocess.run(['samtools', 'view', '-b', '-T', str(ref), '-o', str(bam), str(gmap / 'gmap_test.sam')], capture_output=True)
assert p.returncode == 0, p.stderr
for cap in ['capped', 'no_cap']:
    for mode in ['original', 'low_mem']:
        for fmt in ['SAM', 'BAM']:
            for variant in range(4):
                inp = gmap / 'gmap_test.sam' if fmt == 'SAM' else bam
                flags = ['-s', str(inp), '-f', str(gmap / 'test_genome.fa'), '-x', cap, '-rm', mode, '-b', fmt, '-log', 'log_off', '-vc', '3']
                if variant == 1:
                    flags += ['-e', 'longest_ends', '-sj', 'sj_priority']
                if variant == 2:
                    flags += ['-a', '0', '-m', '0', '-z', '0', '-icm', 'ident_map']
                if variant == 3:
                    flags += ['-c', '80', '-i', '80', '-lde', '3', '-sjt', '5']
                compare(f'gmap-{cap}-{mode}-{fmt}-{variant}', 'collapse', flags)
rng = random.Random(20261004)
for seed in range(options.collapse_cases):
    d = root / f'random-collapse-{seed:03}'
    d.mkdir(exist_ok=True)
    genome = d / 'genome.fa'
    genome.write_text('>chr1\n' + 'A' * 16000 + '\n')
    rows = []
    for t in range(50):
        base = rng.randrange(4) * 3000 + 100
        n = rng.choice([1, 2, 3])
        strand = rng.choice([0, 16])
        starts = [base + x * 250 + rng.choice([0, 1, 5, 10, 20]) for x in range(3)]
        ends = [base + x * 250 + 100 + rng.choice([0, 1, 5, 10, 20]) for x in range(3)]
        if strand == 0:
            starts = starts[-n:]
            ends = ends[-n:]
        else:
            starts = starts[:n]
            ends = ends[:n]
        cigar = ''
        length = 0
        for i, (a, b) in enumerate(zip(starts, ends)):
            if i:
                cigar += str(a - ends[i - 1]) + 'N'
            cigar += str(b - a) + 'M'
            length += b - a
        seq = list('A' * length)
        if t % 7 == 0:
            seq[50] = 'C'
        if t % 11 == 0:
            cigar = '3S' + cigar
            seq = list('TTT') + seq
        name = f'r{t}'
        if seed % 5 == 0 and t % 7 == 0:
            name = 'multi' + str(t % 3)
        rows.append((starts[0], f'{name}\t{strand}\tchr1\t{starts[0]}\t40\t{cigar}\t*\t0\t0\t{''.join(seq)}\t*\n'))
    rows.sort(key=lambda x: x[0])
    sam = d / 'reads.sam'
    sam.write_text(''.join((x[1] for x in rows)))
    flags = ['-s', str(sam), '-f', str(genome), '-x', ['capped', 'no_cap'][seed % 2], '-rm', ['original', 'low_mem'][seed // 2 % 2], '-log', 'log_off', '-c', '80', '-i', '80', '-d', 'merge_dup']
    if seed % 3 == 0:
        flags += ['-sj', 'sj_priority', '-e', 'longest_ends']
    if seed % 3 == 1:
        flags += ['-a', '0', '-m', '0', '-z', '0', '-icm', 'ident_map']
    if seed % 3 == 2:
        flags += ['-lde', '1', '-sjt', '5']
    compare(d.name, 'collapse', flags)
print({s: sum((x['status'] == s for x in results)) for s in ['match', 'mismatch', 'both_failed']})
sys.exit(int(any((x['status'] == 'mismatch' for x in results))))
