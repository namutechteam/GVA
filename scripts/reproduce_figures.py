"""Regenerate final figures without changing supplied exports or results."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, default=ROOT / 'output/figures')
a = p.parse_args()
out = a.output.resolve()
if out == (ROOT / 'figures').resolve():
    p.error('choose a new output directory to preserve supplied figures')
out.mkdir(parents=True, exist_ok=True)
env = dict(os.environ)
env.pop('GVA_PALETTE', None)
env['MPLBACKEND'] = 'Agg'
env.setdefault('MPLCONFIGDIR', str(ROOT / '.matplotlib'))
for stem in [f'make_fig{i}' for i in range(1, 6)] + [f'make_suppfig{i}' for i in range(1, 8)]:
    print(f'Rendering {stem}', flush=True)
    subprocess.run([sys.executable, str(ROOT / 'plotting' / (stem + '.py')),
                    str(ROOT / 'data/GVA_Source_Data.xlsx'), str(out)], check=True, env=env)
print(f'Completed 12 figures in {out}')
