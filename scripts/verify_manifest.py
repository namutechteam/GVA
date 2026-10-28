"""Verify the prepared payload against its SHA-256 manifest."""
import hashlib
from pathlib import Path
import sys
root = Path(__file__).resolve().parents[1]
errors = []
lines = (root / 'verification/files.sha256').read_text().splitlines()
for line in lines:
    digest, name = line.split('  ', 1)
    path = root / name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        errors.append(name)
for name in errors:
    print('FAILED:', name)
print(f'{len(lines)} files checked; {len(errors)} failures')
sys.exit(bool(errors))
