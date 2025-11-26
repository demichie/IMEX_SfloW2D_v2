#!/usr/bin/env python3
# Author: ECP_BREARD
# Cleaner: keep IMEX* solver (NOT *.inp), all *.py, all *.template, README*; delete everything else.
import argparse, fnmatch, os, shutil, sys
from pathlib import Path

def is_solver_to_keep(p: Path) -> bool:
    name = p.name
    if fnmatch.fnmatch(name, 'IMEX*'):
        if name.lower().endswith('.inp'):
            return False
        return True
    return False

def should_keep(p: Path) -> bool:
    name = p.name
    if name.startswith('.'):
        return True
    if is_solver_to_keep(p):
        return True
    if name.lower().startswith('readme'):
        return True
    if name.endswith('.py') or name.endswith('.pyw'):
        return True
    if name.endswith('.template'):
        return True
    return False

def main():
    ap = argparse.ArgumentParser(description='Delete everything in CWD except IMEX* solver (not *.inp), *.py, *.template, README*.')
    ap.add_argument('--apply', action='store_true', help='Actually delete (default: dry-run).')
    ap.add_argument('--limit', type=int, default=100, help='Max items to list (default: 100).')
    ap.add_argument('--include-hidden', action='store_true', help='Also consider hidden files/dirs for deletion.')
    ap.add_argument('--yes', '-y', action='store_true', help='Skip confirmation prompt.')
    args = ap.parse_args()

    cwd = Path.cwd()
    keep, delete = [], []
    for p in sorted(cwd.iterdir(), key=lambda x: x.name.lower()):
        if p.name == __file__:
            keep.append(p); continue
        if (not args.include_hidden) and p.name.startswith('.'):
            keep.append(p); continue
        if should_keep(p):
            keep.append(p); continue
        delete.append(p)

    print('Will KEEP:')
    for k in keep:
        print('  K', k.name + ('/' if k.is_dir() else ''))
    print('\nWill DELETE ({} total, showing up to {}):'.format(len(delete), args.limit))
    for d in delete[:max(0,args.limit)]:
        print('  D', d.name + ('/' if d.is_dir() else ''))
    if len(delete) > args.limit:
        print('  ... ({} more not shown)'.format(len(delete)-args.limit))

    if not args.apply or not delete:
        print('\n[dry-run] Nothing deleted. Use --apply to actually delete.')
        return

    if not args.yes:
        resp = input('\n*** WARNING *** This will permanently delete {} items.\nType DELETE to confirm: '.format(len(delete))).strip()
        if resp != 'DELETE':
            print('Aborted.'); return

    for d in delete:
        try:
            if d.is_dir():
                shutil.rmtree(d)
            else:
                d.unlink()
        except Exception as e:
            print(f'  ! Failed to delete {d}: {e}', file=sys.stderr)
    print('Done.')

if __name__ == '__main__':
    main()
