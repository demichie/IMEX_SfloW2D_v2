#!/usr/bin/env python3
# Author: ECP_BREARD
# height_summary_nc.py — centreline height plots + MP4/GIF animation (1:1 axis scale)

import argparse
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt


def parse_inp(inp_path: Path):
    txt = inp_path.read_text(encoding='utf-8', errors='ignore')

    def ffloat(key, default):
        m = re.search(rf"{key}\s*=\s*([0-9\.\+EeDd-]+)", txt)
        return float(m.group(1).replace('D', 'E')) if m else default

    def fstr(key, default):
        m = re.search(rf'{key}\s*=\s*"(.*?)"', txt)
        return m.group(1) if m else default

    pres = ffloat('PRES', 101300.0)
    t0   = ffloat('T_START', 0.0)
    dt   = ffloat('DT_OUTPUT', 1.0)
    run  = fstr('RUN_NAME', 'run')
    return pres, t0, dt, run


def nc_read_height(nc_path: Path, var_hint: str = None):
    """
    Return (H, time, x, y) where H has shape (nt, ny, nx).
    Tries netCDF4 first, then h5py.
    """
    # Try netCDF4
    try:
        from netCDF4 import Dataset
        ds = Dataset(nc_path.as_posix(), 'r')

        candidates = []
        if var_hint and var_hint in ds.variables:
            candidates = [var_hint]
        if not candidates:
            for name in ds.variables:
                low = name.lower()
                if low in ('h', 'height', 'thickness', 'flow_height', 'hflow', 'sum_h'):
                    candidates.append(name)
        if not candidates:
            for name, v in ds.variables.items():
                dims = getattr(v, 'dimensions', ())
                if len(dims) == 3 and any('time' in d.lower() for d in dims):
                    candidates = [name]
                    break
        if not candidates:
            raise RuntimeError('No suitable height variable found in NetCDF.')

        v = ds.variables[candidates[0]]
        H = np.array(v[:], dtype=float)

        x = y = t = None
        for ax in ('x', 'xc', 'lon', 'longitude'):
            if ax in ds.variables:
                x = np.array(ds.variables[ax][:], dtype=float).reshape(-1)
                break
        for ay in ('y', 'yc', 'lat', 'latitude'):
            if ay in ds.variables:
                y = np.array(ds.variables[ay][:], dtype=float).reshape(-1)
                break
        if 'time' in ds.variables:
            t = np.array(ds.variables['time'][:], dtype=float).reshape(-1)

        ds.close()
        if H.ndim == 2:
            H = H[None, ...]
        return H, t, x, y
    except Exception:
        pass

    # Try h5py
    try:
        import h5py
        f = h5py.File(nc_path.as_posix(), 'r')

        def first_3d(h):
            for k, v in h.items():
                if isinstance(v, h5py.Dataset) and v.ndim == 3:
                    return v
                if isinstance(v, h5py.Group):
                    g = first_3d(v)
                    if g is not None:
                        return g
            return None

        d = first_3d(f)
        if d is None:
            raise RuntimeError('No 3D dataset found (h5py).')

        H = np.array(d[...], dtype=float)
        t = f['time'][...] if 'time' in f.keys() else None
        x = f['x'][...] if 'x' in f.keys() else None
        y = f['y'][...] if 'y' in f.keys() else None
        f.close()

        if H.ndim == 2:
            H = H[None, ...]
        return H, (None if t is None else np.array(t).reshape(-1)), x, y
    except Exception as e:
        raise SystemExit(f'Failed to read height from file: {e}')


def centerline_height(H, x=None, center_fraction=0.5):
    """Extract height along y = center_fraction*(Ny-1)."""
    nt, ny, nx = H.shape
    j = int(round((ny - 1) * float(center_fraction)))
    j = max(0, min(ny - 1, j))
    xs = np.arange(nx) if x is None else np.array(x).reshape(-1)
    line = H[:, j, :]  # (nt, nx)
    return xs, line


def main():
    p = argparse.ArgumentParser(description='Centreline height plots + movie for dam-break runs (1:1 axis scale).')
    p.add_argument('nc', help='NetCDF file, e.g. dambreak2D.nc')
    p.add_argument('--center_fraction', type=float, default=0.5, help='Centreline position in [0,1] (default 0.5).')
    p.add_argument('--fps', type=int, default=12, help='Frames per second for movie.')
    p.add_argument('--animate', action='store_true', help='Write MP4 if ffmpeg available, else GIF.')
    args = p.parse_args()

    nc_path = Path(args.nc)
    if not nc_path.exists():
        raise SystemExit(f'NC file not found: {nc_path}')
    inp = Path('IMEX_SfloW2D.inp')
    if not inp.exists():
        raise SystemExit('IMEX_SfloW2D.inp not found in current directory.')

    _, t0, dt, run = parse_inp(inp)

    H, time, x, _ = nc_read_height(nc_path)
    nt, ny, nx = H.shape
    if time is None or len(time) != nt:
        time = t0 + np.arange(nt) * dt

    xs, line = centerline_height(H, x, args.center_fraction)
    xmin, xmax = float(xs.min()), float(xs.max())

    outdir = Path('frames')
    outdir.mkdir(parents=True, exist_ok=True)

    # Use a stable vertical scale (based on all frames) but still show 1:1 aspect
    ymax_all = max(1e-6, float(np.nanmax(line)) * 1.1)

    for it in range(nt):
        fig, ax = plt.subplots(figsize=(8, 3.6))
        ax.set_title(f'{run} — t = {time[it]:.3f} s  (centreline y={args.center_fraction:.2f})')

        # Fill under current height (darker grey, mild transparency, no outline)
        ax.fill_between(
            xs, 0.0, line[it, :],
            where=(line[it, :] > 0.0),
            color='0.6',
            alpha=0.45,
            step='mid',
            linewidth=0,
            zorder=1
        )

        # Current free-surface line (blue)
        ax.plot(xs, line[it, :], color='C0', lw=2.0, label='height(t)', zorder=3)

        # Domain bounds (vertical lines at x-min and x-max)
        ax.axvline(xmin, color='0.8', lw=1.0, ls='--', zorder=0)
        ax.axvline(xmax, color='0.8', lw=1.0, ls='--', zorder=0)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0.0, ymax_all)

        # Enforce 1:1 units in data space
        ax.set_aspect('equal', adjustable='box')

        ax.set_xlabel('x (m)')
        ax.set_ylabel('flow height h (m)')
        ax.grid(alpha=0.25)
        ax.legend(loc='upper right')
        fig.tight_layout()

        fp = outdir / f'{nc_path.stem}_hcenter_{it:04d}.png'
        fig.savefig(fp, dpi=120)
        plt.close(fig)

    print(f'Generated {nt} frames in {outdir}')

    if args.animate:
        import imageio.v2 as imageio
        frames = sorted(outdir.glob(f'{nc_path.stem}_hcenter_*.png'))
        mp4_path = Path(f'{nc_path.stem}_height_movie.mp4')
        try:
            import imageio_ffmpeg  # ensure writer available
            with imageio.get_writer(str(mp4_path), fps=args.fps, format='FFMPEG') as writer:
                for fp in frames:
                    writer.append_data(imageio.imread(fp))
            print(f'Saved MP4: {mp4_path}')
        except Exception as e:
            gif_path = Path(f'{nc_path.stem}_height_movie.gif')
            imgs = [imageio.imread(fp) for fp in frames]
            imageio.mimsave(gif_path, imgs, duration=1.0 / float(args.fps))
            print(f'FFMPEG not available ({e}). Wrote GIF instead: {gif_path}')
            print("To enable MP4: pip install -U 'imageio[ffmpeg]' imageio-ffmpeg  (and system ffmpeg).")


if __name__ == '__main__':
    main()
