"""
Build a Cesium-ready height tile pyramid from the GEBCO_2026 sub-ice topography GeoTIFFs.

Input : gebco_2026_sub_ice_topo_geotiff/*.tif  (8 files, 90x90 deg each, 21600x21600 int16 m, 15")
Output: tiles/gebco_topo/{z}/{x}/{y}.png   (Cesium GeographicTilingScheme, y counted from north)
        tiles/gebco_topo/meta.json

Each tile is a 257x257 grid of heights at the tile's cell corners (edges shared with neighbours),
so the same tiles serve both as terrain (every 4th sample -> 65x65 heightmap) and as source for
the depth-coloured imagery rendered in the browser.

Height encoding (lossless for integer metres, readable via a browser canvas):
    v = height_m + 32768 ;  R = v >> 8 ;  G = v & 255 ;  B = 0

Level 7 (~20" grid spacing) is the finest level, close to the native 15"; Cesium upsamples beyond it.
Heights are bilinearly interpolated from the pixel-centre registered source; lower levels are
[1 2 1] smoothed and decimated by 2.

Usage: .venv/Scripts/python.exe build_topo_tiles.py
"""
import json
import os
import struct
import time
from multiprocessing import Pool

import numpy as np
from PIL import Image

SRC_DIR = "gebco_2026_sub_ice_topo_geotiff"
OUT_DIR = os.path.join("tiles", "gebco_topo")
MAX_LEVEL = 7
TILE = 256
NATIVE = 21600          # cells per file side (90 deg / 15")
CELLS_PER_DEG = 240
OFFSET = 32768
BAND = 1024             # rows per chunk in the resampling passes, bounds peak memory


def src_path(fx, fy):
    """fx: 0..3 west->east, fy: 0 = northern hemisphere, 1 = southern."""
    w = -180 + 90 * fx
    n, s = (90.0, 0.0) if fy == 0 else (0.0, -90.0)
    return os.path.join(SRC_DIR, f"gebco_2026_sub_ice_n{n}_s{s}_w{w:.1f}_e{w + 90:.1f}_geotiff.tif")


def open_tiff(path):
    """Memory-map an uncompressed, single-band int16, one-row-per-strip GeoTIFF."""
    with open(path, "rb") as f:
        head = f.read(8)
        bo = "<" if head[:2] == b"II" else ">"
        f.seek(struct.unpack(bo + "I", head[4:8])[0])
        tags = {}
        for _ in range(struct.unpack(bo + "H", f.read(2))[0]):
            tag, typ, cnt, val = struct.unpack(bo + "HHI4s", f.read(12))
            tags[tag] = (typ, cnt, val)
        short = lambda t: struct.unpack(bo + "H", tags[t][2][:2])[0]
        assert (short(256), short(257), short(258), short(259), short(339)) == (NATIVE, NATIVE, 16, 1, 2), path
        typ, cnt, val = tags[273]                     # StripOffsets
        f.seek(struct.unpack(bo + "I", val)[0])
        offsets = np.frombuffer(f.read(cnt * 4), dtype=bo + "u4")
        assert np.all(np.diff(offsets) == NATIVE * 2), "strips are not contiguous"
    return np.memmap(path, dtype=bo + "i2", mode="r", offset=int(offsets[0]), shape=(NATIVE, NATIVE))


def padded_source(fx, fy):
    """Source file with a 1-cell border taken from its neighbours (longitude wraps, poles clamp),
    so interpolated values match exactly along file boundaries."""
    src = lambda x, y: open_tiff(src_path(x % 4, y))
    p = np.empty((NATIVE + 2, NATIVE + 2), np.int16)
    p[1:-1, 1:-1] = src(fx, fy)
    p[1:-1, 0] = src(fx - 1, fy)[:, -1]
    p[1:-1, -1] = src(fx + 1, fy)[:, 0]
    p[0, 1:-1] = src(fx, 0)[-1] if fy == 1 else p[1, 1:-1]      # row above: northern file / pole
    p[-1, 1:-1] = src(fx, 1)[0] if fy == 0 else p[-2, 1:-1]     # row below: southern file / pole
    p[0, 0], p[0, -1], p[-1, 0], p[-1, -1] = p[0, 1], p[0, -2], p[-1, 1], p[-1, -2]
    return p


def lerp_weights(size):
    """Corner sample positions (0..size-1 spanning the file) in padded-source index space."""
    u = np.arange(size) * (NATIVE / (size - 1)) + 0.5      # padded index of each corner point
    i0 = np.minimum(np.floor(u).astype(np.int64), NATIVE)
    return i0, (u - i0).astype(np.float32)


def resample(p, size):
    """Separable bilinear resample of the padded source to a size x size corner grid."""
    i0, w = lerp_weights(size)
    cols = np.empty((p.shape[0], size), np.float32)
    for r in range(0, p.shape[0], BAND):
        blk = p[r:r + BAND].astype(np.float32)
        cols[r:r + BAND] = blk[:, i0] * (1 - w) + blk[:, i0 + 1] * w
    out = np.empty((size, size), np.float32)
    for r in range(0, size, BAND):
        a, b, ww = i0[r:r + BAND], i0[r:r + BAND] + 1, w[r:r + BAND, None]
        out[r:r + BAND] = cols[a] * (1 - ww) + cols[b] * ww
    return out


def smooth_decimate(a):
    """[1 2 1]/4 filter in both axes (edge-replicated), then keep every 2nd corner sample."""
    e = np.pad(a, 1, mode="edge")
    r = (e[:, :-2] + 2 * e[:, 1:-1] + e[:, 2:]) * 0.25
    c = (r[:-2] + 2 * r[1:-1] + r[2:]) * 0.25
    return c[::2, ::2].copy()


def save_tile(a, level, x, y):
    v = np.clip(np.rint(a) + OFFSET, 0, 65535).astype(np.uint16)
    rgb = np.dstack([(v >> 8).astype(np.uint8), (v & 255).astype(np.uint8), np.zeros_like(v, np.uint8)])
    d = os.path.join(OUT_DIR, str(level), str(x))
    os.makedirs(d, exist_ok=True)
    Image.fromarray(rgb, "RGB").save(os.path.join(d, f"{y}.png"))


def process_file(args):
    """Tile levels MAX_LEVEL..1 for one 90x90 deg source file."""
    fx, fy = args
    t0 = time.time()
    p = padded_source(fx, fy)
    hmin, hmax = int(p.min()), int(p.max())
    n = 2 ** (MAX_LEVEL - 1)
    arr = resample(p, n * TILE + 1)
    del p

    count = 0
    for level in range(MAX_LEVEL, 0, -1):
        n = 2 ** (level - 1)
        for j in range(n):
            for i in range(n):
                save_tile(arr[j * TILE:j * TILE + TILE + 1, i * TILE:i * TILE + TILE + 1],
                          level, fx * n + i, fy * n + j)
                count += 1
        if level > 1:
            arr = smooth_decimate(arr)
    print(f"  {os.path.basename(src_path(fx, fy))}: {count} tiles, {hmin}..{hmax} m, {time.time() - t0:.0f} s", flush=True)
    return (fx, fy), arr, hmin, hmax


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    jobs = [(fx, fy) for fy in range(2) for fx in range(4)]
    with Pool(len(jobs)) as pool:
        results = pool.map(process_file, jobs)
    level1 = {key: arr for key, arr, _, _ in results}

    # Level 0: two 180x180 deg tiles, each stitched from 2x2 level-1 tiles (shared edges dropped).
    for x in range(2):
        row = lambda fy: np.hstack([level1[(2 * x, fy)], level1[(2 * x + 1, fy)][:, 1:]])
        save_tile(smooth_decimate(np.vstack([row(0), row(1)[1:]])), 0, x, 0)

    meta = {
        "name": "GEBCO_2026 Grid (sub-ice topography)",
        "attribution": "GEBCO Compilation Group (2026) GEBCO 2026 Grid (doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa)",
        "maxLevel": MAX_LEVEL,
        "tileSize": TILE,
        "samples": TILE + 1,
        "encoding": "height_m = R*256 + G - 32768",
        "offset": OFFSET,
        "minHeight": min(r[2] for r in results),
        "maxHeight": max(r[3] for r in results),
    }
    with open(os.path.join(OUT_DIR, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f)
    print(f"Done: {meta['minHeight']}..{meta['maxHeight']} m, {time.time() - t0:.0f} s total")


if __name__ == "__main__":
    main()
