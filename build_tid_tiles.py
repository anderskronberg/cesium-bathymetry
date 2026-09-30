"""
Build a Cesium-ready PNG tile pyramid from the GEBCO_2026 TID (Type Identifier) GeoTIFFs.

Input : gebco_2026_tid_geotiff/*.tif  (8 files, 90x90 deg each, 21600x21600 uint8, 15 arc-second)
Output: tiles/gebco_tid/{z}/{x}/{y}.png   (Cesium GeographicTilingScheme, y counted from north)
        tiles/gebco_tid/meta.json         (class palette/labels + tile availability bitmasks)

Level 8 (~0.0055 deg/px) slightly oversamples the native 15" grid, so no cells are lost.
Lower levels use a block majority vote over the level-8 array (ties broken in favour of direct
measurements), so each zoomed-out pixel shows the dominant source type in its area.

Land (TID 0) and nodata (127) are written transparent; tiles containing only those are skipped.

Usage: .venv/Scripts/python.exe build_tid_tiles.py
"""
import base64
import json
import os
import struct
import time
from multiprocessing import Pool

import numpy as np
from PIL import Image

SRC_DIR = "gebco_2026_tid_geotiff"
OUT_DIR = os.path.join("tiles", "gebco_tid")
MAX_LEVEL = 8
TILE = 256
NATIVE = 21600          # cells per file side (90 deg / 15")
NODATA = 127

# TID classes: code -> (label, hex colour). Source: GEBCO grid documentation, Table 2.
CLASSES = {
    0:  ("Land", None),
    10: ("Singlebeam echo-sounder", "#e41a1c"),
    11: ("Multibeam echo-sounder", "#ff7f00"),
    12: ("Seismic", "#a65628"),
    13: ("Isolated sounding", "#f781bf"),
    14: ("ENC sounding", "#984ea3"),
    15: ("Lidar", "#ffff33"),
    16: ("Optical light sensor", "#00e5a0"),
    17: ("Combination of direct methods", "#b15928"),
    40: ("Predicted – satellite gravity", "#1f3f73"),
    41: ("Interpolated – computer algorithm", "#4a90c2"),
    42: ("Contours from charts", "#66c2a5"),
    43: ("Contours from ENCs", "#1b9e77"),
    44: ("Multi-source, gravity-guided grid", "#7570b3"),
    45: ("Predicted – airborne gravity", "#8dd3c7"),
    46: ("Grounded iceberg draft", "#bebada"),
    47: ("Grounded Argo float", "#80b1d3"),
    48: ("Animal-borne data logger", "#b3de69"),
    70: ("Pre-generated grid (mixed)", "#9e9e9e"),
    71: ("Unknown source", "#616161"),
    72: ("Steering points", "#d0d0d0"),
}
GROUPS = [("Direct measurements", range(10, 18)),
          ("Indirect measurements", range(40, 49)),
          ("Unknown", range(70, 73))]

# Tie-break priority for downsampling (first = highest). Codes not listed rank after these, before land/nodata.
PRIORITY = [11, 15, 10, 16, 17, 12, 14, 13, 70, 72, 71, 44, 42, 43, 41, 45, 46, 47, 48, 40]
TRANSPARENT = (0, NODATA)


def build_luts():
    order = PRIORITY + [c for c in range(256) if c not in PRIORITY and c not in TRANSPARENT] + list(TRANSPARENT)
    code_to_rank = np.zeros(256, np.uint8)
    rank_to_code = np.zeros(256, np.uint8)
    for rank, code in enumerate(order):
        code_to_rank[code] = rank
        rank_to_code[rank] = code
    return code_to_rank, rank_to_code


CODE_TO_RANK, RANK_TO_CODE = build_luts()


def build_palette():
    pal = [255, 0, 255] * 256          # magenta flags any undocumented code
    alpha = [255] * 256
    for code, (_, colour) in CLASSES.items():
        if colour:
            pal[code * 3:code * 3 + 3] = [int(colour[i:i + 2], 16) for i in (1, 3, 5)]
    for code in TRANSPARENT:
        pal[code * 3:code * 3 + 3] = [0, 0, 0]
        alpha[code] = 0
    return pal, bytes(alpha)


PALETTE, ALPHA = build_palette()


def src_path(fx, fy):
    """fx: 0..3 west->east, fy: 0 = northern hemisphere, 1 = southern."""
    w = -180 + 90 * fx
    n, s = (90.0, 0.0) if fy == 0 else (0.0, -90.0)
    return os.path.join(SRC_DIR, f"gebco_2026_tid_n{n}_s{s}_w{w:.1f}_e{w + 90:.1f}_geotiff.tif")


def read_tiff(path):
    """Read an uncompressed, single-band uint8, one-row-per-strip GeoTIFF with plain numpy."""
    with open(path, "rb") as f:
        head = f.read(8)
        bo = "<" if head[:2] == b"II" else ">"
        ifd = struct.unpack(bo + "I", head[4:8])[0]
        f.seek(ifd)
        tags = {}
        for _ in range(struct.unpack(bo + "H", f.read(2))[0]):
            tag, typ, cnt, val = struct.unpack(bo + "HHI4s", f.read(12))
            tags[tag] = (typ, cnt, val)
        width = struct.unpack(bo + "H", tags[256][2][:2])[0]
        height = struct.unpack(bo + "H", tags[257][2][:2])[0]
        compression = struct.unpack(bo + "H", tags[259][2][:2])[0]
        assert (width, height, compression) == (NATIVE, NATIVE, 1), (path, width, height, compression)
        typ, cnt, val = tags[273]                     # StripOffsets
        f.seek(struct.unpack(bo + "I", val)[0])
        offsets = np.frombuffer(f.read(cnt * 4), dtype=bo + "u4")
        assert np.all(np.diff(offsets) == width), "strips are not contiguous"
    return np.fromfile(path, dtype=np.uint8, count=width * height, offset=int(offsets[0])).reshape(height, width)


def block_mode(arr, f):
    """Downsample by factor f: each output pixel is the most frequent code in its f x f block.
    Ties go to the higher-priority code. Computed from the full-resolution array (not cascaded),
    so class area fractions are preserved at every level."""
    h, w = arr.shape
    codes = sorted(np.unique(arr).tolist(), key=lambda c: CODE_TO_RANK[c])
    out = np.empty((h // f, w // f), np.uint8)
    band = max(f, 4096)                            # rows per chunk, bounds peak memory
    for r in range(0, h, band):
        sub = arr[r:r + band]
        best = np.zeros((sub.shape[0] // f, w // f), np.uint32)
        pick = np.zeros_like(best, np.uint8)
        for c in codes:                            # priority order, strict '>' keeps earlier on ties
            cnt = (sub == c).reshape(-1, f, w // f, f).sum(axis=(1, 3), dtype=np.uint32)
            better = cnt > best
            best[better] = cnt[better]
            pick[better] = c
        out[r // f:(r + band) // f] = pick
    return out


def reduce2x2(arr):
    return block_mode(arr, 2)


def save_tile(arr, level, x, y):
    """Write one tile; return False (and write nothing) if it is fully transparent."""
    if np.isin(arr, TRANSPARENT).all():
        return False
    d = os.path.join(OUT_DIR, str(level), str(x))
    os.makedirs(d, exist_ok=True)
    img = Image.fromarray(arr, "P")
    img.putpalette(PALETTE)
    img.save(os.path.join(d, f"{y}.png"), transparency=ALPHA, optimize=False)
    return True


def process_file(args):
    """Tile levels MAX_LEVEL..1 for one 90x90 deg source file."""
    fx, fy = args
    t0 = time.time()
    data = read_tiff(src_path(fx, fy))

    # Nearest-neighbour resample to the level-MAX pixel grid (pixel centres).
    n = 2 ** (MAX_LEVEL - 1)                       # tiles per file side at MAX_LEVEL
    size = n * TILE
    idx = ((np.arange(size) + 0.5) * NATIVE / size).astype(np.int64)
    full = data[idx][:, idx]
    del data

    available = []
    for level in range(MAX_LEVEL, 0, -1):
        arr = full if level == MAX_LEVEL else block_mode(full, 2 ** (MAX_LEVEL - level))
        n = 2 ** (level - 1)
        for j in range(n):
            for i in range(n):
                tile = arr[j * TILE:(j + 1) * TILE, i * TILE:(i + 1) * TILE]
                x, y = fx * n + i, fy * n + j
                if save_tile(tile, level, x, y):
                    available.append((level, x, y))
    print(f"  {os.path.basename(src_path(fx, fy))}: {len(available)} tiles in {time.time() - t0:.0f} s", flush=True)
    return (fx, fy), arr, available


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    jobs = [(fx, fy) for fy in range(2) for fx in range(4)]
    with Pool(len(jobs)) as pool:
        results = pool.map(process_file, jobs)

    level1 = {key: arr for key, arr, _ in results}
    available = [t for _, _, av in results for t in av]

    # Level 0: two 180x180 deg tiles, each stitched from 2x2 level-1 tiles.
    for x in range(2):
        stitched = np.block([[level1[(2 * x, 0)], level1[(2 * x + 1, 0)]],
                             [level1[(2 * x, 1)], level1[(2 * x + 1, 1)]]])
        if save_tile(reduce2x2(stitched), 0, x, 0):
            available.append((0, x, 0))

    # Availability bitmask per level: bit index = y * numX + x, numX = 2^(level+1).
    masks = {}
    for level in range(MAX_LEVEL + 1):
        nx, ny = 2 ** (level + 1), 2 ** level
        bits = np.zeros(nx * ny, np.uint8)
        for lv, x, y in available:
            if lv == level:
                bits[y * nx + x] = 1
        masks[str(level)] = base64.b64encode(np.packbits(bits, bitorder="little").tobytes()).decode()

    meta = {
        "name": "GEBCO_2026 TID grid",
        "attribution": "GEBCO Compilation Group (2026) GEBCO 2026 Grid (doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa)",
        "maxLevel": MAX_LEVEL,
        "tileSize": TILE,
        "classes": {str(c): {"label": l, "color": col} for c, (l, col) in CLASSES.items() if col},
        "groups": [{"name": g, "codes": [c for c in codes if c in CLASSES]} for g, codes in GROUPS],
        "available": masks,
    }
    with open(os.path.join(OUT_DIR, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False)
    print(f"Done: {len(available)} tiles, {time.time() - t0:.0f} s total")


if __name__ == "__main__":
    main()
