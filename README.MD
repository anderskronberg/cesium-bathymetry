# CesiumJS GEBCO bathymetry viewer

| Page | Shows |
|---|---|
| `index.html` | Barbados bathymetry from a Cesium ion terrain asset |
| `index_topo.html` | GEBCO_2026 depth/elevation as 3D terrain with an adjustable colour ramp and hillshade |
| `index_tid.html` | GEBCO_2026 Type Identifier (TID) grid – the source type behind each depth value |

## 1. Cesium ion access token

The pages need a Cesium ion access token (base maps, `index.html` terrain). The repository includes a `config.js` with a public demo token intended for the GitHub Pages deployment. It is not meant for production or other deployments, and it may be revoked or restricted at any time.

To use your own token instead:

1. Create a free account at https://ion.cesium.com and copy a token from https://ion.cesium.com/tokens.
2. Replace the token in `config.js` with your own (`config.example.js` is a clean template).
3. Do not commit your private token. To stop git from picking up your local change to `config.js`, run:

   ```
   git update-index --skip-worktree config.js
   ```

Consider restricting your token to your own URLs in the ion dashboard.

## 2. GEBCO data and tiles

The generated tiles are in the repository up to zoom level 7, for both `tiles/gebco_topo` and `tiles/gebco_tid` (~2.4 GB, ~80 000 files, so the clone is large). The pages work as soon as the repository is cloned.

The GEBCO source GeoTIFFs are not in the repository. You only need them to rebuild the tiles, for example after changing the build scripts:

1. Download the GEBCO_2026 Grid in GeoTIFF format from https://www.gebco.net/data-products/gridded-bathymetry-data:
   - sub-ice topography → extract to `gebco_2026_sub_ice_topo_geotiff/`
   - TID grid → extract to `gebco_2026_tid_geotiff/`
2. Build the tiles (Python 3 with numpy and Pillow):

   ```
   python -m venv .venv
   .venv\Scripts\python -m pip install numpy pillow
   .venv\Scripts\python build_topo_tiles.py    # -> tiles/gebco_topo  levels 0–7 (~2 min)
   .venv\Scripts\python build_tid_tiles.py     # -> tiles/gebco_tid   levels 0–8 (~14 min)
   ```

`build_tid_tiles.py` also builds level 8. That level is git-ignored, and the committed `tiles/gebco_tid/meta.json` sets `maxLevel` to 7. A rebuild writes `maxLevel: 8` to that file. Use the rebuilt file for local work with level 8, and do not commit it.

## 3. Run

The pages must be served over HTTP (they fetch tiles and `meta.json`), not opened as files:

```
python -m http.server 8888
```

Then open http://localhost:8888/index_topo.html (or `index.html`, `index_tid.html`).

## Attribution

GEBCO Compilation Group (2026) GEBCO 2026 Grid (doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa). Not to be used for navigation.
