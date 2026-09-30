# CesiumJS GEBCO bathymetry viewer

| Page | Shows |
|---|---|
| `index.html` | Barbados bathymetry from a Cesium ion terrain asset |
| `index_topo.html` | GEBCO_2026 depth/elevation as 3D terrain with an adjustable colour ramp and hillshade |
| `index_tid.html` | GEBCO_2026 Type Identifier (TID) grid – the source type behind each depth value |

## 1. Cesium ion access token

The pages need a Cesium ion access token (base maps, `index.html` terrain). The token is **not** stored in the repository – use your own:

1. Create a free account at https://ion.cesium.com and copy a token from https://ion.cesium.com/tokens.
2. Copy `config.example.js` to `config.js` in the project root.
3. Replace `your_access_token` in `config.js` with your token.

`config.js` is listed in `.gitignore`, so it is never committed. Consider restricting the token to your own URLs in the ion dashboard.

## 2. GEBCO data and tiles

The GEBCO GeoTIFFs and the generated tiles are not in the repository (several GB).

1. Download the GEBCO_2026 Grid in GeoTIFF format from https://www.gebco.net/data-products/gridded-bathymetry-data:
   - sub-ice topography → extract to `gebco_2026_sub_ice_topo_geotiff/`
   - TID grid → extract to `gebco_2026_tid_geotiff/`
2. Build the tiles (Python 3 with numpy and Pillow):

   ```
   python -m venv .venv
   .venv\Scripts\python -m pip install numpy pillow
   .venv\Scripts\python build_topo_tiles.py    # -> tiles/gebco_topo  (~2 min, 2.2 GB)
   .venv\Scripts\python build_tid_tiles.py     # -> tiles/gebco_tid   (~14 min, 0.4 GB)
   ```

## 3. Run

The pages must be served over HTTP (they fetch tiles and `meta.json`), not opened as files:

```
python -m http.server 8888
```

Then open http://localhost:8888/index_topo.html (or `index.html`, `index_tid.html`).

## Attribution

GEBCO Compilation Group (2026) GEBCO 2026 Grid (doi:10.5285/4f68d5c7-45eb-f999-e063-7086abc036fa). Not to be used for navigation.
