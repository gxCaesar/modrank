#!/usr/bin/env python3
"""Retrieve the two image assets Figure 3a uses, and record the geometry that links them.

WHY THIS EXISTS. The figure shows a real diagnostic slide rather than a drawing, so the figure is
only reproducible if the retrieval is. This script is the retrieval, and it writes a metadata file
recording where the magnified tile sits on the thumbnail so the builder does not carry a
hand-tuned offset that nobody can check.

WHAT IS FETCHED. Case TCGA-HQ-A2OF, diagnostic slide file 5e4ef3f3-5617-45bc-95ee-13fd5c6d0767.
The case is in this study's own 359-case partition, which is the reason it was chosen over any
other bladder slide. Two things are read through the GDC tile service:

  * a whole-slide overview at DeepZoom level 11, three tiles wide, stitched and cropped to the
    left tissue section;
  * one tile at level 14. The slide is 71,959 x 25,018 pixels and level 14 is a factor of eight
    down from full resolution, so a 512-pixel tile there covers 4,096 source pixels, which is the
    patch size the benchmark's own slide pipeline works at.

RIGHTS. GDC open-access data require no authorization and GDC states no restriction on analysis or
publication beyond the prohibition on attempting participant reidentification, so this use is
permitted. That is NOT the same as a transferable licence: no source consulted places an individual
TCGA slide in the public domain or under CC BY, so the figure claims no such licence. The
acknowledgement the source asks for is in the manuscript.

Run it only to regenerate the committed PNGs. The figure build does not call it and does not need
the network.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FILE_ID = "5e4ef3f3-5617-45bc-95ee-13fd5c6d0767"
CASE = "TCGA-HQ-A2OF"
API = "https://api.gdc.cancer.gov/tile/" + FILE_ID
OVERVIEW_LEVEL, OVERVIEW_COLS = 11, (0, 1, 2)
TILE_LEVEL, TILE_COL, TILE_ROW = 14, 4, 3
BG = 0.88                     # anything lighter than this is slide background, not tissue


def get(url, dest):
    r = subprocess.run(["curl", "-s", "--max-time", "60", "-o", dest, "-w", "%{http_code}", url],
                       capture_output=True, text=True)
    if r.stdout.strip() != "200":
        raise RuntimeError("GDC returned %s for %s" % (r.stdout.strip(), url))
    return dest


def main():
    tmp = tempfile.mkdtemp()
    meta = json.loads(open(get("https://api.gdc.cancer.gov/tile/metadata/" + FILE_ID,
                              os.path.join(tmp, "m.json"))).read())
    W, H, TS = int(meta["Width"]), int(meta["Height"]), int(meta["TileSize"])
    print("  slide %s x %s, tile size %s" % (W, H, TS))

    # --- overview: fetch the row of tiles and drop the one-pixel DeepZoom overlap at each seam
    parts = []
    for c in OVERVIEW_COLS:
        im = plt.imread(get("%s?level=%d&x=%d&y=0" % (API, OVERVIEW_LEVEL, c),
                            os.path.join(tmp, "o%d.png" % c)))[..., :3]
        parts.append(im[:, (1 if c else 0):(-1 if c != OVERVIEW_COLS[-1] else None)])
    stitched = np.concatenate(parts, axis=1)
    print("  overview stitched %s" % (stitched.shape,))

    # --- crop to the LEFT tissue section, and RECORD the offset so the builder need not guess
    mask = stitched.mean(axis=2) < BG
    half = mask[:, :stitched.shape[1] // 2]
    cols, rows = np.where(half.any(axis=0))[0], np.where(half.any(axis=1))[0]
    pad = 8
    x0, y0 = max(0, int(cols[0]) - pad), max(0, int(rows[0]) - pad)
    x1 = min(half.shape[1], int(cols[-1]) + pad)
    y1 = min(stitched.shape[0], int(rows[-1]) + pad)
    crop = stitched[y0:y1, x0:x1]
    plt.imsave(os.path.join(HERE, "blca_slide_thumb.png"), crop)
    print("  overview crop %s at offset (%d, %d)" % (crop.shape, x0, y0))

    tile = plt.imread(get("%s?level=%d&x=%d&y=%d" % (API, TILE_LEVEL, TILE_COL, TILE_ROW),
                          os.path.join(tmp, "t.png")))[..., :3][:TS, :TS]
    plt.imsave(os.path.join(HERE, "blca_slide_tile.png"), tile)

    # --- where the tile sits on the CROPPED overview, in cropped-overview pixels
    scale = 2 ** (TILE_LEVEL - OVERVIEW_LEVEL)          # overview pixels per tile-level pixel
    box = {"x": TILE_COL * TS / scale - x0, "y": TILE_ROW * TS / scale - y0,
           "w": TS / scale, "h": TS / scale}
    assert 0 <= box["x"] < crop.shape[1] and 0 <= box["y"] < crop.shape[0], (
        "the magnified tile does not fall inside the cropped overview: %s" % box)
    tissue = float((crop[int(box["y"]):int(box["y"] + box["h"]),
                         int(box["x"]):int(box["x"] + box["w"])].mean(axis=2) < BG).mean())
    assert tissue > 0.5, ("the locator box must land on tissue, not background; it covers %.2f"
                          % tissue)

    json.dump({"case": CASE, "file_id": FILE_ID, "slide_width_px": W, "slide_height_px": H,
               "overview_level": OVERVIEW_LEVEL, "tile_level": TILE_LEVEL,
               "tile_col": TILE_COL, "tile_row": TILE_ROW,
               "tile_covers_source_px": TS * 2 ** (int(np.ceil(np.log2(max(W, H)))) - TILE_LEVEL),
               "crop_offset_px": [x0, y0], "locator_box_px": box,
               "locator_box_tissue_fraction": round(tissue, 3),
               "rights": "GDC open access. Analysis and publication permitted; no transferable "
                         "licence is claimed for the slide itself.",
               "acknowledgement": "The results published here are in whole or part based upon data "
                                  "generated by the TCGA Research Network: "
                                  "https://www.cancer.gov/tcga."},
              open(os.path.join(HERE, "blca_slide_asset.json"), "w"), indent=1)
    print("  locator box %s, tissue fraction %.2f" % (box, tissue))
    return 0


if __name__ == "__main__":
    sys.exit(main())
