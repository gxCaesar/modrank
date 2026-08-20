#!/usr/bin/env python3
"""Retrieve one overview thumbnail per case in the ten-case series, for the slide-array panel.

WHY. The manuscript argues about morphology and, until this panel, never showed any. These ten
cases are not a decorative selection: they are the five the frozen method scores highest and the
five it scores lowest, chosen by a rule fixed before the run, and every one of them is in the
359-case partition. Putting their actual histology side by side is the one place a reader can look
at the tissue and judge for themselves what the image arm is responding to.

HOW THE LEVEL IS CHOSEN. Slides differ in pixel size by more than a factor of three, so a fixed
DeepZoom level would return a 200-pixel thumbnail for one case and an 800-pixel one for the next.
The level is therefore chosen per slide as the smallest that still yields at least TARGET_PX across,
which keeps every cell at the same physical resolution in the figure and keeps the download small.

RIGHTS. As for the single-slide asset: GDC open-access data carry no restriction on analysis or
publication beyond the prohibition on reidentification, which permits this use. No source places an
individual slide under a transferable licence, and none is claimed. The acknowledgement TCGA asks
for is in the manuscript.

Run only to regenerate the committed PNGs. The figure build does not call it and needs no network.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile
import urllib.parse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = os.path.join(HERE, os.environ.get("SLIDE_OUT", "case_series"))
TARGET_PX = int(os.environ.get("SLIDE_PX", "340"))               # each grid cell renders about 0.9 in wide, so this is ~380 dpi
BG = 0.88


def curl(url, dest):
    r = subprocess.run(["curl", "-s", "--max-time", "60", "-o", dest, "-w", "%{http_code}", url],
                       capture_output=True, text=True)
    if r.stdout.strip() != "200":
        raise RuntimeError("GDC returned %s for %s" % (r.stdout.strip(), url))
    return dest


def uuid_for(file_name, tmp):
    # the filter is JSON and must be percent-encoded, or curl silently returns nothing at all
    filt = json.dumps({"op": "in", "content": {"field": "file_name", "value": [file_name]}})
    q = ("https://api.gdc.cancer.gov/files?filters=" + urllib.parse.quote(filt)
         + "&fields=file_id,file_size&format=JSON&size=1")
    hits = json.loads(open(curl(q, os.path.join(tmp, "q.json"))).read())["data"]["hits"]
    if not hits:
        raise RuntimeError("no GDC record for %s" % file_name)
    return hits[0]["file_id"]


def thumbnail(uuid, tmp):
    meta = json.loads(open(curl("https://api.gdc.cancer.gov/tile/metadata/" + uuid,
                                os.path.join(tmp, "m.json"))).read())
    W, H, TS = int(meta["Width"]), int(meta["Height"]), int(meta["TileSize"])
    top = int(math.ceil(math.log2(max(W, H))))
    level = next(l for l in range(top + 1) if W / 2.0 ** (top - l) >= TARGET_PX)
    lw = int(math.ceil(W / 2.0 ** (top - level)))
    lh = int(math.ceil(H / 2.0 ** (top - level)))
    cols, rows = int(math.ceil(lw / TS)), int(math.ceil(lh / TS))
    grid = []
    for r in range(rows):
        row = []
        for c in range(cols):
            im = plt.imread(curl("https://api.gdc.cancer.gov/tile/%s?level=%d&x=%d&y=%d"
                                 % (uuid, level, c, r),
                                 os.path.join(tmp, "t.png")))[..., :3]
            row.append(im[(1 if r else 0):(-1 if r != rows - 1 else None),
                          (1 if c else 0):(-1 if c != cols - 1 else None)])
        h = min(x.shape[0] for x in row)
        grid.append(np.concatenate([x[:h] for x in row], axis=1))
    w = min(x.shape[1] for x in grid)
    return np.concatenate([x[:, :w] for x in grid], axis=0), level, (W, H)


def crop_to_tissue(img, pad=6):
    mask = img.mean(axis=2) < BG
    if not mask.any():
        return img
    cols, rows = np.where(mask.any(axis=0))[0], np.where(mask.any(axis=1))[0]
    return img[max(0, rows[0] - pad):rows[-1] + pad, max(0, cols[0] - pad):cols[-1] + pad]


def main():
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    series = json.load(open(os.path.join(ROOT, "experiments", "20260817-blca-confirm", "results",
                                         "interpret-calibrate-cases.json")))
    cases = series["I_case_studies"]["cases"]
    if os.environ.get("SLIDE_NAMES"):        # a different case list drives the strip
        cases = [{"case_id": c} for c in json.load(open(os.path.join(
            HERE, os.environ["SLIDE_NAMES"])))]
    slides = json.load(open(os.path.join(HERE, os.environ.get("SLIDE_NAMES",
                                                          "case_slide_names.json"))))
    tmp = tempfile.mkdtemp()
    index = []
    for c in cases:
        cid = c["case_id"]
        name = slides[cid][0]
        uuid = uuid_for(name, tmp)
        img, level, dims = thumbnail(uuid, tmp)
        img = crop_to_tissue(img)
        path = os.path.join(OUT, cid + ".png")
        plt.imsave(path, img)
        index.append({"case_id": cid, "file_id": uuid, "level": level,
                      "slide_px": list(dims), "thumb_px": [img.shape[1], img.shape[0]],
                      "bytes": os.path.getsize(path)})
        print("  %-14s level %2d  %5d x %-5d  %6.0f KB"
              % (cid, level, img.shape[1], img.shape[0], os.path.getsize(path) / 1024))
    assert len(index) == len(cases), "fetched %d of %d" % (len(index), len(cases))
    json.dump({"series": "ten-case Series A, five highest and five lowest by the frozen method",
               "rights": "GDC open access. Analysis and publication permitted; no transferable "
                         "licence is claimed for any slide.",
               "cases": index}, open(os.path.join(OUT, "index.json"), "w"), indent=1)
    print("  total %.1f MB" % (sum(r["bytes"] for r in index) / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
