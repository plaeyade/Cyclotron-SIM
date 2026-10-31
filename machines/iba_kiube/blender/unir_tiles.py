"""Une las regiones de render_kiube.py --tiles N en una imagen, con fundido lineal en los solapes."""
import json, sys
import numpy as np
from PIL import Image

meta = json.load(open(sys.argv[1])); out = sys.argv[2]
W, H = meta["W"], meta["H"]
acc = np.zeros((H, W, 3), np.float64); wsum = np.zeros((H, W, 1), np.float64)
for t in meta["tiles"]:
    im = np.asarray(Image.open(t["f"]).convert("RGB"), np.float64)
    x0, x1 = t["x0"], t["x1"]
    top, bot = H - t["y1"], H - t["y0"]                 # Blender cuenta y desde abajo
    h, w = bot - top, x1 - x0
    im = im[:h, :w]
    ov = meta["solape"]
    wx = np.ones(w); wy = np.ones(h)
    if x0 > 0: wx[:2 * ov] = np.linspace(0, 1, 2 * ov)
    if x1 < W: wx[-2 * ov:] = np.minimum(wx[-2 * ov:], np.linspace(1, 0, 2 * ov))
    if top > 0: wy[:2 * ov] = np.linspace(0, 1, 2 * ov)
    if bot < H: wy[-2 * ov:] = np.minimum(wy[-2 * ov:], np.linspace(1, 0, 2 * ov))
    wgt = (wy[:, None] * wx[None, :])[..., None] + 1e-9
    acc[top:bot, x0:x1] += im * wgt; wsum[top:bot, x0:x1] += wgt
Image.fromarray(np.clip(acc / wsum, 0, 255).astype(np.uint8)).save(out, optimize=True)
print("unido →", out, f"{W}x{H}")
