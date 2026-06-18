"""Fast Fig. 2 detection using NumPy."""
import numpy as np
from PIL import Image

img = Image.open('D:/work/agent work/mnt-sim/page_2.png').convert('RGB')
arr = np.array(img, dtype=float)
h, w = arr.shape[:2]

# Colorfulness: max(|r-g|,|r-b|,|g-b|)
rg = np.abs(arr[:,:,0] - arr[:,:,1])
rb = np.abs(arr[:,:,0] - arr[:,:,2])
gb = np.abs(arr[:,:,1] - arr[:,:,2])
colorful = np.maximum(np.maximum(rg, rb), gb) > 40

# Divide into grid
nrows, ncols = 10, 8
max_score, best = 0, (0, 0)
scores = np.zeros((nrows, ncols))
for ri in range(nrows):
    for ci in range(ncols):
        y1, y2 = ri*h//nrows, (ri+1)*h//nrows
        x1, x2 = ci*w//ncols, (ci+1)*w//ncols
        score = colorful[y1:y2, x1:x2].mean()
        scores[ri, ci] = score
        if score > max_score:
            max_score = score
            best = (ri, ci)

print("Colorfulness score grid:")
print(np.array2string(scores, precision=4, suppress_small=True))
print(f"\nBest grid cell: ({best[0]}, {best[1]}), score={max_score:.4f}")

# Crop expanded region around best
ri, ci = best
y1 = max(0, (ri-1)*h//nrows)
y2 = min(h, (ri+2)*h//nrows)
x1 = max(0, (ci-1)*w//ncols)
x2 = min(w, (ci+2)*w//ncols)
print(f"Crop: ({x1},{y1})-({x2},{y2})")
img.crop((x1, y1, x2, y2)).save('D:/work/agent work/mnt-sim/fig2_detected.png')
print("Saved fig2_detected.png")
