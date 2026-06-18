"""Extract and crop Fig 2 right panel properly."""
from PIL import Image
import numpy as np

img = Image.open('D:/work/agent work/mnt-sim/page_2.png').convert('RGB')
arr = np.array(img, dtype=float)
h, w = arr.shape

# The vision analysis confirmed the figure is in the upper portion
# Grid score showed (1,2)=0.209 and (2,2)=0.209
# Typical figure placement: top center
# Let me crop a generous region from the upper portion

# Page 2 at 400 DPI: 3308x4410
# Figure likely occupies the upper 40% and spans most of the width
# But this is a single panel of a TWO-panel figure
# The PDF might have Fig. 1 and Fig. 2 on different pages

# Let me check page 3 as well - Fig. 2 might be on page 3
img3 = Image.open('D:/work/agent work/mnt-sim/page_3.png').convert('RGB')
arr3 = np.array(img3, dtype=float)
h3, w3 = arr3.shape

# Quick colorful check on page 3
rg3 = np.abs(arr3[:,:,0] - arr3[:,:,1])
rb3 = np.abs(arr3[:,:,0] - arr3[:,:,2])
gb3 = np.abs(arr3[:,:,1] - arr3[:,:,2])
colorful3 = np.maximum(np.maximum(rg3, rb3), gb3) > 40

for ri in range(8):
    for ci in range(6):
        y1, y2 = ri*h3//8, (ri+1)*h3//8
        x1, x2 = ci*w3//6, (ci+1)*w3//6
        score = colorful3[y1:y2, x1:x2].mean()
        if score > 0.01:
            print(f'Page 3 ({ri},{ci}): {score:.4f}')

# Also let me just crop the upper third of page 2 carefully
# The figure should be recognized by having black borders + white interior + colored points
# Try: upper 40% of page
crop_top = img.crop((0, 0, w, int(h*0.4)))
crop_top.save('D:/work/agent work/mnt-sim/page2_top.png')
print(f'Saved page2_top.png ({w}x{int(h*0.4)})')

# And try a wider scan of the middle of page 3
crop_mid3 = img3.crop((0, int(h3*0.15), w3, int(h3*0.65)))
crop_mid3.save('D:/work/agent work/mnt-sim/page3_mid.png')
print(f'Saved page3_mid.png')
