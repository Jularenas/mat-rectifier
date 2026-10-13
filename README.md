# Mat Rectifier

Turns a phone photo of a pattern piece on the green cutting mat into a flat image at an exact scale, so a Procreate grid set in centimetres lines up with the mat's printed 1 cm grid.

## Use
1. Keep `index.html` and `grid.js` together in one folder. Open `index.html` in Chrome or Safari; it needs internet the first time to load OpenCV.js.
   - Or serve the folder: `python3 -m http.server 8765` and open http://localhost:8765
2. Drop one or more photos onto the page.
3. **Method: Auto, printed grid** (default). The page finds hundreds of 1 cm grid crossings (cyan dots) and fits the camera perspective, lens bending and any small local unevenness to them.
   - The mat does **not** need to be fully in the photo, and close-up or wide-angle shots are fine.
   - The info line shows the fit error in mm. Under 0.5 mm is excellent; over about 0.8 mm triggers a warning.
4. **Crop**:
   - *Visible mat, on grid lines* (default): the image edges sit exactly on printed grid lines, so a Procreate grid starting at the canvas corner matches with no offset.
   - *Whole mat, outer edge*: crops to the mat's physical edge (only when all 4 corners are in the photo) and reports the measured mat size as a sanity check.
5. Check the result with the **1 cm grid** overlay. The pink lines should sit on the printed lines everywhere. Use **Rotate 90°** if needed.
6. Click **Download** to get `<name>_rect_50pxcm.png`.

If the grid can't be found (very blurry or very small mat in the photo), the page falls back to **Manual: 4 mat corners**. In that mode you drag the 4 points onto the mat corners (with a magnifier), and the photo should be shot from straight above.

## Scale
- Default **50 px per cm**, written into the PNG as **127 DPI** (50 × 2.54), which is what Procreate uses to measure in cm. Change the px/cm under Settings; the DPI follows.

## In Procreate
1. Import the PNG, then check that **Actions → Canvas → Canvas information** shows 127 DPI.
   If it doesn't, create a custom canvas of the same pixel size at 127 DPI and insert the image into it.
2. **Actions → Canvas → Drawing Guide → Edit**: choose 2D Grid, set the units to cm and the grid size to 1 cm.

## Tips for photos
- Use the original camera photo when you can. WhatsApp or screenshot copies are downscaled, and below about 10 photo-px per cm the accuracy drops. The page warns when that happens.
- Leave some green grid visible around the piece on all sides. The grid under the paper is interpolated from the grid around it.
- Keep the mat flat; a wavy mat shows up as a higher fit error.
