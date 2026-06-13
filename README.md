# AmenStyleSPECT — SPECT Brain Render (Amen-style)

A 3D Slicer extension that turns a reconstructed perfusion **SPECT** volume into two
Amen-style 3D renderings:

- **Surface scan** — a shaded cortical surface at the perfusion threshold;
  hypoperfusion shows as dents / holes (a hole is *below threshold*, not zero
  perfusion). Defaults to a marching-cubes isosurface with the Amen positional
  spectrum (pink superior → violet inferior); a GPU volume render is also available.
- **Active scan** — a blue wireframe brain plus the most-active tissue
  (hyperperfusion) rendered opaque red (≥ 85 %) → white (≥ 92 % of cerebellar max).
- **Regional quantification** — a per-region table reading each of the Amen
  regions as a % of the cerebellar reference, classified hypoactive / normal /
  hyperactive, with left/right asymmetry flags and a TSV export.

All thresholds are a **percent of the Maximum Cerebellar Count** (relative,
within-patient scaling). Bands: hyperactive > 80 %, normal 60–80 %, hypoactive < 60 %;
surface dents at 55 %.

> ⚠️ **Research / visualization aid only.** NOT FDA/CE-cleared and NOT for primary
> diagnosis. It renders *relative* within-patient perfusion with fixed percent
> thresholds; it does **not** perform normative-database (z-score) statistics.
> Validate thresholds against reference scans before any interpretive use. SPECT has
> no intrinsic left/right landmark — **confirm orientation on every scan before
> reading laterality.**

---

## Requirements

- **3D Slicer 5.x** (developed and tested on 5.10). Download: <https://download.slicer.org>.
- Nothing else. This is a **pure scripted (Python) module** — no compilation, no CMake
  build, no extra Python packages. It runs the same on **Windows, macOS, and Linux**.

---

## Install on another PC

You only need the `AmenStyleSPECT/SPECTBrainRender/` folder (the rest is for source
control / future Extensions-Index packaging). Copy the whole `AmenStyleSPECT` folder
to the target PC — e.g. unzip `AmenStyleSPECT-<version>.zip` somewhere stable like
`Documents\Slicer\AmenStyleSPECT` (Windows) or `~/Slicer/AmenStyleSPECT` (macOS/Linux).

### Method A — GUI (recommended, always works)

1. Open Slicer → **Edit ▸ Application Settings ▸ Modules**.
2. Under **Additional module paths**, click **Add** and select the folder that contains
   `SPECTBrainRender.py`:
   - Windows: `...\AmenStyleSPECT\SPECTBrainRender`
   - macOS/Linux: `.../AmenStyleSPECT/SPECTBrainRender`
   *(Select the `SPECTBrainRender` folder itself — not the outer `AmenStyleSPECT` folder.)*
3. (Optional) tick **Enable developer mode** to get a **Reload** button while testing.
4. Click **OK** and **restart** Slicer when prompted.
5. The module appears under **Nuclear Medicine ▸ SPECT Brain Render (Amen-style)**
   (use the module search / magnifier if you don't see it).

### Method B — one-paste install helper (convenience)

1. Copy the extension folder to the target PC.
2. Open Slicer → **View ▸ Python Console**.
3. Open `install_helper.py` (next to this README), set the `MODULE_DIR` line at the top
   to the `SPECTBrainRender` path on **that** PC, then paste the whole file into the
   console and press Enter. It registers the path, persists it, and offers to restart.
4. If for any reason the module doesn't appear after restart, use **Method A** — it is
   the authoritative path.

### Updating to a newer version

Replace the folder contents with the new files (keep the same location so the module
path you registered still points at it). Then either restart Slicer, or — in developer
mode — open the module and click **Reload**.

---

## Usage (quick walkthrough)

1. **Input** — select your already-loaded, reconstructed SPECT scalar volume.
2. **Orientation** — click **Show mid-axial slice**, set the L/R · A/P · S/I flips so the
   patient's left is where you expect, then tick **Orientation confirmed**.
3. **Reference (the 100 % anchor)** — keep mode on *Cerebellar max* and either:
   - type the Maximum Cerebellar Count into **Manual reference (counts)**, or
   - **Create cerebellum ROI** → drag it over the cerebellum → **Find cerebellar max**
     → step through candidates with **Jump to #** → **Accept # as 100 %**.
   The bold label shows the resolved reference. Toggle the candidate crosshair / ROI box
   off to declutter.
4. **Preprocess / brain mask** — largest-island masking (drops scalp / sinus / salivary
   uptake) is on by default.
5. **Surface scan** — pick the technique (isosurface mesh, the default, or volume render),
   set the threshold (default 55 %), opacity (default 1.0), and smoothness (default 0.3), pick a
   color scheme (Amen spectrum by default; or solid sculpt, hole-emphasis, perfusion rainbow),
   then **Generate Surface**. Use **Drop to 45 %** for the prognostic view.
6. **Active scan** — set the active threshold / **red** point (default 85 %), the **white** point
   (default 92 % of cerebellar max, the 100 % ceiling), and wireframe options, then **Generate
   Active**. Hyperperfused tissue renders red at the threshold and ramps to white at the white point.
7. **Views / export** — snap the camera to any of the 6 standard orientations, or pick an
   **Output folder** and click **Export Surface 6-view PNG** / **Export Active 6-view PNG**
   to write 6 PNGs + a `montage.png`.
8. **Regional quantification** — in **Regions to include**, tick the lobes (or expand a lobe
   and tick individual regions) you actually want to read — everything is on by default, so
   *Select none* then pick a few, or just untick what you don't need. Click
   **Create selected region ROIs** to drop one named Box ROI per ticked region (approximate
   starting positions) under an *Amen Regions* folder in the Data module, then drag each box
   onto its true structure (SPECT has limited anatomical detail, so this step is manual and
   clinician-judged). *Remove unchecked region ROIs* (off by default) additionally deletes any
   region boxes you've unticked, so the scene holds exactly your selection — it only ever
   touches boxes from this taxonomy, never your cerebellum reference or custom ROIs.
   Choose the region metric (mean / peak) and reference (% of cerebellar **max**, or % of the
   cerebellum-ROI **mean**), set the hypo / hyper / floor / asymmetry thresholds, and click
   **Quantify regions** (it reads whichever region ROIs are in the scene). The table sorts
   most-hypoactive first, tints rows by class, shows each region's **Hot %** (the fraction of
   the ROI above the hyper threshold), and flags left/right asymmetries; **double-click a row**
   to jump the slices to that region. **Export regional table (TSV)** writes `regional_all.tsv`
   and `regional_hypo_hyper.tsv` to the Views/export output folder.

   Quantification reads the **same brain-masked volume the scans render** (largest-island mask, when
   the *Preprocess / brain mask* option is on), so extracranial uptake — scalp, sinus, salivary,
   venous — that is stripped out of the render and never shows as red is **not counted** as region
   signal. (Turn the brain mask off to quantify the raw volume.) *Hypoactive* is judged from the
   region **mean** (low overall perfusion). *Hyperactive* is
   judged from the **hot voxels**: with **Focal hyperactivity** on (default), a region is flagged
   hyperactive when its hot focus reaches **Focal: min hot voxels** (default 3) **or** **Focal: or
   min hot fraction %** of the ROI (whichever is met first) — so a true focal hot spot is not
   diluted to "normal" by the surrounding normal tissue in a loose ROI. This reproduces the older
   "mask at 85 %, average inside the mask" read. The **Hot %** column shows each region's hot
   fraction *and its hot-voxel count in parentheses* — tune the thresholds against those numbers
   (e.g. if a focus reads `4.8% (24)` but isn't flagged, the 24 hot voxels already clear the
   default count of 3; if noise over-flags, raise the count). Turn Focal hyperactivity off to
   classify hyperactivity purely from the region mean/peak.

   Each abnormal region is then given an **ordinal grade** (with cerebellum = 100 % ceiling):
   **hyperperfusion +1 … +4** and **hypoperfusion −1 … −4**, shown in a **Grade** column and shaded
   deep-blue (−4) → green (0) → deep-red (+4).
   - **Hyper (+1…+4)** by the hot-voxel mean as a % of cerebellar max, using the active-scan
     anchors: **+1** red (85–92 %), **+2** white (92–100 %), **+3** above ceiling (100–110 %),
     **+4** (≥ 110 %). The *Grade +2 / white at* and *Grade +4 at* spin boxes set the boundaries.
   - **Hypo (−1…−4)** by the region mean, stepping below the hypoactive threshold (55 %) by the
     *Hypo grade step* (default 10 %): **−1** 45–55 %, **−2** 35–45 %, **−3** 25–35 %, **−4** < 25 %.

   The negative grades are the *surface-scan* (hypoperfusion) read, the positive grades the
   *active-scan* (hyperperfusion) read. The TSV carries a signed `Grade` column.

The **Status / calibration log** prints the resolved reference and the absolute count
value behind each percent threshold on every render, and the full regional table.

### Nodes this creates in the scene
`BrainMasked (Surface)` and `BrainMasked (Active)` (the two scans render from **separate**
masked clones, each with its own volume-rendering display/property node, so they never
overwrite each other), `BlueBrain` + `BlueNodes` (active wireframe), the `AmenSurface` model
(the marching-cubes isosurface surface scan), and the named region Box ROIs under the
*Amen Regions* folder. Your original volume is never modified.

---

## Packaging a distributable zip (maintainer)

From the project root (`slicerAmen/`):

```bash
./AmenStyleSPECT/package.sh           # writes AmenStyleSPECT-<version>.zip
```

The script zips the `AmenStyleSPECT` folder while excluding `.DS_Store`, `__pycache__`,
and any local junk. Copy the resulting zip to other PCs and install per the steps above.

## Going further (optional)

The repo keeps a root `CMakeLists.txt` and per-module `CMakeLists.txt` so the extension
can later be built/submitted to the **Slicer Extensions Index** for one-click install via
the Extensions Manager. That requires building against a matching Slicer and is **not**
needed for the manual install above.

## Files

```
AmenStyleSPECT/
  CMakeLists.txt                         extension metadata (for Extensions Index only)
  README.md                              this file
  package.sh                             build a distributable zip
  install_helper.py                      one-paste path registration for a new PC
  SPECTBrainRender/
    CMakeLists.txt                       scripted-module build rules
    SPECTBrainRender.py                  module (Widget + Logic + self-tests)
    Resources/
      Icons/SPECTBrainRender.png
      UI/SPECTBrainRender.ui
      Presets/AmenColor.vp.json          reference copy of the clinic color preset
```

**Version:** 0.6.0 (pre-release; surface + active renders with a marching-cubes isosurface
surface and the Amen positional spectrum coloring, cerebellar reference, 6-view montage export,
regional hypo/hyper quantification table). DICOM NM import and `CALIBRATION.md` are still to come.
