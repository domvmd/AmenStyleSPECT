# AmenStyleSPECT — SPECT Brain Render (Amen-style)

A 3D Slicer extension that turns a reconstructed perfusion **SPECT** volume into two
Amen-style 3D renderings:

- **Surface scan** — a shaded cortical surface at the perfusion threshold;
  hypoperfusion shows as dents / holes (a hole is *below threshold*, not zero
  perfusion). Defaults to a marching-cubes isosurface with the Amen positional
  spectrum (pink superior → violet inferior); a GPU volume render is also available.
- **Active scan** — a blue wireframe brain plus the most-active tissue
  (hyperperfusion) rendered opaque red (≥ 85 %) → white (≥ 92 % of cerebellar max).
- **Regional quantification** — a per-region table reading each of the **58** Amen
  regions (29 paired L/R, incl. the precuneus as a default-mode-network hub and the
  nucleus accumbens) as a % of the cerebellar reference, classified hypoactive / normal /
  hyperactive, with an ordinal ±4 grade, left/right asymmetry flags and a TSV export. The
  `% Ref` cell is tinted with the Amen 20-step "hot" color scale so the table reads like
  the active render.
- **Automatic atlas quantification** *(v0.7)* — register the SPECT to a bundled MNI
  **HMPAO** template, warp the **AAL3** atlas into patient space, and read all 58 regions
  from their exact anatomical masks in **one click** — with a registration-quality verdict
  (PASS / REVIEW / FAIL) and a manual landmark rescue for atypical cases. The warped atlas
  is **named** (hover any region in the slice views to see it) and toggles on/off over the
  SPECT. No MRI required.

All thresholds are a **percent of the Maximum Cerebellar Count** (relative,
within-patient scaling), following the Amen reading convention. Bands: hot/hyperactive
≥ 85 %, normal 55–85 %, hypoactive < 55 %; surface dents at 55 %.

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
   most-hypoactive first, tints each region's **% Ref** cell with the Amen 20-step hot scale
   (blue → green → yellow → red → white, so it matches the active render) and shades the
   **Class / Grade** cells by severity, shows each region's **Hot %** (the fraction of
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

   Each abnormal region is then given an **ordinal grade** — the **Amen 20-step color scale**
   (each step = **5 %** of the cerebellar reference) read as signed deviation from the normal
   window: **hyperperfusion +1 … +4** and **hypoperfusion −1 … −4**, shown in a **Grade** column
   and shaded deep-blue (−4) → green (0) → deep-red (+4).
   - **Hyper (+1…+4)** by the hot-voxel mean, stepping up from the 85 % hot line:
     **+1** 85–90 %, **+2** 90–95 %, **+3** 95–100 %, **+4** ≥ 100 % (at/above the cerebellar
     ceiling).
   - **Hypo (−1…−4)** by the region mean, stepping down from the 55 % hypo line:
     **−1** 50–55 %, **−2** 45–50 %, **−3** 40–45 %, **−4** < 40 %.

   The single *Grade step (% per color)* spin box (default 5) sets the width of every step
   (so it stays a true 20-step scale). The active-scan render's white point (92 %) is set
   separately in the Active section and is independent of the grade ladder. The negative grades
   are the *surface-scan* (hypoperfusion) read, the positive grades the *active-scan*
   (hyperperfusion) read. The TSV carries a signed `Grade` column.

9. **Automatic atlas quantification** *(optional — replaces the manual box-dragging step)* —
   in **Regional quantification**, click **Register atlas to this SPECT (auto)**: it registers
   the scan to the bundled MNI HMPAO template, warps the **AAL3** atlas into the patient's space
   (~30–90 s), auto-fills the cerebellar reference, and reports a **PASS / REVIEW / FAIL**
   registration-quality verdict (containment, cerebellum/brain ratio, template-match NCC). The
   warped atlas is laid over the SPECT as a **named** colored overlay — **hover any region in the
   slice views** to read it in the Data Probe — and the **Show SPECT** / **Show atlas overlay**
   checkboxes + **Atlas opacity** slider let you blend or isolate either layer. Then
   click **Quantify from atlas (exact masks)** to read every region from its precise atlas voxels —
   the same table / report / TSV as the manual path, with no box-dragging. *Verify the auto-set
   reference, and confirm orientation first.*
   - On **REVIEW / FAIL**: rescue it — **Place patient landmarks** (6 named points appear at
     approximate spots, each carrying a placement hint), drag each onto its true location. Five are
     silhouette **edges** (frontal/occipital poles, vertex, L/R temporal); the **cerebellum** point
     is the **centre** of the bright cerebellar blob and has the most leverage. Then
     **Re-align atlas from landmarks**. If it still fails, fall back to the manual
     **Create selected region ROIs** workflow — the verdict is telling you the registration can't
     be trusted for this scan.
   - **Seed boxes from atlas** positions the manual region boxes at their atlas centroids (better
     starting points than the fixed scaffold) if you prefer the adjustable-box workflow.
   The atlas and template are **bundled** with the extension — nothing is downloaded.

The **Status / calibration log** prints the resolved reference and the absolute count
value behind each percent threshold on every render, and the full regional table.

### Nodes this creates in the scene
`BrainMasked (Surface)` and `BrainMasked (Active)` (the two scans render from **separate**
masked clones, each with its own volume-rendering display/property node, so they never
overwrite each other), `BlueBrain` + `BlueNodes` (active wireframe), the `AmenSurface` model
(the marching-cubes isosurface surface scan), the named region Box ROIs under the
*Amen Regions* folder, and — for the atlas pipeline — the warped `AAL3 (in patient)` labelmap
(with an `AAL3 region names` color table so it hover-labels each region) and the
`Patient landmarks` fiducials. Your original volume is never modified.

---

## License & attribution

This extension is free software under the **GNU General Public License v3** (see the bundled
`LICENSE`). It redistributes two third-party assets under their own terms:

- **AAL3 atlas** (`Resources/Atlas/AAL3v1.nii` + `AAL3v1.nii.txt`) — Automated Anatomical
  Labelling atlas 3, GIN-IMN (Bordeaux), released as copyright freeware under the **GNU GPL**.
  Cite: *Rolls, Huang, Lin, Feng & Joliot, "Automated anatomical labelling atlas 3," NeuroImage
  2020;206:116189*; original AAL: *Tzourio-Mazoyer et al., NeuroImage 2002;15:273–289*.
- **MNI HMPAO perfusion template** (`Resources/Atlas/SPECT_HMPAO_template.nii.gz`) — the Tc-99m
  HMPAO SPECT template distributed with **SPM** (Statistical Parametric Mapping; GPL), converted
  to NIfTI, used as the spatial-normalization target.

> The **AAL3 → Amen-taxonomy region grouping is a draft clinical mapping** (e.g. Broca's = pars
> opercularis only; the "Deep Limbic / Thalamic" region pools several small midbrain nuclei).
> Review it against your own conventions before any interpretive use.

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
  LICENSE                                GNU GPL v3
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
      AAL3v1.nii.txt                     AAL3 region lookup (label -> name)
      Atlas/AAL3v1.nii                   AAL3 atlas (MNI, GPL)
      Atlas/SPECT_HMPAO_template.nii.gz  MNI HMPAO perfusion template (from SPM)
```

**Version:** 0.9.0 (open-source, GPL v3). Aligns the grading and taxonomy to the published Amen
reading convention: the ordinal grade is now the **20-step / 5 %-per-color scale** (bands
hypo < 55 %, normal 55–85 %, hot ≥ 85 %), the regional table's `% Ref` cell is tinted with the
**Amen hot color LUT**, the **nucleus accumbens** is restored as its own region (**58 regions /
29 paired**), the warped atlas now **hover-labels** each region and **toggles on/off** over the
SPECT with an opacity slider, and the landmark-rescue points carry **placement hints**. Builds on
0.7.0 (automatic AAL3 atlas pipeline: in-module registration to a bundled MNI HMPAO template,
exact-mask quantification, PASS/REVIEW/FAIL QC, landmark rescue) and 0.6.0 (marching-cubes surface
scan + Amen positional spectrum, active scan, cerebellar reference, 6-view montage export).
DICOM NM import and `CALIBRATION.md` are still to come.
