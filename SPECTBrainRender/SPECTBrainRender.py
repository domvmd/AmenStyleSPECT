import logging
import os

import slicer
import vtk
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule,
    ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleTest,
    ScriptedLoadableModuleWidget,
)
from slicer.util import VTKObservationMixin

try:
    import qt
    import ctk
except ImportError:
    pass


# Reference (the count value that "100%" maps to). Cerebellar max is the Amen
# default; the other two are conveniences / fallbacks.
REF_CEREBELLUM = "Cerebellar max (Amen default)"
REF_ROBUST_MAX = "Robust whole-brain max"
REF_RAW_MAX = "Raw max"

# Amen rainbow color scale as (fraction-of-reference, (r, g, b)). Derived from
# the clinic preset "Brain SPECT.vp.json" whose cerebellar max was 870 counts,
# so the absolute control points are stored here as fractions of that reference
# and rescaled to each patient's reference at render time.
AMEN_COLOR_FRACTIONS = [
    (0.0843, (0.000000, 0.000000, 0.000000)),  # black
    (0.2302, (0.250980, 0.000000, 0.501961)),  # purple
    (0.3761, (0.000000, 0.000000, 1.000000)),  # blue
    (0.5221, (0.000000, 1.000000, 0.000000)),  # green
    (0.8139, (1.000000, 1.000000, 0.000000)),  # yellow
    (0.9599, (1.000000, 0.752941, 0.000000)),  # orange
    (1.2472, (1.000000, 0.011765, 0.000000)),  # red
]

# Surface rendering technique. The isosurface mesh (marching cubes at the
# perfusion threshold) is the research-recommended match for the current Amen
# 3D surface look: low-perfusion regions fall inward as crisp solid dents and
# the mesh shades like a sculpted brain. The GPU volume render is kept for
# comparison.
SURFACE_TECH_MESH = "Isosurface mesh (marching cubes)"
SURFACE_TECH_VR = "Volume render (GPU ray-cast)"

# Surface color schemes.
SURFACE_SCHEME_SPECTRUM = "Amen spectrum (positional rainbow)"
SURFACE_SCHEME_SOLID = "Solid sculpt (shape only)"
SURFACE_SCHEME_HOLE = "Hole-emphasis (purple->green->yellow)"
SURFACE_SCHEME_RAINBOW = "Perfusion rainbow (tint by uptake)"

# Single-material color for the "solid sculpt" surface: a warm neutral tissue tone.
SURFACE_SOLID_COLOR = (0.83, 0.79, 0.72)

# Amen "spectrum" surface palette. Sampling the reference Amen healthy surface
# image shows the hue is a smooth spectrum mapped to ANATOMICAL POSITION (an
# inferior->superior rainbow that rotates with the brain), while strong
# directional lighting carries the shape/shadow - NOT a perfusion map (hue stays
# constant within a region while only brightness varies). These RGB stops, keyed
# by normalized position along the chosen axis (0 = one end, 1 = the other),
# reproduce that look. With the axis/invert below: violet at the bottom
# (inferior) through to pink at the top (superior).
SURFACE_SPECTRUM = [
    (0.00, (0.59, 0.35, 0.89)),  # violet
    (0.18, (0.43, 0.55, 0.90)),  # blue
    (0.36, (0.35, 0.80, 0.85)),  # cyan
    (0.52, (0.50, 0.95, 0.55)),  # green
    (0.68, (0.95, 0.92, 0.40)),  # yellow
    (0.84, (0.98, 0.62, 0.38)),  # orange
    (1.00, (0.96, 0.45, 0.62)),  # pink
]
# RAS axis the spectrum runs along (2 = S/I) and whether to put 0 at the high end.
SURFACE_SPECTRUM_AXIS = 2
SURFACE_SPECTRUM_INVERT = False  # inferior (low z) -> spectrum 0 -> violet;
                                 # superior (high z) -> spectrum 1 -> pink

# Hole-emphasis: shadow / recessed tissue near holes (lower perfusion, just above
# the surface threshold) reads purple; brighter raised cortex (higher perfusion)
# reads green -> yellow. Stored as (fraction-of-reference, (r, g, b)).
HOLE_COLOR_FRACTIONS = [
    (0.45, (0.42, 0.05, 0.65)),  # deep purple  - shadow / near holes
    (0.58, (0.32, 0.18, 0.80)),  # indigo
    (0.70, (0.10, 0.65, 0.45)),  # green
    (0.82, (0.45, 0.85, 0.15)),  # yellow-green
    (0.95, (0.95, 0.92, 0.10)),  # yellow       - light tones
    (1.15, (1.00, 1.00, 0.55)),  # bright highlight
]

# Canonical Amen view orientations: name, camera direction (RAS: x=R, y=A, z=S),
# view-up. Inferior (underside) is critical for the temporal lobes.
VIEW_DIRECTIONS = [
    ("Left", (-1, 0, 0), (0, 0, 1)),
    ("Right", (1, 0, 0), (0, 0, 1)),
    ("Anterior", (0, 1, 0), (0, 0, 1)),
    ("Posterior", (0, -1, 0), (0, 0, 1)),
    ("Superior", (0, 0, 1), (0, 1, 0)),
    ("Inferior", (0, 0, -1), (0, 1, 0)),
]

# Node names (kept compatible with the clinic's standalone scripts).
MASKED_VOLUME_NAME = "BrainMasked"   # legacy name / default mask clone
# Surface and active scans are independent renderings, each with its OWN masked
# render-volume clone (and therefore its own volume-rendering display node), so
# generating one never overwrites the other's transfer functions and the
# surface smoothing never mutates the active volume.
SURFACE_RENDER_NAME = "BrainMasked (Surface)"
ACTIVE_RENDER_NAME = "BrainMasked (Active)"
WIRE_MODEL_NAME = "BlueBrain"
WIRE_NODES_NAME = "BlueNodes"
SURFACE_MODEL_NAME = "AmenSurface"  # marching-cubes fallback surface

#
# Regional quantification (ports the clinic's regional_table.py)
# ------------------------------------------------------------------
# Each region's MEAN or PEAK counts are expressed as a % of the cerebellar
# reference (the Maximum Cerebellar Count, or the cerebellum-region mean) and
# classified hypoactive (<60%) / normal / hyperactive (>80%).
REGION_METRIC_MEAN = "Region mean"
REGION_METRIC_PEAK = "Region peak (hottest voxel)"
REGION_REF_MAX = "% of cerebellar MAX"
REGION_REF_MEAN = "% of cerebellar-region MEAN"

DEFAULT_HYPO_PCT = 55.0          # region MEAN < this % of reference = hypoactive
DEFAULT_HYPER_PCT = 85.0         # voxels >= this % of reference are "hot"
DEFAULT_REGION_FLOOR_PCT = 5.0   # ignore voxels below this % of cereb. max (CSF/background)
DEFAULT_ASYMMETRY_PCT = 15.0     # |L-R| above this % of their mean is flagged
# Focal hyperactivity: a region is hyperactive if it contains a hot focus
# (voxels >= hyper % of the cerebellar max), regardless of the diluting mean.
# A focus counts when it has at least DEFAULT_MIN_HOT_VOXELS hot voxels OR fills
# at least DEFAULT_MIN_HOT_FRAC_PCT of the ROI - whichever is met first. The
# COUNT is the intuitive, ROI-size-independent dial (a small intense focus in a
# loose box is a small FRACTION but still a real cluster of voxels). Mirrors the
# clinic's old "mask at 85%, average inside the mask" read.
DEFAULT_MIN_HOT_VOXELS = 3
DEFAULT_MIN_HOT_FRAC_PCT = 10.0

# Ordinal perfusion grade (Amen-style), all as % of cerebellar max (100% ceiling).
# HYPERperfusion (+1..+4), graded on the hot-voxel mean using the active-scan
# red/white anchors:  +1 red (hyper threshold .. white), +2 white (white ..
# ceiling), +3 above ceiling (100 .. top), +4 (>= top).
DEFAULT_HYPER_WHITE_PCT = 92.0   # +2 onset = active-scan white point
DEFAULT_HYPER_TOP_PCT = 110.0    # +4 onset (the +3 onset is the 100% ceiling)
# HYPOperfusion (-1..-4), graded on the region mean stepping below the hypo
# threshold by this much per level (-1 = hypo..-step, -2 = -step..-2step, ...).
DEFAULT_HYPO_STEP_PCT = 10.0

REGION_FOLDER_NAME = "Amen Regions"  # subject-hierarchy folder for the ROI scaffold

# Per-class display colors for the result table and the ROI boxes.
CLASS_COLORS = {
    "Hypoactive": (0.20, 0.45, 1.00),   # cold = blue
    "Normal": (0.55, 0.78, 0.45),       # green
    "Hyperactive": (1.00, 0.28, 0.22),  # hot = red
}
# Ordinal grade colors: deep blue (-4) -> green (0) -> deep red (+4).
GRADE_COLORS = {
    -4: (0.04, 0.10, 0.50), -3: (0.13, 0.28, 0.80), -2: (0.33, 0.52, 1.00),
    -1: (0.62, 0.78, 1.00), 0: (0.55, 0.78, 0.45), 1: (1.00, 0.78, 0.45),
    2: (1.00, 0.52, 0.28), 3: (1.00, 0.28, 0.18), 4: (0.72, 0.04, 0.04),
}


def gradeLabel(grade):
    """'+2', '-3', or '0' (normal) for an ordinal grade integer."""
    return "%+d" % grade if grade else "0"


def rowDisplayColor(row):
    """Color for a region's ROI box and table cell, by ordinal grade."""
    return GRADE_COLORS.get(row.get("grade", 0), (0.7, 0.7, 0.7))

# Canonical Amen-style region taxonomy. Each spec is
#   (group, base name, paired?, clinical note, fy, fz, lat)
# and expands to ONE midline region (paired=False) or TWO (L/R) regions.
# fy/fz/lat are NORMALIZED anatomical fractions of the brain bounding box used
# ONLY to drop an APPROXIMATE starting Box ROI that the clinician then drags
# onto the true structure (SPECT has limited anatomical detail):
#   fy : 0 = posterior, 1 = anterior
#   fz : 0 = inferior,  1 = superior
#   lat: lateral offset from midline; L sits at 0.5-lat, R at 0.5+lat. In the
#        display RAS frame the patient's LEFT is the low-x side (matches the
#        Left/Right camera presets) - so confirm orientation before reading L/R.
# Revised taxonomy (user, 2026-06-11): 27 base regions, EVERY one paired L/R = 54 ROIs.
REGION_SPECS = [
    # 1. Frontal lobe (7) - paired L/R
    ("Frontal Lobe", "Prefrontal Pole", True, "", 0.95, 0.45, 0.18),
    ("Frontal Lobe", "Anterior Lateral PFC", True, "", 0.88, 0.52, 0.30),
    ("Frontal Lobe", "Dorsolateral PFC", True, "", 0.80, 0.74, 0.30),
    ("Frontal Lobe", "Posterior Frontal", True, "", 0.62, 0.72, 0.30),
    ("Frontal Lobe", "Inferior Orbital PFC", True, "", 0.88, 0.20, 0.18),
    ("Frontal Lobe", "Medial PFC", True, "", 0.86, 0.55, 0.07),
    ("Frontal Lobe", "Inferior Frontal Gyrus (Broca's)", True, "Broca's area (dominant hemisphere)", 0.72, 0.42, 0.34),
    # 2. Cingulate (5) - paired L/R, near midline
    ("Cingulate Gyrus", "Anterior Cingulate - Dorsal", True, "", 0.72, 0.62, 0.07),
    ("Cingulate Gyrus", "Anterior Cingulate - Genu", True, "", 0.80, 0.48, 0.07),
    ("Cingulate Gyrus", "Anterior Cingulate - Ventral", True, "vACC - predicts SSRI response", 0.78, 0.36, 0.07),
    ("Cingulate Gyrus", "Middle Cingulate", True, "", 0.52, 0.68, 0.07),
    ("Cingulate Gyrus", "Posterior Cingulate", True, "first to drop in Alzheimer's", 0.32, 0.58, 0.07),
    # 3. Temporal (6) - paired L/R
    ("Temporal Lobe", "Lateral Temporal - Anterior", True, "", 0.70, 0.30, 0.38),
    ("Temporal Lobe", "Lateral Temporal - Middle", True, "", 0.52, 0.30, 0.40),
    ("Temporal Lobe", "Lateral Temporal - Posterior", True, "", 0.34, 0.32, 0.38),
    ("Temporal Lobe", "Medial Temporal - Amygdala", True, "amygdala", 0.66, 0.28, 0.18),
    ("Temporal Lobe", "Medial Temporal - Hippocampus", True, "hippocampus", 0.50, 0.30, 0.18),
    ("Temporal Lobe", "Medial Temporal - Posterior", True, "parahippocampal gyrus", 0.38, 0.30, 0.18),
    # 4. Parietal (2) - paired L/R
    ("Parietal Lobe", "Superior Parietal", True, "", 0.28, 0.82, 0.28),
    ("Parietal Lobe", "Inferior Parietal", True, "", 0.30, 0.62, 0.36),
    # 5. Occipital (2) - paired L/R
    ("Occipital Lobe", "Primary Visual Cortex", True, "calcarine", 0.08, 0.42, 0.12),
    ("Occipital Lobe", "Visual Association Cortex", True, "", 0.14, 0.52, 0.30),
    # 6. Subcortical (4) - paired L/R
    ("Subcortical", "Basal Ganglia (Caudate/Putamen)", True, "", 0.55, 0.50, 0.18),
    ("Subcortical", "Thalamus", True, "", 0.45, 0.52, 0.09),
    ("Subcortical", "Insular Cortex", True, "", 0.58, 0.45, 0.30),
    ("Subcortical", "Deep Limbic / Thalamic", True, "deep limbic system", 0.52, 0.42, 0.09),
    # 7. Cerebellum (1) - paired L/R, reference standard, anchors 100%
    ("Cerebellum", "Cerebellum", True, "reference standard (~100%); medial/lateral", 0.18, 0.18, 0.22),
]

# Box ROI size as a fraction of each brain-bounding-box extent, by group.
GROUP_BOX_FRAC = {
    "Subcortical": 0.075,
    "Cingulate Gyrus": 0.085,
    "Temporal Lobe": 0.10,
}
DEFAULT_BOX_FRAC = 0.11


def expandRegionSpecs(specs):
    """Expand REGION_SPECS into a flat list of region dicts (one per ROI)."""
    regions = []
    for group, base, paired, note, fy, fz, lat in specs:
        if paired:
            regions.append(dict(group=group, base=base, side="L", name="L " + base,
                                note=note, fx=0.5 - lat, fy=fy, fz=fz))
            regions.append(dict(group=group, base=base, side="R", name="R " + base,
                                note=note, fx=0.5 + lat, fy=fy, fz=fz))
        else:
            regions.append(dict(group=group, base=base, side="M", name=base,
                                note=note, fx=0.5, fy=fy, fz=fz))
    return regions


AMEN_REGIONS = expandRegionSpecs(REGION_SPECS)
AMEN_REGIONS_BY_NAME = {r["name"]: r for r in AMEN_REGIONS}


def _qcount(getter):
    """PythonQt exposes parameterless Qt count getters (e.g. topLevelItemCount,
    childCount) as int PROPERTIES rather than methods; tolerate either form."""
    return getter() if callable(getter) else int(getter)


#
# SPECTBrainRender
#
class SPECTBrainRender(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "SPECT Brain Render (Amen-style)"
        self.parent.categories = ["Nuclear Medicine"]
        self.parent.dependencies = []
        self.parent.contributors = ["Dominic Velasco, MD (Brain SPECT Clinic)"]
        self.parent.helpText = """
<p><b>SPECT Brain Render</b> reproduces an Amen-style perfusion read from a
reconstructed SPECT volume, using the Maximum Cerebellar Count as the 100%
reference:</p>
<ul>
<li><b>Surface scan</b> &mdash; a solid shaded brain surface at the threshold
(% of reference); hypoperfusion falls inward as dents / holes (a hole is below
threshold, NOT zero perfusion). Choose the technique: a marching-cubes
isosurface (recommended &mdash; closest to the Amen look) or the GPU volume
render. Color it as a single sculpted tone or tint it by perfusion.</li>
<li><b>Active scan</b> &mdash; blue wireframe brain (anatomy you can see through)
plus the most-active tissue (&ge;80% of reference) volume-rendered opaque
red&rarr;pink&rarr;white.</li>
<li><b>Regional quantification</b> &mdash; one named Box ROI per region; each
region's mean / peak uptake is read as a % of the cerebellar reference and
classified hypoactive / normal / hyperactive, with left/right asymmetry flags
and a TSV export.</li>
</ul>
<p>Bands relative to cerebellar max: hyperactive &gt;80%, normal 60&ndash;80%,
hypoactive &lt;60%.</p>
"""
        self.parent.acknowledgementText = """
Visualization convention inspired by the Amen Clinics display style and ported
from the clinic's Slicer Python workflow. Research / visualization aid only
&mdash; NOT FDA/CE-cleared, NOT for primary diagnosis. Renders RELATIVE
within-patient perfusion with fixed percent thresholds; does NOT perform
normative-database (z-score) statistics. SPECT has no intrinsic L/R landmark
&mdash; confirm orientation on every scan before reading laterality.
"""


#
# SPECTBrainRenderWidget
#
class SPECTBrainRenderWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self.ui = None
        self._updatingFromCode = False
        self._baselineIJKToRAS = None
        # The accepted 100% reference (cerebellar max). 0 = not set yet.
        self.cerebellarMax = 0.0
        # Candidates from the last cerebellar-max search: list of dicts.
        self._candidates = []

    # ---- lifecycle -------------------------------------------------------
    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/SPECTBrainRender.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)
        uiWidget.setMRMLScene(slicer.mrmlScene)
        # The .ui has no <connections> wiring mrmlSceneChanged -> setMRMLScene,
        # so set the scene on every qMRML selector explicitly.
        self.ui.inputVolumeSelector.setMRMLScene(slicer.mrmlScene)
        self.ui.cerebellumRoiSelector.setMRMLScene(slicer.mrmlScene)

        self.logic = SPECTBrainRenderLogic()

        self.ui.referenceModeComboBox.addItems([REF_CEREBELLUM, REF_ROBUST_MAX, REF_RAW_MAX])
        self.ui.referenceModeComboBox.setCurrentText(REF_CEREBELLUM)
        self.ui.surfaceTechniqueComboBox.addItems([SURFACE_TECH_MESH, SURFACE_TECH_VR])
        self.ui.surfaceTechniqueComboBox.setCurrentText(SURFACE_TECH_MESH)
        self.ui.surfaceColorSchemeComboBox.addItems(
            [SURFACE_SCHEME_SPECTRUM, SURFACE_SCHEME_SOLID,
             SURFACE_SCHEME_HOLE, SURFACE_SCHEME_RAINBOW])
        self.ui.surfaceColorSchemeComboBox.setCurrentText(SURFACE_SCHEME_SPECTRUM)

        # Input
        self.ui.inputVolumeSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onInputVolumeChanged)
        self.ui.loadDicomButton.connect("clicked(bool)", self.onLoadDicom)

        # Orientation
        for cb in (self.ui.flipLRCheckBox, self.ui.flipAPCheckBox, self.ui.flipSICheckBox):
            cb.connect("toggled(bool)", self.onOrientationFlipChanged)
        self.ui.showOrientationButton.connect("clicked(bool)", self.onShowOrientation)

        # Reference
        self.ui.referenceModeComboBox.connect("currentTextChanged(QString)", self.onReferenceModeChanged)
        self.ui.createRoiButton.connect("clicked(bool)", self.onCreateCerebellumRoi)
        self.ui.findCerebellarMaxButton.connect("clicked(bool)", self.onFindCerebellarMax)
        self.ui.jumpToCandidateButton.connect("clicked(bool)", self.onJumpToCandidate)
        self.ui.acceptCandidateButton.connect("clicked(bool)", self.onAcceptCandidate)
        self.ui.manualReferenceSpinBox.connect("valueChanged(double)", self.onManualReferenceChanged)
        self.ui.showReticleCheckBox.connect("toggled(bool)", self.onToggleReticle)
        self.ui.showRoiCheckBox.connect("toggled(bool)", self.onToggleRoi)

        # Surface
        self.ui.surfaceTechniqueComboBox.connect("currentTextChanged(QString)", self.onSurfaceParamChanged)
        self.ui.surfaceThresholdSlider.connect("valueChanged(double)", self.onSurfaceParamChanged)
        self.ui.surfaceOpacitySlider.connect("valueChanged(double)", self.onSurfaceParamChanged)
        self.ui.surfaceSmoothnessSlider.connect("valueChanged(double)", self.onSurfaceParamChanged)
        self.ui.surfaceColorSchemeComboBox.connect("currentTextChanged(QString)", self.onSurfaceParamChanged)
        self.ui.keepPresetOpacityCheckBox.connect("toggled(bool)", self.onSurfaceParamChanged)
        self.ui.prognostic45Button.connect("clicked(bool)", self.onPrognostic45)
        self.ui.generateSurfaceButton.connect("clicked(bool)", self.onGenerateSurface)

        # Active
        self.ui.generateActiveButton.connect("clicked(bool)", self.onGenerateActive)

        # Views / export
        import os
        self.ui.outputDirButton.directory = os.path.expanduser("~")
        self._viewButtons = {
            "Superior": self.ui.viewSuperiorButton,
            "Inferior": self.ui.viewInferiorButton,
            "Anterior": self.ui.viewAnteriorButton,
            "Posterior": self.ui.viewPosteriorButton,
            "Left": self.ui.viewLeftButton,
            "Right": self.ui.viewRightButton,
        }
        for name, btn in self._viewButtons.items():
            btn.connect("clicked(bool)", lambda checked, n=name: self.onSnapView(n))
        self.ui.exportSurfaceButton.connect("clicked(bool)", self.onExportSurface)
        self.ui.exportActiveButton.connect("clicked(bool)", self.onExportActive)

        # Regional quantification
        self.ui.regionMetricComboBox.addItems([REGION_METRIC_MEAN, REGION_METRIC_PEAK])
        self.ui.regionRefModeComboBox.addItems([REGION_REF_MAX, REGION_REF_MEAN])
        self.ui.regionRefModeComboBox.setCurrentText(REGION_REF_MAX)
        self._configureRegionTable()
        self._updatingTree = False
        self._buildRegionTree()
        self.ui.regionSelectTreeWidget.connect("itemChanged(QTreeWidgetItem*,int)", self.onRegionTreeItemChanged)
        self.ui.selectAllRegionsButton.connect("clicked(bool)", lambda c: self._setAllRegions(qt.Qt.Checked))
        self.ui.selectNoneRegionsButton.connect("clicked(bool)", lambda c: self._setAllRegions(qt.Qt.Unchecked))
        self.ui.createRegionScaffoldButton.connect("clicked(bool)", self.onCreateRegionScaffold)
        self.ui.quantifyRegionsButton.connect("clicked(bool)", self.onQuantifyRegions)
        self.ui.exportRegionTsvButton.connect("clicked(bool)", self.onExportRegionTsv)
        self.ui.regionTableWidget.connect("cellDoubleClicked(int,int)", self.onRegionRowDoubleClicked)
        self._regionRows = []

        self.layout.addStretch(1)
        self.updateInputDependentEnabling()
        self.onReferenceModeChanged(self.ui.referenceModeComboBox.currentText)
        logging.info("SPECTBrainRender widget setup complete (Stage 4 - regional quantification).")

    def cleanup(self):
        self.removeObservers()

    # ---- helpers ---------------------------------------------------------
    def currentVolumeNode(self):
        return self.ui.inputVolumeSelector.currentNode()

    def updateInputDependentEnabling(self):
        hasVolume = self.currentVolumeNode() is not None
        for w in (
            self.ui.flipLRCheckBox, self.ui.flipAPCheckBox, self.ui.flipSICheckBox,
            self.ui.orientationConfirmedCheckBox, self.ui.showOrientationButton,
            self.ui.createRoiButton, self.ui.findCerebellarMaxButton,
            self.ui.generateSurfaceButton, self.ui.prognostic45Button,
            self.ui.surfaceThresholdSlider, self.ui.generateActiveButton,
            self.ui.viewSuperiorButton, self.ui.viewInferiorButton,
            self.ui.viewAnteriorButton, self.ui.viewPosteriorButton,
            self.ui.viewLeftButton, self.ui.viewRightButton,
            self.ui.exportSurfaceButton, self.ui.exportActiveButton,
            self.ui.createRegionScaffoldButton, self.ui.quantifyRegionsButton,
            self.ui.exportRegionTsvButton,
        ):
            w.setEnabled(hasVolume)

    def log(self, text, replace=False):
        logging.info(text)
        box = self.ui.statusTextEdit
        box.setPlainText(text) if replace else box.appendPlainText(text)

    def preprocessOptions(self):
        return dict(
            keepLargestIsland=self.ui.brainMaskCheckBox.checked,
            islandPct=self.ui.islandThresholdSpinBox.value,
            islandErode=int(self.ui.islandErodeSpinBox.value),
            robustPercentile=self.ui.robustPercentileSpinBox.value,
        )

    def currentReference(self, volumeNode, stats=None):
        """Resolve the 100% reference value for the chosen mode."""
        mode = self.ui.referenceModeComboBox.currentText
        return self.logic.getReferenceValue(volumeNode, mode, self.cerebellarMax, stats)

    def refreshReferenceLabel(self, volumeNode=None):
        node = volumeNode or self.currentVolumeNode()
        if node is None:
            self.ui.referenceValueLabel.setText("Reference (100%): (no volume)")
            return
        try:
            ref = self.currentReference(node)
        except Exception as exc:
            self.ui.referenceValueLabel.setText(f"Reference (100%): unavailable ({exc})")
            return
        mode = self.ui.referenceModeComboBox.currentText
        self.ui.referenceValueLabel.setText(
            f"Reference (100%) = {ref:.1f} counts  [{mode}]"
        )

    # ---- input callbacks -------------------------------------------------
    def onInputVolumeChanged(self, node):
        self.updateInputDependentEnabling()
        self._updatingFromCode = True
        for cb in (self.ui.flipLRCheckBox, self.ui.flipAPCheckBox, self.ui.flipSICheckBox):
            cb.setChecked(False)
        self.ui.orientationConfirmedCheckBox.setChecked(False)
        self._updatingFromCode = False
        self._baselineIJKToRAS = None
        if node is None:
            return
        mat = vtk.vtkMatrix4x4()
        node.GetIJKToRASMatrix(mat)
        self._baselineIJKToRAS = mat
        try:
            stats, _ = self.logic.computeStatistics(node, **{
                k: self.preprocessOptions()[k] for k in ("robustPercentile",)})
        except Exception as exc:
            self.log(f"Could not compute statistics: {exc}", replace=True)
            return
        self.refreshReferenceLabel(node)
        self.log(self.logic.formatReport(node, stats, self.cerebellarMax,
                                         self.ui.referenceModeComboBox.currentText), replace=True)
        self.log(
            "ORIENTATION: SPECT has no intrinsic L/R landmark. Confirm which side "
            "is the patient's LEFT and tick 'Orientation confirmed' before reading "
            "laterality."
        )
        self.onShowOrientation()

    def onLoadDicom(self):
        slicer.util.infoDisplay(
            "DICOM import (with NM multiframe validation/repair) arrives in a "
            "later stage. For now, select an already-loaded scalar volume node.",
            windowTitle="SPECT Brain Render",
        )

    # ---- orientation -----------------------------------------------------
    def onOrientationFlipChanged(self, _checked):
        if self._updatingFromCode:
            return
        node = self.currentVolumeNode()
        if node is None or self._baselineIJKToRAS is None:
            return
        self._updatingFromCode = True
        self.ui.orientationConfirmedCheckBox.setChecked(False)
        self._updatingFromCode = False
        self.logic.applyOrientationFlips(
            node, self._baselineIJKToRAS,
            self.ui.flipLRCheckBox.checked, self.ui.flipAPCheckBox.checked,
            self.ui.flipSICheckBox.checked,
        )
        self.onShowOrientation()

    def onShowOrientation(self):
        node = self.currentVolumeNode()
        if node is None:
            return
        slicer.util.setSliceViewerLayers(background=node, fit=True)
        slicer.app.applicationLogic().FitSliceToAll()

    # ---- reference / cerebellar max --------------------------------------
    def onReferenceModeChanged(self, mode):
        isCereb = (mode == REF_CEREBELLUM)
        for w in (self.ui.cerebellumRoiSelector, self.ui.createRoiButton,
                  self.ui.findCerebellarMaxButton, self.ui.candidateRankSpinBox,
                  self.ui.jumpToCandidateButton, self.ui.acceptCandidateButton,
                  self.ui.manualReferenceSpinBox):
            w.setEnabled(isCereb)
        self.refreshReferenceLabel()

    def onCreateCerebellumRoi(self):
        node = self.currentVolumeNode()
        if node is None:
            return
        roi = self.logic.createCerebellumRoi(node)
        self.ui.cerebellumRoiSelector.setCurrentNode(roi)
        self.log(
            "Created a Box ROI placed over the posterior fossa. Drag it in the "
            "slice views to loosely cover the CEREBELLUM, then click 'Find "
            "cerebellar max'."
        )

    def onFindCerebellarMax(self):
        node = self.currentVolumeNode()
        roi = self.ui.cerebellumRoiSelector.currentNode()
        if node is None or roi is None:
            self.log("Select a volume and a cerebellum ROI first.")
            return
        topN = int(self.ui.candidateRankSpinBox.maximum)
        try:
            self._candidates = self.logic.cerebellarMaxCandidates(node, roi, topN)
        except Exception as exc:
            self.log(f"Cerebellar max search failed: {exc}")
            return
        if not self._candidates:
            self.log("No voxels inside the ROI - is the box over the brain?")
            return
        lines = ["=== Cerebellar max candidates (hottest voxels in ROI) ===",
                 "Rank   Value        RAS (R, A, S) mm"]
        for c in self._candidates:
            lines.append("#%-3d %9.1f    (%7.1f, %7.1f, %7.1f)" %
                         (c["rank"], c["value"], c["ras"][0], c["ras"][1], c["ras"][2]))
        lines.append("Step through with 'Jump to #', confirm a TRUE cerebellar "
                     "voxel (not sinus/scalp), then 'Accept # as 100%'.")
        self.log("\n".join(lines))
        self.ui.candidateRankSpinBox.value = 1
        self.onJumpToCandidate()

    def onJumpToCandidate(self):
        if not self._candidates:
            return
        rank = int(self.ui.candidateRankSpinBox.value)
        cand = next((c for c in self._candidates if c["rank"] == rank), None)
        if cand is None:
            return
        self.logic.jumpToRAS(cand["ras"], showCrosshair=self.ui.showReticleCheckBox.checked)
        self.log("Jumped to #%d  value=%.1f" % (rank, cand["value"]))

    def onAcceptCandidate(self):
        if not self._candidates:
            return
        rank = int(self.ui.candidateRankSpinBox.value)
        cand = next((c for c in self._candidates if c["rank"] == rank), None)
        if cand is None:
            return
        self.cerebellarMax = float(cand["value"])
        self._updatingFromCode = True
        self.ui.manualReferenceSpinBox.value = self.cerebellarMax
        self._updatingFromCode = False
        self.log("Accepted #%d as Maximum Cerebellar Count = %.1f counts (100%% "
                 "reference)." % (rank, self.cerebellarMax))
        self.refreshReferenceLabel()

    def onManualReferenceChanged(self, value):
        if self._updatingFromCode:
            return
        self.cerebellarMax = float(value)
        self.refreshReferenceLabel()

    def onToggleReticle(self, checked):
        # Show / hide the crosshair "reticle" used when inspecting candidates.
        self.logic.setCrosshairVisible(checked)
        self.log("Candidate reticle " + ("shown." if checked else "hidden."))

    def onToggleRoi(self, checked):
        roi = self.ui.cerebellumRoiSelector.currentNode()
        if roi and roi.GetDisplayNode():
            roi.GetDisplayNode().SetVisibility(checked)
            self.log("Cerebellum ROI box " + ("shown." if checked else "hidden."))

    # ---- surface ---------------------------------------------------------
    def _requireReference(self, node):
        mode = self.ui.referenceModeComboBox.currentText
        if mode == REF_CEREBELLUM and self.cerebellarMax <= 0:
            self.log(
                "No cerebellar max set yet. Either enter it in 'Manual reference', "
                "use 'Find cerebellar max', or switch the reference mode. Falling "
                "back to robust whole-brain max for now."
            )
        return self.currentReference(node)

    def _surfaceVisible(self):
        """True if either surface representation (VR display node or mesh model)
        is currently shown - so a param/technique change live-regenerates."""
        dn = self.logic.surfaceDisplayNode
        if dn and slicer.mrmlScene.IsNodePresent(dn) and dn.GetVisibility():
            return True
        mn = self.logic.surfaceModelNode
        if (mn and slicer.mrmlScene.IsNodePresent(mn) and mn.GetDisplayNode()
                and mn.GetDisplayNode().GetVisibility()):
            return True
        return False

    def onSurfaceParamChanged(self, _value=None):
        if self._updatingFromCode:
            return
        if self._surfaceVisible():
            self._generateSurface(live=True)

    def onPrognostic45(self):
        self._updatingFromCode = True
        self.ui.surfaceThresholdSlider.value = 45
        self._updatingFromCode = False
        self._generateSurface(live=False)

    def onGenerateSurface(self):
        self._generateSurface(live=False)

    def _generateSurface(self, live):
        node = self.currentVolumeNode()
        if node is None:
            return
        if not self.ui.orientationConfirmedCheckBox.checked:
            self.log("WARNING: orientation not confirmed - do not read laterality yet.")
        ref = self._requireReference(node)
        opts = self.preprocessOptions()
        threshold = self.ui.surfaceThresholdSlider.value
        technique = self.ui.surfaceTechniqueComboBox.currentText
        scheme = self.ui.surfaceColorSchemeComboBox.currentText
        opacity = self.ui.surfaceOpacitySlider.value
        sigma = self.ui.surfaceSmoothnessSlider.value
        brainOnly = "Brain-only (largest island)." if opts["keepLargestIsland"] else ""
        prefix = "Live update" if live else "Surface scan"
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            if technique == SURFACE_TECH_MESH:
                result = self.logic.generateSurfaceMesh(
                    node, threshold, ref, opacity=opacity, smoothingSigma=sigma,
                    colorScheme=scheme, keepLargestIsland=opts["keepLargestIsland"],
                    islandPct=opts["islandPct"], islandErode=opts["islandErode"])
            else:
                result = self.logic.generateSurfaceVR(
                    node, threshold, ref, opacity=opacity, smoothingSigma=sigma,
                    colorScheme=scheme, keepLargestIsland=opts["keepLargestIsland"],
                    islandPct=opts["islandPct"], islandErode=opts["islandErode"],
                    keepPresetOpacity=self.ui.keepPresetOpacityCheckBox.checked)
        finally:
            qt.QApplication.restoreOverrideCursor()
        if technique == SURFACE_TECH_MESH:
            self.log("%s: marching-cubes isosurface at %.0f%% of reference (%.1f counts), "
                     "%s, %d vertices, opacity %.2f. %s" % (
                         prefix, threshold, threshold / 100.0 * ref, scheme,
                         result["points"], opacity, brainOnly))
        else:
            self.log("%s: %s VR surface, opacity %.2f, ramp at %.0f%% of reference "
                     "(%.1f counts). %s" % (
                         prefix, scheme, opacity, threshold, threshold / 100.0 * ref, brainOnly))
        self._frame3D()

    # ---- active ----------------------------------------------------------
    def onGenerateActive(self):
        node = self.currentVolumeNode()
        if node is None:
            return
        if not self.ui.orientationConfirmedCheckBox.checked:
            self.log("WARNING: orientation not confirmed - do not read laterality yet.")
        ref = self._requireReference(node)
        opts = self.preprocessOptions()
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            result = self.logic.generateActive(
                node, ref,
                activePct=self.ui.activePctSpinBox.value,
                whitePct=self.ui.activeWhitePctSpinBox.value,
                meshPct=self.ui.meshPctSpinBox.value,
                nodeRadius=self.ui.nodeRadiusSpinBox.value,
                decimate=self.ui.decimateSpinBox.value,
                wireOpacity=self.ui.wireOpacitySpinBox.value,
                keepLargestIsland=opts["keepLargestIsland"],
                islandPct=opts["islandPct"], islandErode=opts["islandErode"],
            )
        finally:
            qt.QApplication.restoreOverrideCursor()
        A = self.ui.activePctSpinBox.value / 100.0 * ref
        self.log("Active scan: blue wireframe brain + opaque tissue red >= %.0f%% "
                 "(%.1f counts), white >= %.0f%%. Wire mesh %d pts." %
                 (self.ui.activePctSpinBox.value, A, self.ui.activeWhitePctSpinBox.value,
                  result["meshPoints"]))
        self._frame3D(black=True)

    # ---- views / export --------------------------------------------------
    def onSnapView(self, name):
        node = self.currentVolumeNode()
        if node is not None:
            self.logic.snapToView(node, name)

    def onExportSurface(self):
        self._exportMontage("surface")

    def onExportActive(self):
        self._exportMontage("active")

    def _exportMontage(self, kind):
        node = self.currentVolumeNode()
        if node is None:
            return
        outdir = self.ui.outputDirButton.directory
        if not outdir:
            self.log("Set an output folder first (Views / export section).")
            return
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            # Regenerate the requested scan so the montage always shows it.
            if kind == "surface":
                self._generateSurface(live=False)
                montage, _ = self.logic.captureSixViews(node, outdir, "surface_views", "SURFACE")
            else:
                self.onGenerateActive()
                montage, _ = self.logic.captureSixViews(node, outdir, "active_views", "ACTIVE")
        except Exception as exc:
            qt.QApplication.restoreOverrideCursor()
            self.log("Montage export failed: %s" % exc)
            return
        qt.QApplication.restoreOverrideCursor()
        self.log("Saved 6 views + montage.png to: %s" % montage)
        try:
            qt.QDesktopServices.openUrl(qt.QUrl.fromLocalFile(montage))
        except Exception:
            pass

    # ---- regional quantification -----------------------------------------
    def _configureRegionTable(self):
        t = self.ui.regionTableWidget
        headers = ["Region", "Group", "% Ref", "Hot %", "Class", "Grade", "L-R %"]
        t.setColumnCount(len(headers))
        t.setHorizontalHeaderLabels(headers)
        try:
            t.horizontalHeader().setStretchLastSection(False)
            t.horizontalHeader().setSectionResizeMode(0, qt.QHeaderView.Stretch)
            for c in range(1, len(headers)):
                t.horizontalHeader().setSectionResizeMode(c, qt.QHeaderView.ResizeToContents)
        except Exception:
            pass
        t.verticalHeader().setVisible(False)

    def _buildRegionTree(self):
        """One checkable top-level item per lobe/group, with the group's regions
        as checkable children. Everything ticked by default."""
        t = self.ui.regionSelectTreeWidget
        t.clear()
        t.setColumnCount(1)
        self._updatingTree = True
        try:
            groupItems = {}
            order = []
            for reg in AMEN_REGIONS:
                g = reg["group"]
                if g not in groupItems:
                    gi = qt.QTreeWidgetItem(t, [g])
                    gi.setFlags(gi.flags() | qt.Qt.ItemIsUserCheckable)
                    gi.setCheckState(0, qt.Qt.Checked)
                    gi.setExpanded(False)
                    groupItems[g] = gi
                    order.append(g)
                label = reg["name"] + ("  - " + reg["note"] if reg["note"] else "")
                ci = qt.QTreeWidgetItem(groupItems[g], [label])
                ci.setFlags(ci.flags() | qt.Qt.ItemIsUserCheckable)
                ci.setCheckState(0, qt.Qt.Checked)
                ci.setData(0, qt.Qt.UserRole, reg["name"])
        finally:
            self._updatingTree = False

    def onRegionTreeItemChanged(self, item, _col):
        if self._updatingTree:
            return
        self._updatingTree = True
        try:
            if _qcount(item.childCount) > 0:
                st = item.checkState(0)
                if st != qt.Qt.PartiallyChecked:
                    for i in range(_qcount(item.childCount)):
                        item.child(i).setCheckState(0, st)
            else:
                parent = item.parent()
                if parent is not None:
                    n = _qcount(parent.childCount)
                    checked = sum(1 for i in range(n)
                                  if parent.child(i).checkState(0) == qt.Qt.Checked)
                    parent.setCheckState(0, qt.Qt.Checked if checked == n
                                         else (qt.Qt.Unchecked if checked == 0
                                               else qt.Qt.PartiallyChecked))
        finally:
            self._updatingTree = False

    def _setAllRegions(self, state):
        t = self.ui.regionSelectTreeWidget
        self._updatingTree = True
        try:
            for g in range(_qcount(t.topLevelItemCount)):
                gi = t.topLevelItem(g)
                gi.setCheckState(0, state)
                for c in range(_qcount(gi.childCount)):
                    gi.child(c).setCheckState(0, state)
        finally:
            self._updatingTree = False

    def _checkedRegionNames(self):
        names = []
        t = self.ui.regionSelectTreeWidget
        for g in range(_qcount(t.topLevelItemCount)):
            gi = t.topLevelItem(g)
            for c in range(_qcount(gi.childCount)):
                ci = gi.child(c)
                if ci.checkState(0) == qt.Qt.Checked:
                    names.append(ci.data(0, qt.Qt.UserRole))
        return names

    def _regionParams(self):
        return dict(
            metric=self.ui.regionMetricComboBox.currentText,
            refMode=self.ui.regionRefModeComboBox.currentText,
            hypoPct=self.ui.hypoPctSpinBox.value,
            hyperPct=self.ui.hyperPctSpinBox.value,
            floorPct=self.ui.regionFloorSpinBox.value,
            asymmetryPct=self.ui.asymmetryPctSpinBox.value,
            focalHyper=self.ui.focalHyperCheckBox.checked,
            minHotVoxels=int(self.ui.minHotVoxelsSpinBox.value),
            minHotFrac=self.ui.minHotFracSpinBox.value,
            hyperWhitePct=self.ui.hyperWhiteSpinBox.value,
            hyperTopPct=self.ui.hyperTopSpinBox.value,
            hypoStepPct=self.ui.hypoStepSpinBox.value,
        )

    def onCreateRegionScaffold(self):
        node = self.currentVolumeNode()
        if node is None:
            return
        names = self._checkedRegionNames()
        if not names:
            self.log("No regions selected - tick at least one lobe or region in "
                     "'Regions to include', then 'Create selected region ROIs'.")
            return
        if not self.ui.orientationConfirmedCheckBox.checked:
            self.log("WARNING: orientation not confirmed - the scaffold's L/R sides "
                     "follow the current display orientation.")
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            n = self.logic.createRegionRoiScaffold(node, names=names)
            removed = (self.logic.removeUnselectedRegionRois(names)
                       if self.ui.removeUncheckedCheckBox.checked else 0)
        finally:
            qt.QApplication.restoreOverrideCursor()
        msg = ("Created/repositioned %d of %d region Box ROIs under the '%s' folder "
               "(Data module)." % (n, len(AMEN_REGIONS), REGION_FOLDER_NAME))
        if removed:
            msg += " Removed %d unticked region ROI(s)." % removed
        msg += (" These are APPROXIMATE starting positions - drag each box onto its "
                "true structure, then 'Quantify regions'. The cerebellum reference "
                "ROI for '%s' is separate (use 'Create cerebellum ROI')." % REGION_REF_MEAN)
        self.log(msg)

    def onQuantifyRegions(self):
        node = self.currentVolumeNode()
        if node is None:
            return
        if not self.ui.orientationConfirmedCheckBox.checked:
            self.log("WARNING: orientation not confirmed - do not read laterality yet.")
        refRoi = self.ui.cerebellumRoiSelector.currentNode()
        p = self._regionParams()
        if p["refMode"] == REGION_REF_MEAN and refRoi is None:
            self.log("Select a cerebellum ROI (the reference) for '%s' mode, or "
                     "switch to '%s'." % (REGION_REF_MEAN, REGION_REF_MAX))
            return
        # robust whole-brain max as a fallback denominator when no cerebellar max set.
        try:
            stats, _ = self.logic.computeStatistics(node, self.ui.robustPercentileSpinBox.value)
            robust = stats["robustMax"]
        except Exception:
            robust = 0.0
        opts = self.preprocessOptions()
        qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
        try:
            rows, meta = self.logic.quantifyRegions(
                node, refRoi, self.cerebellarMax, robustReference=robust,
                maskToBrain=opts["keepLargestIsland"], islandPct=opts["islandPct"],
                islandErode=opts["islandErode"], **p)
        except Exception as exc:
            self.log("Regional quantification failed: %s" % exc)
            return
        finally:
            qt.QApplication.restoreOverrideCursor()
        if not rows:
            self.log("No region ROIs found. Click 'Create selected region ROIs' (or "
                     "place & rename a Box ROI per region), position them, then "
                     "'Quantify regions'.")
            return
        self._regionRows = rows
        if self.ui.colorRoiByClassCheckBox.checked:
            self.logic.colorRegionRois(rows)
        self._populateRegionTable(rows)
        self.log(self.logic.formatRegionReport(rows, meta), replace=True)

    def _populateRegionTable(self, rows):
        t = self.ui.regionTableWidget
        t.setRowCount(len(rows))
        for i, row in enumerate(rows):
            asym = "" if row["asym"] is None else ("%+.0f%%%s" % (
                row["asym"], " *" if row["asymFlag"] else ""))
            clsTxt = row["cls"] + (" (focal)" if row["focal"] else "")
            grd = gradeLabel(row["grade"]) if row["grade"] else ""
            cells = [row["name"], row["group"], "%.1f%%" % row["pct"],
                     "%.0f%% (%d)" % (row["hotFrac"], row["hotN"]), clsTxt, grd, asym]
            cr, cg, cb = rowDisplayColor(row)
            bg = qt.QColor(int(cr * 255), int(cg * 255), int(cb * 255))
            bg.setAlpha(90)
            for c, text in enumerate(cells):
                it = qt.QTableWidgetItem(text)
                it.setData(qt.Qt.UserRole, row["name"])
                if c in (2, 3, 5, 6):
                    it.setTextAlignment(qt.Qt.AlignRight | qt.Qt.AlignVCenter)
                if c in (4, 5):
                    it.setBackground(bg)
                if c == 3 and row["focal"]:
                    it.setBackground(qt.QColor(255, 120, 90, 130))
                if c == 6 and row["asymFlag"]:
                    it.setBackground(qt.QColor(255, 210, 80, 110))
                t.setItem(i, c, it)

    def onRegionRowDoubleClicked(self, rowIndex, _col):
        if rowIndex < 0 or rowIndex >= len(self._regionRows):
            return
        name = self._regionRows[rowIndex]["name"]
        roi = slicer.mrmlScene.GetFirstNodeByName(name)
        if roi and roi.IsA("vtkMRMLMarkupsROINode"):
            c = [0.0, 0.0, 0.0]
            roi.GetCenter(c)
            self.logic.jumpToRAS(c, showCrosshair=self.ui.showReticleCheckBox.checked)
            self.log("Jumped to %s." % name)

    def onExportRegionTsv(self):
        if not self._regionRows:
            self.log("Quantify regions first, then export.")
            return
        outdir = self.ui.outputDirButton.directory
        if not outdir:
            self.log("Set an output folder first (Views / export section).")
            return
        p = self._regionParams()
        meta = dict(metric=p["metric"],
                    reflabel="cerebellar mean" if p["refMode"] == REGION_REF_MEAN else "cerebellar max")
        full, abn = self.logic.exportRegionTsv(self._regionRows, meta, outdir)
        nAbn = sum(1 for r in self._regionRows if r["cls"] in ("Hypoactive", "Hyperactive"))
        self.log("Saved regional tables:\n  %s (all %d)\n  %s (%d hypo/hyper)" % (
            full, len(self._regionRows), abn, nAbn))
        try:
            qt.QDesktopServices.openUrl(qt.QUrl.fromLocalFile(outdir))
        except Exception:
            pass

    # ---- view framing ----------------------------------------------------
    def _frame3D(self, black=True):
        lm = slicer.app.layoutManager()
        lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
        w = lm.threeDWidget(0)
        if not w:
            return
        v = w.threeDView()
        vn = v.mrmlViewNode()
        # The 3D background is always black for both scans.
        vn.SetBackgroundColor(0, 0, 0)
        vn.SetBackgroundColor2(0, 0, 0)
        v.resetFocalPoint()
        v.resetCamera()
        try:
            cam = v.renderWindow().GetRenderers().GetFirstRenderer().GetActiveCamera()
            cam.Zoom(2.8)
            v.forceRender()
        except Exception:
            pass


#
# SPECTBrainRenderLogic
#
class SPECTBrainRenderLogic(ScriptedLoadableModuleLogic):
    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.surfaceDisplayNode = None   # VR surface display node
        self.surfaceModelNode = None     # marching-cubes surface model node
        self.activeDisplayNode = None

    # ---- statistics ------------------------------------------------------
    def computeStatistics(self, volumeNode, robustPercentile=99.5):
        import numpy as np
        arr = slicer.util.arrayFromVolume(volumeNode).astype(float)
        rawMax = float(arr.max()) if arr.size else 0.0
        nonzero = arr[arr > 0]
        provMax = float(np.percentile(nonzero, robustPercentile)) if nonzero.size else rawMax
        # Light mask just for robust mean / robust max reporting.
        binary = arr > 0.12 * provMax if provMax > 0 else arr > 0
        maskArr = self._largestConnectedComponent(binary)
        vals = arr[maskArr] if maskArr.any() else (nonzero if nonzero.size else arr.ravel())
        robustMax = float(np.percentile(vals, robustPercentile)) if vals.size else rawMax
        mean = float(vals.mean()) if vals.size else 0.0
        return dict(rawMax=rawMax, robustMax=robustMax, mean=mean,
                    robustPercentile=robustPercentile,
                    maskVoxels=int(maskArr.sum()), totalVoxels=int(arr.size)), maskArr

    def getReferenceValue(self, volumeNode, mode, cerebellarMax, stats=None):
        if mode == REF_CEREBELLUM and cerebellarMax and cerebellarMax > 0:
            return float(cerebellarMax)
        if stats is None:
            stats, _ = self.computeStatistics(volumeNode)
        if mode == REF_RAW_MAX:
            return stats["rawMax"]
        # REF_ROBUST_MAX, or cerebellum-but-unset -> robust whole-brain max.
        return stats["robustMax"]

    def _largestConnectedComponent(self, binary):
        import numpy as np
        if not binary.any():
            return binary
        try:
            from scipy import ndimage
            labeled, n = ndimage.label(binary)
            if n <= 1:
                return binary
            counts = np.bincount(labeled.ravel())
            counts[0] = 0
            return labeled == int(counts.argmax())
        except Exception:
            pass
        from vtk.util import numpy_support
        shape = binary.shape
        img = vtk.vtkImageData()
        img.SetDimensions(shape[2], shape[1], shape[0])
        flat = np.ascontiguousarray(binary.astype(np.uint8)).ravel(order="C")
        img.GetPointData().SetScalars(
            numpy_support.numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR))
        cc = vtk.vtkImageConnectivityFilter()
        cc.SetInputData(img)
        cc.SetExtractionModeToLargestRegion()
        cc.SetScalarRange(1, 1)
        cc.Update()
        return numpy_support.vtk_to_numpy(cc.GetOutput().GetPointData().GetScalars()).reshape(shape) > 0

    # ---- orientation -----------------------------------------------------
    def applyOrientationFlips(self, volumeNode, baselineMatrix, flipLR, flipAP, flipSI):
        flip = vtk.vtkMatrix4x4()
        flip.Identity()
        flip.SetElement(0, 0, -1.0 if flipLR else 1.0)
        flip.SetElement(1, 1, -1.0 if flipAP else 1.0)
        flip.SetElement(2, 2, -1.0 if flipSI else 1.0)
        result = vtk.vtkMatrix4x4()
        vtk.vtkMatrix4x4.Multiply4x4(flip, baselineMatrix, result)
        volumeNode.SetIJKToRASMatrix(result)
        volumeNode.Modified()

    # ---- brain mask (largest island) ------------------------------------
    def _largestIslandMask(self, volumeNode, thresholdValue, erode=0):
        """Boolean brain-mask array (same shape as arrayFromVolume): the largest
        connected component above thresholdValue, with an optional morphological
        opening (erode -> keep-largest -> dilate) to snap thin bridges so
        extracranial uptake touching the brain is dropped too. The single source
        of truth for both the render mask and the regional-quantification mask."""
        from vtk.util import numpy_support
        th = vtk.vtkImageThreshold()
        th.SetInputData(volumeNode.GetImageData())
        th.ThresholdByUpper(thresholdValue)
        th.SetInValue(1)
        th.SetOutValue(0)
        th.SetOutputScalarTypeToUnsignedChar()
        th.Update()
        binimg = th.GetOutput()
        if erode > 0:
            k = 2 * int(erode) + 1
            er = vtk.vtkImageDilateErode3D()
            er.SetInputData(binimg)
            er.SetDilateValue(0)
            er.SetErodeValue(1)
            er.SetKernelSize(k, k, k)
            er.Update()
            binimg = er.GetOutput()
        cc = vtk.vtkImageConnectivityFilter()
        cc.SetInputData(binimg)
        cc.SetExtractionModeToLargestRegion()
        cc.SetScalarRange(1, 1)
        cc.SetLabelModeToConstantValue()
        cc.SetLabelConstantValue(1)
        cc.SetLabelScalarTypeToUnsignedChar()
        cc.Update()
        out = cc.GetOutput()
        if erode > 0:
            k = 2 * int(erode) + 1
            di = vtk.vtkImageDilateErode3D()
            di.SetInputData(out)
            di.SetDilateValue(1)
            di.SetErodeValue(0)
            di.SetKernelSize(k, k, k)
            di.Update()
            out = di.GetOutput()
        m = numpy_support.vtk_to_numpy(out.GetPointData().GetScalars())
        shape = slicer.util.arrayFromVolume(volumeNode).shape
        return m.reshape(shape) > 0

    def largestIslandVolume(self, volumeNode, thresholdValue, erode=0,
                            maskedName=MASKED_VOLUME_NAME):
        """Clone `maskedName` = largest connected component above thresholdValue.
        Optional morphological opening (erode -> keep-largest -> dilate) snaps
        thin bridges so extracranial uptake touching the brain is dropped too.
        Ported from the clinic's surface_scan.py / active_scan.py."""
        import numpy as np
        maskBool = self._largestIslandMask(volumeNode, thresholdValue, erode)
        masked = slicer.mrmlScene.GetFirstNodeByName(maskedName)
        if not masked or not masked.IsA("vtkMRMLScalarVolumeNode"):
            masked = slicer.modules.volumes.logic().CloneVolume(
                slicer.mrmlScene, volumeNode, maskedName)
        src = slicer.util.arrayFromVolume(volumeNode)
        ma = slicer.util.arrayFromVolume(masked)
        ma[:] = np.where(maskBool, src, 0)
        slicer.util.arrayFromVolumeModified(masked)
        # Keep masked aligned with current orientation of the source.
        ijk = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(ijk)
        masked.SetIJKToRASMatrix(ijk)
        return masked

    def _hideOtherVR(self, keepVol):
        for vrdn in slicer.util.getNodesByClass("vtkMRMLVolumeRenderingDisplayNode"):
            n = vrdn.GetVolumeNode()
            if n and n is not keepVol:
                vrdn.SetVisibility(False)

    def _smoothVolumeInPlace(self, volumeNode, sigma):
        """Gaussian-smooth a volume's scalars in place. Only ever called on the
        render clone (never the original), so the surface reads solid instead of
        showing noisy squiggles through a semi-transparent render."""
        if not sigma or sigma <= 0:
            return
        g = vtk.vtkImageGaussianSmooth()
        g.SetInputData(volumeNode.GetImageData())
        g.SetStandardDeviations(sigma, sigma, sigma)
        g.SetRadiusFactors(2, 2, 2)
        g.Update()
        volumeNode.GetImageData().DeepCopy(g.GetOutput())
        volumeNode.Modified()

    def _resolveRenderVolume(self, volumeNode, keepLargestIsland, islandPct, islandErode,
                             reference, renderName=MASKED_VOLUME_NAME):
        """Returns the per-scan masked render clone named `renderName` (never the
        original, so it can be smoothed safely). Each scan passes its own name so
        the surface and active renderings stay fully independent. Island mask on
        -> drop scalp/sinus/salivary uptake; off -> a plain copy of the source."""
        if keepLargestIsland:
            thrv = (reference if reference > 0 else volumeNode.GetImageData().GetScalarRange()[1]) * islandPct / 100.0
            masked = self.largestIslandVolume(volumeNode, thrv, islandErode, maskedName=renderName)
        else:
            masked = slicer.mrmlScene.GetFirstNodeByName(renderName)
            if not masked or not masked.IsA("vtkMRMLScalarVolumeNode"):
                masked = slicer.modules.volumes.logic().CloneVolume(
                    slicer.mrmlScene, volumeNode, renderName)
            slicer.util.arrayFromVolume(masked)[:] = slicer.util.arrayFromVolume(volumeNode)
            slicer.util.arrayFromVolumeModified(masked)
            ijk = vtk.vtkMatrix4x4()
            volumeNode.GetIJKToRASMatrix(ijk)
            masked.SetIJKToRASMatrix(ijk)
        self._hideOtherVR(masked)
        return masked

    def _volumeRenderingDisplayNode(self, renderVol):
        vrLogic = slicer.modules.volumerendering.logic()
        try:
            dn = vrLogic.GetFirstVolumeRenderingDisplayNode(renderVol)
        except Exception:
            dn = None
        if dn is None:
            dn = vrLogic.CreateDefaultVolumeRenderingNodes(renderVol)
        dn.SetVisibility(True)
        dn.SetCroppingEnabled(False)
        return dn

    def applySurfaceColor(self, volumeProperty, reference, scheme=SURFACE_SCHEME_HOLE):
        points = HOLE_COLOR_FRACTIONS if scheme == SURFACE_SCHEME_HOLE else AMEN_COLOR_FRACTIONS
        c = volumeProperty.GetRGBTransferFunction(0)
        c.RemoveAllPoints()
        for frac, (r, g, b) in points:
            c.AddRGBPoint(frac * reference, r, g, b)
        return len(points)

    # ---- surface scan (volume render) ------------------------------------
    def generateSurfaceVR(self, volumeNode, surfacePct, reference, opacity=0.9,
                          smoothingSigma=1.0, colorScheme=SURFACE_SCHEME_HOLE,
                          keepLargestIsland=True, islandPct=30.0, islandErode=2,
                          keepPresetOpacity=False):
        if reference <= 0:
            raise ValueError("Reference (100%) must be > 0.")
        renderVol = self._resolveRenderVolume(volumeNode, keepLargestIsland, islandPct,
                                              islandErode, reference, renderName=SURFACE_RENDER_NAME)
        # Smooth the render clone so the surface reads solid, not squiggly.
        self._smoothVolumeInPlace(renderVol, smoothingSigma)
        lo, hi = renderVol.GetImageData().GetScalarRange()

        # Hide active-scan models and the mesh surface (the VR is now the surface).
        for nm in (WIRE_MODEL_NAME, WIRE_NODES_NAME, SURFACE_MODEL_NAME):
            mn = slicer.mrmlScene.GetFirstNodeByName(nm)
            if mn and mn.GetDisplayNode():
                mn.GetDisplayNode().SetVisibility(False)

        dn = self._volumeRenderingDisplayNode(renderVol)
        self.surfaceDisplayNode = dn
        vprop = dn.GetVolumePropertyNode().GetVolumeProperty()

        colorPoints = self.applySurfaceColor(vprop, reference, colorScheme)
        peak = max(0.05, min(1.0, opacity))
        T = reference * surfacePct / 100.0
        so = vprop.GetScalarOpacity()
        so.RemoveAllPoints()
        if keepPresetOpacity:
            # Soft preset look (no hard defect threshold).
            so.AddPoint(lo, 0.0)
            so.AddPoint(reference * 0.38, 0.0)
            so.AddPoint(reference * 0.72, 0.49 * peak)
            so.AddPoint(reference, peak)
            so.AddPoint(hi, min(1.0, peak + 0.07))
        else:
            # Sharp threshold so hypoperfusion reads as dents / holes; the ramp
            # reaches the (high) peak opacity just above threshold = solid shell.
            so.AddPoint(lo, 0.0)
            so.AddPoint(T * 0.95, 0.0)
            so.AddPoint(T, peak)
            so.AddPoint(hi, min(1.0, peak + 0.07))
        go = vprop.GetGradientOpacity()
        go.RemoveAllPoints()
        go.AddPoint(lo, 1.0)
        go.AddPoint(hi, 1.0)
        vprop.ShadeOn()
        vprop.SetDiffuse(0.65)
        vprop.SetAmbient(0.35)
        vprop.SetSpecular(0.25)  # lower specular = fewer shiny squiggle highlights
        vprop.SetSpecularPower(20)
        vprop.SetInterpolationTypeToLinear()
        return dict(displayNode=dn, renderVolume=renderVol, threshold=T,
                    colorPoints=colorPoints)

    # ---- surface scan (marching-cubes isosurface) ------------------------
    def generateSurfaceMesh(self, volumeNode, surfacePct, reference, opacity=0.95,
                            smoothingSigma=1.0, colorScheme=SURFACE_SCHEME_SOLID,
                            keepLargestIsland=True, islandPct=30.0, islandErode=2,
                            smoothIterations=18, passBand=0.10, decimate=0.0):
        """Amen-style surface as a SHADED MARCHING-CUBES ISOSURFACE at the chosen
        perfusion threshold. The isovalue separates above-threshold (inside) from
        below-threshold (outside) voxels, so low-perfusion regions fall inward as
        crisp solid holes / dents - the authentic Amen mechanism (and the research-
        recommended match for the current Amen look). The mesh is windowed-sinc
        smoothed and Gouraud-shaded with a matte tissue material so it reads as a
        solid sculpted brain, not a translucent cloud. With a color scheme other
        than 'solid', each vertex is tinted by the underlying perfusion."""
        if reference <= 0:
            raise ValueError("Reference (100%) must be > 0.")
        renderVol = self._resolveRenderVolume(volumeNode, keepLargestIsland, islandPct,
                                              islandErode, reference, renderName=SURFACE_RENDER_NAME)
        # Pre-smooth the render clone so the isosurface is not stair-stepped.
        self._smoothVolumeInPlace(renderVol, smoothingSigma)
        T = reference * surfacePct / 100.0

        # Only one surface at a time: hide the VR surface and the active models.
        if self.surfaceDisplayNode and slicer.mrmlScene.IsNodePresent(self.surfaceDisplayNode):
            self.surfaceDisplayNode.SetVisibility(False)
        for nm in (WIRE_MODEL_NAME, WIRE_NODES_NAME):
            mn = slicer.mrmlScene.GetFirstNodeByName(nm)
            if mn and mn.GetDisplayNode():
                mn.GetDisplayNode().SetVisibility(False)

        mc = vtk.vtkFlyingEdges3D()
        mc.SetInputData(renderVol.GetImageData())
        mc.SetValue(0, T)
        mc.ComputeNormalsOff()  # recomputed after smoothing + RAS transform
        mc.Update()
        conn = vtk.vtkPolyDataConnectivityFilter()
        conn.SetInputConnection(mc.GetOutputPort())
        conn.SetExtractionModeToLargestRegion()  # the brain, not stray islands
        conn.Update()
        sm = vtk.vtkWindowedSincPolyDataFilter()
        sm.SetInputConnection(conn.GetOutputPort())
        sm.SetNumberOfIterations(int(smoothIterations))
        sm.SetPassBand(passBand)         # ~0.1: smooth without shrinking the dents
        sm.NormalizeCoordinatesOn()
        sm.NonManifoldSmoothingOn()
        sm.Update()
        last = sm
        if decimate and decimate > 0:
            de = vtk.vtkDecimatePro()
            de.SetInputConnection(last.GetOutputPort())
            de.SetTargetReduction(decimate)
            de.PreserveTopologyOn()
            de.Update()
            last = de
        ijkPoly = last.GetOutput()  # still in IJK index space

        # Perfusion tint (hole / perfusion-rainbow) samples the underlying volume,
        # so it is done in IJK space - where the mesh and the image data share
        # index coordinates - BEFORE the RAS transform. The positional spectrum is
        # applied later, on the RAS mesh, since it keys off anatomical position.
        colored, colorPoints = False, 0
        if colorScheme in (SURFACE_SCHEME_HOLE, SURFACE_SCHEME_RAINBOW):
            colorPoints = self._colorMeshByPerfusion(ijkPoly, renderVol.GetImageData(),
                                                     reference, colorScheme)
            colored = colorPoints > 0

        i2r = vtk.vtkMatrix4x4()
        renderVol.GetIJKToRASMatrix(i2r)
        tr = vtk.vtkTransform()
        tr.SetMatrix(i2r)
        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetTransform(tr)
        tf.SetInputData(ijkPoly)
        tf.Update()
        nrm = vtk.vtkPolyDataNormals()
        nrm.SetInputConnection(tf.GetOutputPort())
        nrm.SetFeatureAngle(60)
        nrm.SplittingOff()      # keep vertex count 1:1 so the RGB array survives
        nrm.ConsistencyOn()
        nrm.Update()
        poly = nrm.GetOutput()

        # Positional spectrum keys off anatomical (RAS) position, so it is applied
        # here on the transformed mesh.
        if colorScheme == SURFACE_SCHEME_SPECTRUM:
            colorPoints = self._colorMeshBySpectrum(poly)
            colored = colorPoints > 0

        model = slicer.mrmlScene.GetFirstNodeByName(SURFACE_MODEL_NAME)
        if not model or not model.IsA("vtkMRMLModelNode"):
            model = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", SURFACE_MODEL_NAME)
        model.SetAndObservePolyData(poly)
        if not model.GetDisplayNode():
            model.CreateDefaultDisplayNodes()
        d = model.GetDisplayNode()
        d.SetVisibility(True)
        d.SetOpacity(max(0.05, min(1.0, opacity)))
        d.SetColor(*SURFACE_SOLID_COLOR)
        d.SetScalarVisibility(colored)
        if colored:
            d.SetActiveScalarName("SurfaceRGB")
            try:
                d.SetScalarRangeFlag(slicer.vtkMRMLDisplayNode.UseDirectMapping)
            except Exception:
                pass
        # Strong directional lighting carries the shape: diffuse-dominant so the
        # surface brightens where it faces the light and the dents drop into
        # shadow, low specular (matte, not glossy), smooth Gouraud shading (needs
        # the point normals computed above, else faceted). Tuned to the reference
        # contrast (highlights ~0.9, recesses ~0.3 of peak).
        d.SetAmbient(0.28)
        d.SetDiffuse(0.82)
        d.SetSpecular(0.06)
        d.SetPower(8)
        d.SetLighting(True)
        try:
            d.SetInterpolation(1)  # 1 = Gouraud
        except Exception:
            pass
        try:
            d.SetVisibility2D(False)
        except Exception:
            pass
        self.surfaceModelNode = model
        return dict(modelNode=model, renderVolume=renderVol, threshold=T,
                    points=int(poly.GetNumberOfPoints()), colored=colored,
                    colorPoints=colorPoints)

    def _attachDirectRGB(self, poly, rgb):
        """Attach an (n,3) uint8 array as the mesh's active 'SurfaceRGB' scalars;
        Slicer renders a 3-component active scalar as direct RGB."""
        import numpy as np
        from vtk.util import numpy_support
        arr = numpy_support.numpy_to_vtk(np.ascontiguousarray(rgb.astype(np.uint8)),
                                         deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
        arr.SetName("SurfaceRGB")
        poly.GetPointData().SetScalars(arr)
        poly.GetPointData().SetActiveScalars("SurfaceRGB")

    def _colorMeshByPerfusion(self, poly, imageData, reference, scheme):
        """Tint each vertex by the underlying perfusion (sampled in IJK index
        space) through the chosen perfusion scheme, as a direct-RGB array.
        Returns the number of color control points."""
        import numpy as np
        from vtk.util import numpy_support
        points = HOLE_COLOR_FRACTIONS if scheme == SURFACE_SCHEME_HOLE else AMEN_COLOR_FRACTIONS
        ctf = vtk.vtkColorTransferFunction()
        for frac, (r, g, b) in points:
            ctf.AddRGBPoint(frac * reference, r, g, b)
        probe = vtk.vtkProbeFilter()
        probe.SetInputData(poly)
        probe.SetSourceData(imageData)
        probe.Update()
        sampled = probe.GetOutput().GetPointData().GetScalars()
        n = poly.GetNumberOfPoints()
        rgb = np.zeros((n, 3), dtype=np.uint8)
        if sampled is not None and n:
            vals = numpy_support.vtk_to_numpy(sampled).astype(float).ravel()
            for idx in range(n):
                c = ctf.GetColor(float(vals[idx]))
                rgb[idx] = (int(c[0] * 255), int(c[1] * 255), int(c[2] * 255))
        self._attachDirectRGB(poly, rgb)
        return len(points)

    def _colorMeshBySpectrum(self, poly, axis=SURFACE_SPECTRUM_AXIS,
                             invert=SURFACE_SPECTRUM_INVERT):
        """Tint each vertex by a smooth spectrum keyed to its anatomical (RAS)
        position along `axis` (default S/I), as a direct-RGB array - reproducing
        the Amen surface's positional rainbow. Lighting (set on the display node)
        carries the shape on top of this. Returns the number of spectrum stops."""
        import numpy as np
        from vtk.util import numpy_support
        n = poly.GetNumberOfPoints()
        if not n:
            return 0
        pts = numpy_support.vtk_to_numpy(poly.GetPoints().GetData()).astype(float)
        coord = pts[:, axis]
        lo, hi = float(coord.min()), float(coord.max())
        s = (coord - lo) / (hi - lo) if hi > lo else np.zeros(n)
        if invert:
            s = 1.0 - s
        ctf = vtk.vtkColorTransferFunction()
        for frac, (r, g, b) in SURFACE_SPECTRUM:
            ctf.AddRGBPoint(frac, r, g, b)
        rgb = np.zeros((n, 3), dtype=np.uint8)
        for idx in range(n):
            c = ctf.GetColor(float(s[idx]))
            rgb[idx] = (int(c[0] * 255), int(c[1] * 255), int(c[2] * 255))
        self._attachDirectRGB(poly, rgb)
        return len(SURFACE_SPECTRUM)

    # ---- active scan (wireframe + hot volume) ----------------------------
    def generateActive(self, volumeNode, reference, activePct=80.0, whitePct=92.0, meshPct=15.0,
                       nodeRadius=0.75, decimate=0.85, wireOpacity=0.5,
                       meshColor=(0.20, 0.45, 1.0), nodeColor=(0.55, 0.75, 1.0),
                       meshLineWidth=1.0, keepLargestIsland=True, islandPct=30.0, islandErode=2):
        if reference <= 0:
            raise ValueError("Reference (100%) must be > 0.")
        renderVol = self._resolveRenderVolume(volumeNode, keepLargestIsland, islandPct,
                                              islandErode, reference, renderName=ACTIVE_RENDER_NAME)
        lo, hi = renderVol.GetImageData().GetScalarRange()
        A = reference * activePct / 100.0       # red onset (>= activePct % of cereb max)
        Wpt = reference * whitePct / 100.0       # white onset (>= whitePct %)
        top = max(reference, hi)                 # cerebellum (100%) is the ceiling

        # Hide the surface-scan model fallback if present.
        mn = slicer.mrmlScene.GetFirstNodeByName(SURFACE_MODEL_NAME)
        if mn and mn.GetDisplayNode():
            mn.GetDisplayNode().SetVisibility(False)

        meshPoints = self._blueBrainMesh(
            renderVol, reference * meshPct / 100.0, decimate, nodeRadius,
            meshColor, nodeColor, meshLineWidth, wireOpacity)

        dn = self._volumeRenderingDisplayNode(renderVol)
        self.activeDisplayNode = dn
        vp = dn.GetVolumePropertyNode().GetVolumeProperty()
        so = vp.GetScalarOpacity()
        so.RemoveAllPoints()
        so.AddPoint(0, 0.0)
        so.AddPoint(A * 0.97, 0.0)
        so.AddPoint(A, 0.9)
        so.AddPoint(top, 1.0)
        # red (>= activePct) ramps to white (>= whitePct), white to the ceiling.
        c = vp.GetRGBTransferFunction(0)
        c.RemoveAllPoints()
        c.AddRGBPoint(0, 1.0, 0.0, 0.0)
        c.AddRGBPoint(A, 1.0, 0.0, 0.0)
        c.AddRGBPoint(max(Wpt, A + 1e-3), 1.0, 1.0, 1.0)
        c.AddRGBPoint(top, 1.0, 1.0, 1.0)
        go = vp.GetGradientOpacity()
        go.RemoveAllPoints()
        go.AddPoint(0, 1.0)
        go.AddPoint(top, 1.0)
        vp.ShadeOn()
        vp.SetAmbient(0.35)
        vp.SetDiffuse(0.65)
        vp.SetSpecular(0.30)
        vp.SetSpecularPower(25)
        vp.SetInterpolationTypeToLinear()
        return dict(displayNode=dn, renderVolume=renderVol, active=A, meshPoints=meshPoints)

    def _blueBrainMesh(self, volumeNode, iso, decimate, nodeRadius, meshColor,
                       nodeColor, meshLineWidth, wireOpacity):
        g = vtk.vtkImageGaussianSmooth()
        g.SetInputData(volumeNode.GetImageData())
        g.SetStandardDeviations(1.5, 1.5, 1.5)
        g.SetRadiusFactors(2, 2, 2)
        g.Update()
        mc = vtk.vtkFlyingEdges3D()
        mc.SetInputConnection(g.GetOutputPort())
        mc.SetValue(0, iso)
        mc.ComputeNormalsOff()
        mc.Update()
        conn = vtk.vtkPolyDataConnectivityFilter()
        conn.SetInputConnection(mc.GetOutputPort())
        conn.SetExtractionModeToLargestRegion()
        conn.Update()
        sm = vtk.vtkWindowedSincPolyDataFilter()
        sm.SetInputConnection(conn.GetOutputPort())
        sm.SetNumberOfIterations(15)
        sm.SetPassBand(0.12)
        sm.NormalizeCoordinatesOn()
        sm.Update()
        up = sm.GetOutputPort()
        if decimate and decimate > 0:
            de = vtk.vtkDecimatePro()
            de.SetInputConnection(up)
            de.SetTargetReduction(decimate)
            de.PreserveTopologyOn()
            de.Update()
            up = de.GetOutputPort()
        i2r = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(i2r)
        t = vtk.vtkTransform()
        t.SetMatrix(i2r)
        tf = vtk.vtkTransformPolyDataFilter()
        tf.SetTransform(t)
        tf.SetInputConnection(up)
        tf.Update()
        nrm = vtk.vtkPolyDataNormals()
        nrm.SetInputConnection(tf.GetOutputPort())
        nrm.SetFeatureAngle(60)
        nrm.Update()
        cl = vtk.vtkCleanPolyData()
        cl.SetInputConnection(nrm.GetOutputPort())
        cl.Update()
        poly = cl.GetOutput()

        m = slicer.mrmlScene.GetFirstNodeByName(WIRE_MODEL_NAME)
        if not m or not m.IsA("vtkMRMLModelNode"):
            m = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", WIRE_MODEL_NAME)
        m.SetAndObservePolyData(poly)
        if not m.GetDisplayNode():
            m.CreateDefaultDisplayNodes()
        d = m.GetDisplayNode()
        d.SetColor(*meshColor)
        d.SetScalarVisibility(False)
        d.SetRepresentation(1)  # 1 = wireframe
        d.SetLineWidth(meshLineWidth)
        d.SetOpacity(wireOpacity)
        d.SetLighting(True)
        d.SetVisibility(True)
        try:
            d.SetVisibility2D(False)
        except Exception:
            pass

        # Optional solid spheres at each wire vertex. Radius <= 0 -> no nodes
        # (just the bare wireframe); hide any nodes from a previous run.
        existingNodes = slicer.mrmlScene.GetFirstNodeByName(WIRE_NODES_NAME)
        if nodeRadius and nodeRadius > 0:
            sph = vtk.vtkSphereSource()
            sph.SetRadius(nodeRadius)
            sph.SetThetaResolution(8)
            sph.SetPhiResolution(8)
            gl = vtk.vtkGlyph3D()
            gl.SetInputData(poly)
            gl.SetSourceConnection(sph.GetOutputPort())
            gl.SetScaleModeToDataScalingOff()
            gl.Update()
            mnode = existingNodes
            if not mnode or not mnode.IsA("vtkMRMLModelNode"):
                mnode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLModelNode", WIRE_NODES_NAME)
            mnode.SetAndObservePolyData(gl.GetOutput())
            if not mnode.GetDisplayNode():
                mnode.CreateDefaultDisplayNodes()
            dn2 = mnode.GetDisplayNode()
            dn2.SetColor(*nodeColor)
            dn2.SetScalarVisibility(False)
            dn2.SetRepresentation(2)  # 2 = surface (solid spheres)
            dn2.SetVisibility(True)
            try:
                dn2.SetVisibility2D(False)
            except Exception:
                pass
        elif existingNodes and existingNodes.GetDisplayNode():
            existingNodes.GetDisplayNode().SetVisibility(False)
        return int(poly.GetNumberOfPoints())

    # ---- cerebellar max --------------------------------------------------
    def createCerebellumRoi(self, volumeNode):
        """Place a Box ROI roughly over the posterior fossa (posterior-inferior)."""
        import numpy as np
        bounds = [0.0] * 6
        volumeNode.GetRASBounds(bounds)
        cx = (bounds[0] + bounds[1]) / 2.0
        # Posterior third in A, inferior third in S.
        ay = bounds[2] + 0.22 * (bounds[3] - bounds[2])
        sz = bounds[4] + 0.28 * (bounds[5] - bounds[4])
        roi = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", "R")
        roi.SetCenter(cx, ay, sz)
        size = [0.45 * (bounds[1] - bounds[0]),
                0.30 * (bounds[3] - bounds[2]),
                0.25 * (bounds[5] - bounds[4])]
        roi.SetSize(*size)
        roi.CreateDefaultDisplayNodes()
        return roi

    def cerebellarMaxCandidates(self, volumeNode, roiNode, topN=10):
        """Return the topN hottest voxels inside the (axis-aligned) Box ROI as a
        list of {rank, value, ijk, ras}. Ported from cerebellar_max.py."""
        import numpy as np
        import itertools
        arr = slicer.util.arrayFromVolume(volumeNode)
        rasToIjk = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(rasToIjk)
        ijkToRas = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(ijkToRas)
        bounds = [0.0] * 6
        roiNode.GetRASBounds(bounds)
        if not (bounds[1] > bounds[0] and bounds[3] > bounds[2] and bounds[5] > bounds[4]):
            raise RuntimeError("ROI bounds look degenerate - is the Box ROI defined?")
        ijk = np.array([rasToIjk.MultiplyPoint([x, y, z, 1.0])[:3]
                        for x, y, z in itertools.product(bounds[0:2], bounds[2:4], bounds[4:6])])
        nk, nj, ni = arr.shape
        i0, j0, k0 = np.floor(ijk.min(0)).astype(int)
        i1, j1, k1 = np.ceil(ijk.max(0)).astype(int)
        i0, i1 = max(i0, 0), min(i1, ni - 1)
        j0, j1 = max(j0, 0), min(j1, nj - 1)
        k0, k1 = max(k0, 0), min(k1, nk - 1)
        sub = arr[k0:k1 + 1, j0:j1 + 1, i0:i1 + 1]
        if sub.size == 0:
            return []
        flat = sub.ravel()
        N = int(min(topN, flat.size))
        order = np.argpartition(flat, -N)[-N:]
        order = order[np.argsort(flat[order])[::-1]]
        out = []
        for rank, fi in enumerate(order, 1):
            kk, jj, ii = np.unravel_index(fi, sub.shape)
            gi, gj, gk = int(ii + i0), int(jj + j0), int(kk + k0)
            r = ijkToRas.MultiplyPoint([gi, gj, gk, 1.0])[:3]
            out.append(dict(rank=rank, value=float(sub[kk, jj, ii]),
                            ijk=(gi, gj, gk), ras=(float(r[0]), float(r[1]), float(r[2]))))
        return out

    def jumpToRAS(self, ras, showCrosshair=True):
        lm = slicer.app.layoutManager()
        if lm:
            for nm in ("Red", "Yellow", "Green"):
                w = lm.sliceWidget(nm)
                if w:
                    w.mrmlSliceNode().JumpSliceByCentering(ras[0], ras[1], ras[2])
        for ch in slicer.util.getNodesByClass("vtkMRMLCrosshairNode"):
            if showCrosshair:
                ch.SetCrosshairRAS(ras)
                ch.SetCrosshairMode(slicer.vtkMRMLCrosshairNode.ShowBasic)
            else:
                ch.SetCrosshairMode(slicer.vtkMRMLCrosshairNode.NoCrosshair)

    def setCrosshairVisible(self, visible):
        mode = (slicer.vtkMRMLCrosshairNode.ShowBasic if visible
                else slicer.vtkMRMLCrosshairNode.NoCrosshair)
        for ch in slicer.util.getNodesByClass("vtkMRMLCrosshairNode"):
            ch.SetCrosshairMode(mode)

    # ---- regional quantification (ports regional_table.py) ---------------
    def regionRoiStats(self, roiNode, arr, rasToIjk, floor, hotThreshold=None):
        """Tissue (> floor) stats inside an ROI's axis-aligned bounding box, as
        dict(mean, peak, nvox, hotN, hotFrac, hotMean) - or None if no tissue.
        When hotThreshold (counts) is given, also counts/averages the "hot"
        voxels >= it (the focal-hyperactivity test). Same AABB approach as the
        clinic's regional_table._roi_stats / cerebellar_max."""
        import numpy as np
        import itertools
        bounds = [0.0] * 6
        roiNode.GetRASBounds(bounds)
        if not (bounds[1] > bounds[0] and bounds[3] > bounds[2] and bounds[5] > bounds[4]):
            return None
        ijk = np.array([rasToIjk.MultiplyPoint([x, y, z, 1.0])[:3]
                        for x, y, z in itertools.product(bounds[0:2], bounds[2:4], bounds[4:6])])
        nk, nj, ni = arr.shape
        i0, j0, k0 = np.floor(ijk.min(0)).astype(int)
        i1, j1, k1 = np.ceil(ijk.max(0)).astype(int)
        i0, i1 = max(i0, 0), min(i1, ni - 1)
        j0, j1 = max(j0, 0), min(j1, nj - 1)
        k0, k1 = max(k0, 0), min(k1, nk - 1)
        sub = arr[k0:k1 + 1, j0:j1 + 1, i0:i1 + 1]
        tissue = sub[sub > floor]
        if tissue.size == 0:
            return None
        res = dict(mean=float(tissue.mean()), peak=float(tissue.max()),
                   nvox=int(tissue.size), hotN=0, hotFrac=0.0, hotMean=None)
        if hotThreshold is not None:
            hot = tissue[tissue >= hotThreshold]
            res["hotN"] = int(hot.size)
            res["hotFrac"] = float(hot.size) / float(tissue.size)
            res["hotMean"] = float(hot.mean()) if hot.size else None
        return res

    def cerebellumMean(self, volumeNode, roiNode, floor):
        """Mean counts inside the cerebellum reference ROI (for '% of cerebellar
        mean' mode). Returns 0.0 if unavailable."""
        if roiNode is None:
            return 0.0
        arr = slicer.util.arrayFromVolume(volumeNode)
        r2i = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(r2i)
        res = self.regionRoiStats(roiNode, arr, r2i, floor)
        return float(res["mean"]) if res else 0.0

    def createRegionRoiScaffold(self, volumeNode, names=None, sizeScale=1.0):
        """Drop (or reposition) one approximately-placed, named Box ROI for each
        selected AMEN_REGIONS entry, collected under a subject-hierarchy folder.
        names=None creates all regions; otherwise only those whose name is listed.
        Returns the count created/updated. The clinician drags each box onto the
        true structure; positions are only a starting scaffold."""
        selected = None if names is None else set(names)
        bounds = [0.0] * 6
        volumeNode.GetRASBounds(bounds)
        ex = (bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
        sh = slicer.mrmlScene.GetSubjectHierarchyNode()
        folderItem = 0
        if sh is not None:
            folderItem = sh.GetItemChildWithName(sh.GetSceneItemID(), REGION_FOLDER_NAME)
            if not folderItem:
                folderItem = sh.CreateFolderItem(sh.GetSceneItemID(), REGION_FOLDER_NAME)
        n = 0
        for reg in AMEN_REGIONS:
            if selected is not None and reg["name"] not in selected:
                continue
            cx = bounds[0] + reg["fx"] * ex[0]
            cy = bounds[2] + reg["fy"] * ex[1]
            cz = bounds[4] + reg["fz"] * ex[2]
            frac = GROUP_BOX_FRAC.get(reg["group"], DEFAULT_BOX_FRAC) * sizeScale
            size = (frac * ex[0], frac * ex[1], frac * ex[2])
            roi = slicer.mrmlScene.GetFirstNodeByName(reg["name"])
            if not roi or not roi.IsA("vtkMRMLMarkupsROINode"):
                roi = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", reg["name"])
                roi.CreateDefaultDisplayNodes()
            roi.SetCenter(cx, cy, cz)
            roi.SetSize(*size)
            dn = roi.GetDisplayNode()
            if dn:
                dn.SetHandlesInteractive(True)
                try:
                    dn.SetFillOpacity(0.15)
                except Exception:
                    pass
            if sh is not None and folderItem:
                sh.SetItemParent(sh.GetItemByDataNode(roi), folderItem)
            n += 1
        return n

    def removeUnselectedRegionRois(self, keepNames):
        """Delete any taxonomy-named region Box ROI that is NOT in keepNames, so
        the scene holds exactly the selection. Only touches AMEN_REGIONS names -
        never the cerebellum reference or any custom ROI. Returns count removed."""
        keep = set(keepNames)
        removed = 0
        for reg in AMEN_REGIONS:
            if reg["name"] in keep:
                continue
            node = slicer.mrmlScene.GetFirstNodeByName(reg["name"])
            if node and node.IsA("vtkMRMLMarkupsROINode"):
                slicer.mrmlScene.RemoveNode(node)
                removed += 1
        return removed

    def quantifyRegions(self, volumeNode, referenceRoiNode, cerebellarMax,
                        metric=REGION_METRIC_MEAN, refMode=REGION_REF_MAX,
                        hypoPct=DEFAULT_HYPO_PCT, hyperPct=DEFAULT_HYPER_PCT,
                        floorPct=DEFAULT_REGION_FLOOR_PCT,
                        asymmetryPct=DEFAULT_ASYMMETRY_PCT, robustReference=0.0,
                        focalHyper=True, minHotVoxels=DEFAULT_MIN_HOT_VOXELS,
                        minHotFrac=DEFAULT_MIN_HOT_FRAC_PCT,
                        maskToBrain=True, islandPct=30.0, islandErode=2,
                        hyperWhitePct=DEFAULT_HYPER_WHITE_PCT,
                        hyperTopPct=DEFAULT_HYPER_TOP_PCT,
                        hypoStepPct=DEFAULT_HYPO_STEP_PCT):
        """Quantify every Markups ROI (except the cerebellum reference) against
        the cerebellar reference and classify hypo/normal/hyper. Returns
        (rows, meta). Each row dict: name, group, side, base, note, value, pct,
        cls, focal, nvox, hotN, hotFrac, hotMeanPct, asym, asymFlag.

        With maskToBrain on (default), stats are computed on the SAME largest-
        island brain mask the surface/active scans render from, so extracranial
        uptake (scalp / sinus / salivary / venous) that is masked OUT of the
        render is not counted - the table's hot voxels then match the red voxels.

        HYPOACTIVE uses the region MEAN (it is satisfying for low-perfusion
        regions). HYPERACTIVE additionally uses the FOCAL test: a region is
        hyperactive if it holds a hot focus - at least `minHotVoxels` voxels OR
        `minHotFrac` % of the ROI at/above hyperPct % of the cerebellar max - so
        a true hot focus is not diluted to "normal" by surrounding normal tissue
        in a loose box. Mirrors the clinic's old mask-at-85%-and-average read.
        Ports the clinic's regional_table()."""
        import numpy as np
        r2i = vtk.vtkMatrix4x4()
        volumeNode.GetRASToIJKMatrix(r2i)

        # Floor is always a fraction of the cerebellar MAX (CSF/background cut).
        floorBase = cerebellarMax if cerebellarMax and cerebellarMax > 0 else robustReference
        floor = floorBase * floorPct / 100.0

        # Brain-mask the stats volume so only brain voxels are read (matches the
        # render). Work on a COPY - never mutate the source volume's array.
        maskBase = cerebellarMax if cerebellarMax and cerebellarMax > 0 else robustReference
        masked = bool(maskToBrain and maskBase and maskBase > 0)
        if masked:
            brainMask = self._largestIslandMask(
                volumeNode, maskBase * islandPct / 100.0, int(islandErode))
            arr = np.array(slicer.util.arrayFromVolume(volumeNode), dtype=float)
            arr[~brainMask] = 0.0
        else:
            arr = slicer.util.arrayFromVolume(volumeNode)

        usePeak = (metric == REGION_METRIC_PEAK)
        if refMode == REGION_REF_MEAN:
            refSt = (self.regionRoiStats(referenceRoiNode, arr, r2i, floor)
                     if referenceRoiNode is not None else None)
            refval = refSt["mean"] if refSt else 0.0
            reflabel = "cerebellar mean"
        else:
            refval = cerebellarMax if cerebellarMax and cerebellarMax > 0 else robustReference
            reflabel = "cerebellar max"
        if not refval or refval <= 0:
            raise ValueError("Reference value unavailable (set the cerebellar max, "
                             "or select a cerebellum ROI for '%s')." % REGION_REF_MEAN)

        # "Hot" voxels are defined against the cerebellar MAX (the same anchor as
        # the active-scan red voxels), so the table agrees with the render.
        hotRef = cerebellarMax if cerebellarMax and cerebellarMax > 0 else refval
        hotThreshold = hotRef * hyperPct / 100.0

        refName = referenceRoiNode.GetName() if referenceRoiNode is not None else None
        rows = []
        for roi in slicer.util.getNodesByClass("vtkMRMLMarkupsROINode"):
            if referenceRoiNode is not None and roi is referenceRoiNode:
                continue
            nm = roi.GetName()
            if refName is not None and nm == refName:
                continue
            st = self.regionRoiStats(roi, arr, r2i, floor, hotThreshold=hotThreshold)
            if st is None:
                continue
            val = st["peak"] if usePeak else st["mean"]
            pct = val / refval * 100.0
            hotFracPct = st["hotFrac"] * 100.0
            hotMeanPct = (st["hotMean"] / refval * 100.0) if st["hotMean"] is not None else None
            focal = bool(focalHyper and (st["hotN"] >= minHotVoxels or hotFracPct >= minHotFrac))
            if focal or pct > hyperPct:
                cls = "Hyperactive"
            elif pct < hypoPct:
                cls = "Hypoactive"
            else:
                cls = "Normal"
            # Ordinal grade: hyper +1..+4 (hot-voxel mean: red/white/ceiling/top),
            # hypo -1..-4 (region mean stepping below the hypo threshold), normal 0.
            grade = 0
            if cls == "Hyperactive":
                hv = hotMeanPct if hotMeanPct is not None else pct
                grade = 1 + (hv >= hyperWhitePct) + (hv >= 100.0) + (hv >= hyperTopPct)
            elif cls == "Hypoactive":
                steps = int((hypoPct - pct) / hypoStepPct) if hypoStepPct > 0 else 0
                grade = -min(4, 1 + max(0, steps))
            schema = AMEN_REGIONS_BY_NAME.get(nm)
            rows.append(dict(
                name=nm,
                group=schema["group"] if schema else "(custom)",
                side=schema["side"] if schema else "?",
                base=schema["base"] if schema else nm,
                note=schema["note"] if schema else "",
                value=val, pct=pct, cls=cls, grade=grade, focal=focal,
                nvox=st["nvox"], hotN=st["hotN"], hotFrac=hotFracPct, hotMeanPct=hotMeanPct,
                asym=None, asymFlag=False))

        # L/R asymmetry for paired regions present on both sides (mean-based).
        byBaseSide = {}
        for row in rows:
            if row["side"] in ("L", "R"):
                byBaseSide.setdefault(row["base"], {})[row["side"]] = row
        for base, sides in byBaseSide.items():
            if "L" in sides and "R" in sides:
                lv, rv = sides["L"]["value"], sides["R"]["value"]
                denom = (lv + rv) / 2.0
                if denom > 0:
                    asym = (lv - rv) / denom * 100.0  # +ve = left higher
                    flag = abs(asym) > asymmetryPct
                    for s in ("L", "R"):
                        sides[s]["asym"] = asym
                        sides[s]["asymFlag"] = flag

        rows.sort(key=lambda r: r["pct"])  # most hypo first
        meta = dict(refval=refval, reflabel=reflabel, floor=floor, metric=metric,
                    hypoPct=hypoPct, hyperPct=hyperPct, asymmetryPct=asymmetryPct,
                    focalHyper=focalHyper, minHotVoxels=minHotVoxels, minHotFrac=minHotFrac,
                    hotThreshold=hotThreshold, masked=masked,
                    hyperWhitePct=hyperWhitePct, hyperTopPct=hyperTopPct, hypoStepPct=hypoStepPct,
                    nHypo=sum(1 for r in rows if r["cls"] == "Hypoactive"),
                    nHyper=sum(1 for r in rows if r["cls"] == "Hyperactive"),
                    gradeCounts={g: sum(1 for r in rows if r["grade"] == g)
                                 for g in range(-4, 5) if g != 0},
                    nFocal=sum(1 for r in rows if r["focal"]),
                    nAsym=sum(1 for r in rows if r["asymFlag"] and r["side"] == "L"),
                    nRegions=len(rows))
        return rows, meta

    def colorRegionRois(self, rows):
        """Tint each region's Box ROI by its ordinal grade."""
        for row in rows:
            roi = slicer.mrmlScene.GetFirstNodeByName(row["name"])
            if roi and roi.GetDisplayNode():
                c = rowDisplayColor(row)
                dn = roi.GetDisplayNode()
                dn.SetSelectedColor(*c)
                dn.SetColor(*c)

    def exportRegionTsv(self, rows, meta, outputDir):
        """Write the full table and a hypo/hyper-only TSV (matches the clinic's
        regional_hypo_hyper.tsv). Returns (fullPath, abnormalPath)."""
        import os
        os.makedirs(outputDir, exist_ok=True)
        metricCol = "Peak_counts" if meta["metric"] == REGION_METRIC_PEAK else "Mean_counts"
        header = ("Region\tGroup\tSide\t%s\tPct_of_%s\tHot_voxels\tHot_pct_of_ROI\tHot_mean_pct\t"
                  "Classification\tGrade\tFocal\tLR_asymmetry_pct\tVoxels\n" % (
                      metricCol, meta["reflabel"].replace(" ", "_")))

        def fmt(row):
            asym = "" if row["asym"] is None else "%.1f" % row["asym"]
            hotMean = "" if row["hotMeanPct"] is None else "%.1f" % row["hotMeanPct"]
            return "%s\t%s\t%s\t%.1f\t%.1f\t%d\t%.1f\t%s\t%s\t%s\t%s\t%s\t%d\n" % (
                row["name"], row["group"], row["side"], row["value"], row["pct"],
                row["hotN"], row["hotFrac"], hotMean, row["cls"], gradeLabel(row["grade"]),
                "yes" if row["focal"] else "", asym, row["nvox"])

        full = os.path.join(outputDir, "regional_all.tsv")
        with open(full, "w") as f:
            f.write(header)
            for row in rows:
                f.write(fmt(row))
        abn = os.path.join(outputDir, "regional_hypo_hyper.tsv")
        with open(abn, "w") as f:
            f.write(header)
            for row in rows:
                if row["cls"] in ("Hypoactive", "Hyperactive"):
                    f.write(fmt(row))
        return full, abn

    def formatRegionReport(self, rows, meta):
        valhdr = "PeakCts" if meta["metric"] == REGION_METRIC_PEAK else "MeanCts"
        lines = [
            "=== REGIONAL QUANTIFICATION (Amen-style) ===",
            "Metric = region %s   Reference (100%%) = %.1f counts (%s)   Volume = %s" % (
                "peak" if meta["metric"] == REGION_METRIC_PEAK else "mean",
                meta["refval"], meta["reflabel"],
                "brain-masked (matches render)" if meta.get("masked") else "RAW (incl. extracranial)"),
            "Hypoactive = MEAN < %.0f%%.   Hyperactive = MEAN > %.0f%% OR %sa focal hot "
            "focus: >= %d hot voxels OR >= %.0f%% of ROI, hot = >= %.0f%% (%.1f counts)." % (
                meta["hypoPct"], meta["hyperPct"],
                "" if meta["focalHyper"] else "[focal OFF] ",
                meta["minHotVoxels"], meta["minHotFrac"], meta["hyperPct"], meta["hotThreshold"]),
            "Ordinal grade (%% of cereb max, 100%% ceiling): hyper +1 red %.0f-%.0f%%, "
            "+2 white %.0f-100%%, +3 100-%.0f%%, +4 >= %.0f%%;  hypo -1..-4 step %.0f%% "
            "below %.0f%%." % (
                meta["hyperPct"], meta["hyperWhitePct"], meta["hyperWhitePct"],
                meta["hyperTopPct"], meta["hyperTopPct"], meta["hypoStepPct"], meta["hypoPct"]),
            "%-44s %8s  %6s  %11s   %-14s %5s %6s" % (
                "Region", valhdr, "%Ref", "Hot%(vox)", "Class", "Grade", "L-R%"),
            "-" * 104,
        ]
        for row in rows:
            asym = "" if row["asym"] is None else "%+.0f" % row["asym"]
            mark = " *" if row["asymFlag"] else ""
            clsTxt = row["cls"] + ("*" if row["focal"] else "")
            grd = gradeLabel(row["grade"]) if row["grade"] else ""
            hot = "%.1f%%(%d)" % (row["hotFrac"], row["hotN"])
            lines.append("%-44s %8.1f  %5.1f%%  %11s   %-14s %5s %6s%s" % (
                row["name"], row["value"], row["pct"], hot, clsTxt, grd, asym, mark))
        lines.append("-" * 104)
        gc = meta["gradeCounts"]
        hyperBreak = " / ".join("+%d:%d" % (g, gc[g]) for g in (1, 2, 3, 4) if gc[g])
        hypoBreak = " / ".join("%d:%d" % (g, gc[g]) for g in (-1, -2, -3, -4) if gc[g])
        lines.append("%d region(s): %d hypoactive [%s], %d hyperactive [%s] (%d focal*), "
                     "%d L/R asymmetr%s >%.0f%% (marked *)." % (
                         meta["nRegions"], meta["nHypo"], hypoBreak or "-",
                         meta["nHyper"], hyperBreak or "-", meta["nFocal"], meta["nAsym"],
                         "y" if meta["nAsym"] == 1 else "ies", meta["asymmetryPct"]))
        # Surface the clinically-loaded regions when abnormal.
        callouts = [r for r in rows if r["note"] and r["cls"] != "Normal"]
        for r in callouts:
            label = "%s (%s)" % (r["cls"].lower(), gradeLabel(r["grade"]))
            lines.append("  NOTE  %s is %s - %s." % (r["name"], label, r["note"]))
        return "\n".join(lines)

    # ---- view presets + montage -----------------------------------------
    def _renderVolumeForViews(self, volumeNode):
        # Both scan render volumes share the brain's geometry, so either gives the
        # same bounds; prefer a currently-visible one, else any, else the source.
        names = (SURFACE_RENDER_NAME, ACTIVE_RENDER_NAME, MASKED_VOLUME_NAME)
        nodes = [n for n in (slicer.mrmlScene.GetFirstNodeByName(nm) for nm in names)
                 if n and n.IsA("vtkMRMLScalarVolumeNode")]
        for n in nodes:
            vrLogic = slicer.modules.volumerendering.logic()
            dn = vrLogic.GetFirstVolumeRenderingDisplayNode(n)
            if dn and dn.GetVisibility():
                return n
        return nodes[0] if nodes else volumeNode

    def brainBounds(self, volumeNode):
        """RAS bounding box of the nonzero (brain) voxels."""
        import numpy as np
        import itertools
        a = slicer.util.arrayFromVolume(volumeNode)
        ks, js, iis = np.where(a > 0)
        if iis.size == 0:
            b = [0.0] * 6
            volumeNode.GetRASBounds(b)
            return b
        i2r = vtk.vtkMatrix4x4()
        volumeNode.GetIJKToRASMatrix(i2r)
        pts = np.array([i2r.MultiplyPoint([float(i), float(j), float(k), 1.0])[:3]
                        for i, j, k in itertools.product(
                            [int(iis.min()), int(iis.max())],
                            [int(js.min()), int(js.max())],
                            [int(ks.min()), int(ks.max())])])
        mn, mx = pts.min(0), pts.max(0)
        return [mn[0], mx[0], mn[1], mx[1], mn[2], mx[2]]

    def _aimCamera(self, ren, cam, bounds, direction, viewUp):
        C = [(bounds[0] + bounds[1]) / 2.0, (bounds[2] + bounds[3]) / 2.0,
             (bounds[4] + bounds[5]) / 2.0]
        diag = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
        cam.SetFocalPoint(*C)
        cam.SetPosition(C[0] + direction[0] * diag * 2.0,
                        C[1] + direction[1] * diag * 2.0,
                        C[2] + direction[2] * diag * 2.0)
        cam.SetViewUp(*viewUp)
        ren.ResetCamera(bounds)

    def snapToView(self, volumeNode, name):
        lm = slicer.app.layoutManager()
        if not lm or not lm.threeDWidget(0):
            return
        view = lm.threeDWidget(0).threeDView()
        ren = view.renderWindow().GetRenderers().GetFirstRenderer()
        cam = ren.GetActiveCamera()
        d, up = next(((dd, uu) for nm, dd, uu in VIEW_DIRECTIONS if nm == name),
                     ((0, 0, 1), (0, 1, 0)))
        self._aimCamera(ren, cam, self.brainBounds(self._renderVolumeForViews(volumeNode)), d, up)
        view.forceRender()

    def captureSixViews(self, volumeNode, outputDir, subdir, scanLabel=""):
        """Snap the 6 canonical views, save a PNG of each, and assemble a 2x3
        montage. Ported from the clinic's capture_views.py. Returns (montage, paths)."""
        import os
        lm = slicer.app.layoutManager()
        if not lm or not lm.threeDWidget(0):
            raise RuntimeError("No 3D view available to capture.")
        view = lm.threeDWidget(0).threeDView()
        ren = view.renderWindow().GetRenderers().GetFirstRenderer()
        cam = ren.GetActiveCamera()
        vn = view.mrmlViewNode()
        bounds = self.brainBounds(self._renderVolumeForViews(volumeNode))

        vn.SetBackgroundColor(0, 0, 0)
        vn.SetBackgroundColor2(0, 0, 0)
        try:
            vn.SetOrientationMarkerType(slicer.vtkMRMLAbstractViewNode.OrientationMarkerTypeAxes)
            vn.SetOrientationMarkerSize(slicer.vtkMRMLAbstractViewNode.OrientationMarkerSizeMedium)
        except Exception:
            pass
        ca = view.cornerAnnotation()
        ca.GetTextProperty().SetColor(1, 1, 1)
        try:
            ca.SetMaximumFontSize(28)
        except Exception:
            pass
        for ch in slicer.util.getNodesByClass("vtkMRMLCrosshairNode"):
            ch.SetCrosshairMode(slicer.vtkMRMLCrosshairNode.NoCrosshair)

        outdir = os.path.join(outputDir, subdir)
        os.makedirs(outdir, exist_ok=True)
        paths = []
        for name, d, up in VIEW_DIRECTIONS:
            self._aimCamera(ren, cam, bounds, d, up)
            ca.SetText(2, (scanLabel + " - " + name) if scanLabel else name)  # 2 = upper-left
            view.forceRender()
            w2i = vtk.vtkWindowToImageFilter()
            w2i.SetInput(view.renderWindow())
            w2i.ReadFrontBufferOff()
            w2i.Update()
            p = os.path.join(outdir, "%s.png" % name)
            wr = vtk.vtkPNGWriter()
            wr.SetFileName(p)
            wr.SetInputConnection(w2i.GetOutputPort())
            wr.Write()
            paths.append(p)
        ca.SetText(2, "")
        montage = os.path.join(outdir, "montage.png")
        self.assembleMontage(paths, montage)
        view.forceRender()
        return montage, paths

    def assembleMontage(self, imagePaths, montagePath):
        """2 rows x 3 cols: top row = first 3 views, bottom row = last 3."""
        imgs = []
        for p in imagePaths:
            r = vtk.vtkPNGReader()
            r.SetFileName(p)
            r.Update()
            o = vtk.vtkImageData()
            o.DeepCopy(r.GetOutput())
            imgs.append(o)

        def rowcat(three):
            ap = vtk.vtkImageAppend()
            ap.SetAppendAxis(0)
            for im in three:
                ap.AddInputData(im)
            ap.Update()
            out = vtk.vtkImageData()
            out.DeepCopy(ap.GetOutput())
            return out

        top, bot = rowcat(imgs[0:3]), rowcat(imgs[3:6])
        apv = vtk.vtkImageAppend()
        apv.SetAppendAxis(1)
        apv.AddInputData(bot)   # bottom row first (low Y), top row second
        apv.AddInputData(top)
        apv.Update()
        wr = vtk.vtkPNGWriter()
        wr.SetFileName(montagePath)
        wr.SetInputData(apv.GetOutput())
        wr.Write()
        return montagePath

    # ---- reporting -------------------------------------------------------
    def formatReport(self, volumeNode, stats, cerebellarMax, mode):
        dims = volumeNode.GetImageData().GetDimensions() if volumeNode.GetImageData() else (0, 0, 0)
        sp = volumeNode.GetSpacing()
        rng = volumeNode.GetImageData().GetScalarRange() if volumeNode.GetImageData() else (0, 0)
        ref = self.getReferenceValue(volumeNode, mode, cerebellarMax, stats)
        lines = [
            f"=== LOAD / CALIBRATION REPORT: {volumeNode.GetName()} ===",
            f"Dimensions (i,j,k): {dims}   Spacing (mm): ({sp[0]:.2f}, {sp[1]:.2f}, {sp[2]:.2f})",
            f"Scalar range: [{rng[0]:.1f}, {rng[1]:.1f}]   raw max: {stats['rawMax']:.1f}   "
            f"robust max: {stats['robustMax']:.1f}   mean: {stats['mean']:.1f}",
            f"Reference mode: {mode}",
            f"Cerebellar max set: {cerebellarMax:.1f}" if cerebellarMax > 0 else
            "Cerebellar max set: (not set - use Find cerebellar max or enter manually)",
            f"REFERENCE (100%) = {ref:.1f} counts",
            "Absolute count value behind each % of reference:",
            f"   45% (prognostic surface) = {0.45 * ref:.1f}",
            f"   55% (surface threshold)  = {0.55 * ref:.1f}",
            f"   60% (hypo/normal cut)    = {0.60 * ref:.1f}",
            f"   80% (active / hyper cut) = {0.80 * ref:.1f}",
            f"   90% (pink)               = {0.90 * ref:.1f}",
            f"  100% (white)              = {1.00 * ref:.1f}",
        ]
        return "\n".join(lines)


#
# SPECTBrainRenderTest
#
class SPECTBrainRenderTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_ReferenceAndStats()
        self.test_IslandMask()
        self.test_CerebellarCandidates()
        self.test_SurfaceVR()
        self.test_SurfaceMesh()
        self.test_ActiveScan()
        self.test_SeparateRenderVolumes()
        self.test_Montage()
        self.test_RegionSchema()
        self.test_RegionQuantify()
        self.test_RegionSelection()
        self.test_FocalHyperactive()
        self.test_QuantifyMasksExtracranial()
        self.test_OrdinalGrade()

    def _makePhantom(self):
        """Ellipsoid 'brain' (~600), a cerebellum-like hot region posterior-inferior
        (~870), a cortical hot spot (~760), a cold dent (~250), plus a single
        hot-pixel artifact (9000)."""
        import numpy as np
        nz, ny, nx = 48, 64, 64
        zz, yy, xx = np.mgrid[0:nz, 0:ny, 0:nx].astype(float)
        cz, cy, cx = nz / 2.0, ny / 2.0, nx / 2.0
        ell = ((xx - cx) / 26.0) ** 2 + ((yy - cy) / 26.0) ** 2 + ((zz - cz) / 18.0) ** 2
        arr = np.zeros((nz, ny, nx), dtype=float)
        arr[ell <= 1.0] = 600.0
        # cerebellum-like region: posterior (low y) and inferior (low z)
        cereb = ((xx - cx) ** 2 / 12.0 ** 2 + (yy - (cy - 16)) ** 2 / 8.0 ** 2 +
                 (zz - (cz - 12)) ** 2 / 6.0 ** 2) <= 1.0
        arr[cereb] = 870.0
        # cortical hot spot
        hot = ((xx - (cx + 12)) ** 2 + (yy - (cy + 10)) ** 2 + (zz - (cz + 4)) ** 2) < 6.0 ** 2
        arr[hot] = 760.0
        cold = ((xx - (cx - 14)) ** 2 + (yy - (cy - 2)) ** 2 + (zz - (cz + 2)) ** 2) < 5.0 ** 2
        arr[cold & (ell <= 1.0)] = 250.0
        arr += np.abs(np.sin(xx + yy + zz)) * 15.0
        arr[1, 1, 1] = 9000.0
        vol = slicer.util.addVolumeFromArray(arr.astype("float32"))
        vol.SetName("PhantomSPECT")
        vol.SetSpacing(3.0, 3.0, 3.0)
        return vol

    def test_ReferenceAndStats(self):
        self.delayDisplay("Reference + robust statistics")
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        stats, _ = logic.computeStatistics(vol)
        self.assertEqual(stats["rawMax"], 9000.0)
        self.assertLess(stats["robustMax"], 2000.0)  # ignores hot pixel
        # cerebellum reference uses the manual value
        self.assertEqual(logic.getReferenceValue(vol, REF_CEREBELLUM, 870.0, stats), 870.0)
        self.assertEqual(logic.getReferenceValue(vol, REF_RAW_MAX, 0.0, stats), 9000.0)
        print(f"TEST robustMax={stats['robustMax']:.1f} ref_cereb=870 ref_raw={stats['rawMax']:.0f}")
        self.delayDisplay("Reference/stats passed")

    def test_IslandMask(self):
        self.delayDisplay("Largest-island brain mask")
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        masked = logic.largestIslandVolume(vol, 0.30 * 870.0, erode=2)
        self.assertTrue(masked.IsA("vtkMRMLScalarVolumeNode"))
        import numpy as np
        ma = slicer.util.arrayFromVolume(masked)
        self.assertGreater(int((ma > 0).sum()), 1000)
        # hot-pixel artifact at (1,1,1) must be masked out (not connected to brain)
        self.assertEqual(float(ma[1, 1, 1]), 0.0)
        print(f"TEST island voxels={int((ma>0).sum())} artifact_masked={ma[1,1,1]==0}")
        self.delayDisplay("Island mask passed")

    def test_CerebellarCandidates(self):
        self.delayDisplay("Cerebellar max candidate finder")
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        roi = logic.createCerebellumRoi(vol)
        cands = logic.cerebellarMaxCandidates(vol, roi, topN=10)
        self.assertGreater(len(cands), 0)
        self.assertEqual(cands[0]["rank"], 1)
        # the ROI sits over the cerebellum region (~870), top candidate near it
        self.assertGreater(cands[0]["value"], 500.0)
        print(f"TEST cerebellar top value={cands[0]['value']:.1f} n={len(cands)}")
        self.delayDisplay("Cerebellar candidates passed")

    def test_SurfaceVR(self):
        self.delayDisplay("Surface volume render")
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        res = logic.generateSurfaceVR(vol, 55.0, 870.0, opacity=0.9)
        dn = res["displayNode"]
        self.assertIsNotNone(dn)
        vp = dn.GetVolumePropertyNode().GetVolumeProperty()
        # default scheme is hole-emphasis (6 points)
        self.assertEqual(vp.GetRGBTransferFunction(0).GetSize(), len(HOLE_COLOR_FRACTIONS))
        self.assertEqual(res["colorPoints"], len(HOLE_COLOR_FRACTIONS))
        # sharp opacity ramp has 4 control points and honors the peak opacity
        self.assertEqual(vp.GetScalarOpacity().GetSize(), 4)
        self.assertAlmostEqual(vp.GetScalarOpacity().GetValue(res["threshold"]), 0.9, delta=0.03)
        # rainbow scheme yields 7 points
        res2 = logic.generateSurfaceVR(vol, 55.0, 870.0, colorScheme=SURFACE_SCHEME_RAINBOW)
        self.assertEqual(res2["colorPoints"], len(AMEN_COLOR_FRACTIONS))
        print(f"TEST surfaceVR holePts={res['colorPoints']} rainbowPts={res2['colorPoints']} "
              f"peak@T={vp.GetScalarOpacity().GetValue(res['threshold']):.2f} "
              f"threshold={res['threshold']:.1f}")
        self.delayDisplay("Surface VR passed")

    def test_SurfaceMesh(self):
        self.delayDisplay("Surface isosurface mesh (marching cubes)")
        slicer.mrmlScene.Clear()
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        res = logic.generateSurfaceMesh(vol, 55.0, 870.0, opacity=0.95)
        model = slicer.mrmlScene.GetFirstNodeByName(SURFACE_MODEL_NAME)
        self.assertIsNotNone(model)
        self.assertTrue(model.IsA("vtkMRMLModelNode"))
        self.assertIsNotNone(model.GetPolyData())
        self.assertGreater(res["points"], 100)
        d = model.GetDisplayNode()
        self.assertTrue(d.GetVisibility())
        self.assertAlmostEqual(d.GetOpacity(), 0.95, delta=0.02)
        # point normals present -> smooth Gouraud shading, not faceted
        self.assertIsNotNone(model.GetPolyData().GetPointData().GetNormals())
        # solid scheme (logic default): single material, no scalar coloring
        self.assertFalse(res["colored"])
        self.assertFalse(d.GetScalarVisibility())
        # a lower threshold keeps more (above-threshold) tissue -> >= as many verts
        resLow = logic.generateSurfaceMesh(vol, 40.0, 870.0)
        self.assertGreaterEqual(resLow["points"], 100)
        # perfusion-tint scheme: per-vertex 3-component RGB array attached + scalars on
        res2 = logic.generateSurfaceMesh(vol, 55.0, 870.0, colorScheme=SURFACE_SCHEME_HOLE)
        self.assertTrue(res2["colored"])
        poly2 = slicer.mrmlScene.GetFirstNodeByName(SURFACE_MODEL_NAME).GetPolyData()
        rgb = poly2.GetPointData().GetArray("SurfaceRGB")
        self.assertIsNotNone(rgb)
        self.assertEqual(rgb.GetNumberOfComponents(), 3)
        self.assertEqual(rgb.GetNumberOfTuples(), poly2.GetNumberOfPoints())
        self.assertTrue(model.GetDisplayNode().GetScalarVisibility())
        # positional spectrum scheme: direct-RGB array that actually varies along
        # the S/I axis (a true rainbow, not one flat color)
        resSpec = logic.generateSurfaceMesh(vol, 55.0, 870.0, colorScheme=SURFACE_SCHEME_SPECTRUM)
        self.assertTrue(resSpec["colored"])
        polyS = slicer.mrmlScene.GetFirstNodeByName(SURFACE_MODEL_NAME).GetPolyData()
        rgbS = polyS.GetPointData().GetArray("SurfaceRGB")
        self.assertIsNotNone(rgbS)
        self.assertEqual(rgbS.GetNumberOfComponents(), 3)
        import numpy as _np
        from vtk.util import numpy_support as _ns
        cols = _ns.vtk_to_numpy(rgbS)
        self.assertGreater(int(_np.unique(cols, axis=0).shape[0]), 5)  # many distinct hues
        # mutual exclusion: VR hides the mesh; mesh hides the VR
        rvr = logic.generateSurfaceVR(vol, 55.0, 870.0)
        self.assertFalse(slicer.mrmlScene.GetFirstNodeByName(
            SURFACE_MODEL_NAME).GetDisplayNode().GetVisibility())
        self.assertTrue(rvr["displayNode"].GetVisibility())
        logic.generateSurfaceMesh(vol, 55.0, 870.0)
        self.assertFalse(rvr["displayNode"].GetVisibility())
        self.assertTrue(slicer.mrmlScene.GetFirstNodeByName(
            SURFACE_MODEL_NAME).GetDisplayNode().GetVisibility())
        print("TEST surfaceMesh points=%d opacity=%.2f solidColored=%s holeColored=%s "
              "normals=%s" % (res["points"], d.GetOpacity(), res["colored"], res2["colored"],
                              model.GetPolyData().GetPointData().GetNormals() is not None))
        self.delayDisplay("Surface mesh passed")

    def test_ActiveScan(self):
        self.delayDisplay("Active scan (wireframe + hot volume)")
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        res = logic.generateActive(vol, 870.0, activePct=85.0)
        self.assertGreater(res["meshPoints"], 50)
        blue = slicer.mrmlScene.GetFirstNodeByName(WIRE_MODEL_NAME)
        nodes = slicer.mrmlScene.GetFirstNodeByName(WIRE_NODES_NAME)
        self.assertIsNotNone(blue.GetPolyData())
        self.assertEqual(blue.GetDisplayNode().GetRepresentation(), 1)  # wireframe
        self.assertGreater(nodes.GetPolyData().GetNumberOfPoints(), 0)
        vp = res["displayNode"].GetVolumePropertyNode().GetVolumeProperty()
        ctf = vp.GetRGBTransferFunction(0)
        self.assertEqual(ctf.GetSize(), 4)  # red, red, white, white
        # red at the active threshold (85%), white at the white point (92%)
        red = [ctf.GetColor(870.0 * 0.85)[k] for k in range(3)]
        white = [ctf.GetColor(870.0 * 0.92)[k] for k in range(3)]
        self.assertAlmostEqual(red[0], 1.0, delta=0.05)
        self.assertLess(red[2], 0.2)                       # red has ~no blue
        self.assertGreater(min(white), 0.9)                # white at 92%
        # node radius 0 -> the sphere "nodes" are hidden (bare wireframe)
        res0 = logic.generateActive(vol, 870.0, activePct=85.0, nodeRadius=0.0)
        nodes0 = slicer.mrmlScene.GetFirstNodeByName(WIRE_NODES_NAME)
        self.assertFalse(nodes0.GetDisplayNode().GetVisibility())
        print(f"TEST active meshPoints={res['meshPoints']} active@={res['active']:.1f} "
              f"red@85={tuple(round(x,2) for x in red)} white@92={tuple(round(x,2) for x in white)}")
        self.delayDisplay("Active scan passed")

    def test_SeparateRenderVolumes(self):
        self.delayDisplay("Surface / active use separate render volumes")
        slicer.mrmlScene.Clear()
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        rs = logic.generateSurfaceVR(vol, 55.0, 870.0)
        ra = logic.generateActive(vol, 870.0, activePct=85.0)
        sVol, aVol = rs["renderVolume"], ra["renderVolume"]
        # distinct render volumes with the expected dedicated names
        self.assertIsNot(sVol, aVol)
        self.assertEqual(sVol.GetName(), SURFACE_RENDER_NAME)
        self.assertEqual(aVol.GetName(), ACTIVE_RENDER_NAME)
        # distinct VR display nodes - neither overwrote the other
        sdn, adn = rs["displayNode"], ra["displayNode"]
        self.assertIsNot(sdn, adn)
        # both render volumes still present after both renders (independent)
        self.assertIsNotNone(slicer.mrmlScene.GetFirstNodeByName(SURFACE_RENDER_NAME))
        self.assertIsNotNone(slicer.mrmlScene.GetFirstNodeByName(ACTIVE_RENDER_NAME))
        # surface keeps its color scheme; active keeps red->pink->white (4 pts)
        sPts = sdn.GetVolumePropertyNode().GetVolumeProperty().GetRGBTransferFunction(0).GetSize()
        aPts = adn.GetVolumePropertyNode().GetVolumeProperty().GetRGBTransferFunction(0).GetSize()
        self.assertEqual(aPts, 4)
        self.assertIn(sPts, (len(HOLE_COLOR_FRACTIONS), len(AMEN_COLOR_FRACTIONS)))
        # regenerating the surface must not delete the active volume (and vice versa)
        logic.generateSurfaceVR(vol, 55.0, 870.0)
        self.assertIsNotNone(slicer.mrmlScene.GetFirstNodeByName(ACTIVE_RENDER_NAME))
        print(f"TEST separate surfaceVol='{sVol.GetName()}' activeVol='{aVol.GetName()}' "
              f"distinctVR={sdn is not adn} surfPts={sPts} actPts={aPts}")
        self.delayDisplay("Separate render volumes passed")

    def _writeSolidPng(self, path, w, h, val):
        import numpy as np
        from vtk.util import numpy_support
        img = vtk.vtkImageData()
        img.SetDimensions(w, h, 1)
        arr = np.zeros((h * w, 3), dtype=np.uint8)
        arr[:] = val
        img.GetPointData().SetScalars(
            numpy_support.numpy_to_vtk(arr, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR))
        wr = vtk.vtkPNGWriter()
        wr.SetFileName(path)
        wr.SetInputData(img)
        wr.Write()

    def test_Montage(self):
        self.delayDisplay("6-view montage assembly")
        import os
        import tempfile
        logic = SPECTBrainRenderLogic()
        d = tempfile.mkdtemp()
        W, H = 30, 20
        paths = []
        for i, (name, _, _) in enumerate(VIEW_DIRECTIONS):
            p = os.path.join(d, name + ".png")
            self._writeSolidPng(p, W, H, (i * 40) % 256)
            paths.append(p)
        out = os.path.join(d, "montage.png")
        logic.assembleMontage(paths, out)
        self.assertTrue(os.path.exists(out))
        r = vtk.vtkPNGReader()
        r.SetFileName(out)
        r.Update()
        dims = r.GetOutput().GetDimensions()
        self.assertEqual(dims[0], W * 3)  # 3 columns
        self.assertEqual(dims[1], H * 2)  # 2 rows
        # brainBounds sanity on the phantom
        b = logic.brainBounds(self._makePhantom())
        self.assertGreater(b[1] - b[0], 0)
        print(f"TEST montage dims={dims} expected=({W * 3},{H * 2}) brainBoundsOK=True")
        self.delayDisplay("Montage assembly passed")

    def test_RegionSchema(self):
        self.delayDisplay("Region taxonomy schema")
        # Group counts implied by the enumerated Amen list.
        from collections import Counter
        groups = Counter(r["group"] for r in AMEN_REGIONS)
        self.assertEqual(groups["Frontal Lobe"], 14)    # 7 paired
        self.assertEqual(groups["Cingulate Gyrus"], 10)  # 5 paired
        self.assertEqual(groups["Temporal Lobe"], 12)    # 6 paired
        self.assertEqual(groups["Parietal Lobe"], 4)     # 2 paired
        self.assertEqual(groups["Occipital Lobe"], 4)    # 2 paired
        self.assertEqual(groups["Subcortical"], 8)       # 4 paired
        self.assertEqual(groups["Cerebellum"], 2)        # 1 paired
        self.assertEqual(len(AMEN_REGIONS), 54)  # 27 base regions x L/R
        # all regions paired L/R; names unique
        self.assertTrue(all(r["side"] in ("L", "R") for r in AMEN_REGIONS))
        self.assertEqual(len({r["base"] for r in AMEN_REGIONS}), 27)
        names = [r["name"] for r in AMEN_REGIONS]
        self.assertEqual(len(names), len(set(names)))
        self.assertIn("L Dorsolateral PFC", AMEN_REGIONS_BY_NAME)
        self.assertIn("R Dorsolateral PFC", AMEN_REGIONS_BY_NAME)
        self.assertIn("L Thalamus", AMEN_REGIONS_BY_NAME)
        self.assertIn("R Thalamus", AMEN_REGIONS_BY_NAME)
        # clinical note carried through to both sides of the paired region
        self.assertIn("SSRI", AMEN_REGIONS_BY_NAME["L Anterior Cingulate - Ventral"]["note"])
        self.assertIn("SSRI", AMEN_REGIONS_BY_NAME["R Anterior Cingulate - Ventral"]["note"])
        print(f"TEST schema regions={len(AMEN_REGIONS)} bases=27 groups={dict(groups)}")
        self.delayDisplay("Region schema passed")

    def test_RegionQuantify(self):
        self.delayDisplay("Regional quantification + scaffold")
        import os
        import tempfile
        slicer.mrmlScene.Clear()  # isolate: only the scaffold ROIs in the scene
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        n = logic.createRegionRoiScaffold(vol)
        self.assertEqual(n, 54)
        self.assertEqual(len(slicer.util.getNodesByClass("vtkMRMLMarkupsROINode")), 54)

        rows, meta = logic.quantifyRegions(vol, None, 870.0,
                                           metric=REGION_METRIC_MEAN, refMode=REGION_REF_MAX)
        self.assertGreater(len(rows), 0)
        # every row classified, pct = value / reference * 100
        valid = {"Hypoactive", "Normal", "Hyperactive"}
        self.assertTrue(all(r["cls"] in valid for r in rows))
        self.assertEqual(meta["nHypo"] + meta["nHyper"]
                         + sum(1 for r in rows if r["cls"] == "Normal"), len(rows))
        r0 = rows[0]
        self.assertAlmostEqual(r0["pct"], r0["value"] / 870.0 * 100.0, delta=0.01)
        # rows sorted most-hypo first
        self.assertLessEqual(rows[0]["pct"], rows[-1]["pct"])
        # at least one paired base present on both sides -> asymmetry computed correctly
        byname = {r["name"]: r for r in rows}
        checked = False
        for spec in REGION_SPECS:
            base = spec[1]
            L, R = byname.get("L " + base), byname.get("R " + base)
            if L and R:
                denom = (L["value"] + R["value"]) / 2.0
                self.assertAlmostEqual(L["asym"], (L["value"] - R["value"]) / denom * 100.0, delta=0.01)
                self.assertEqual(L["asym"], R["asym"])
                checked = True
                break
        self.assertTrue(checked, "expected at least one paired base on both sides")

        # '% of cerebellar MEAN' with no reference ROI must raise
        with self.assertRaises(ValueError):
            logic.quantifyRegions(vol, None, 870.0, refMode=REGION_REF_MEAN)

        # TSV export writes the full table and a hypo/hyper-only table
        d = tempfile.mkdtemp()
        full, abn = logic.exportRegionTsv(rows, meta, d)
        self.assertTrue(os.path.exists(full) and os.path.exists(abn))
        with open(full) as f:
            self.assertEqual(sum(1 for _ in f), len(rows) + 1)  # header + rows
        print(f"TEST quantify rows={len(rows)} hypo={meta['nHypo']} hyper={meta['nHyper']} "
              f"asymFlagged={meta['nAsym']} ref={meta['refval']:.1f}({meta['reflabel']})")
        self.delayDisplay("Regional quantification passed")

    def test_RegionSelection(self):
        self.delayDisplay("Selective region scaffold")
        logic = SPECTBrainRenderLogic()
        slicer.mrmlScene.Clear()
        vol = self._makePhantom()
        # create only a chosen subset
        subset = ["L Dorsolateral PFC", "R Dorsolateral PFC", "L Thalamus", "R Thalamus"]
        n = logic.createRegionRoiScaffold(vol, names=subset)
        self.assertEqual(n, len(subset))
        present = {r.GetName() for r in slicer.util.getNodesByClass("vtkMRMLMarkupsROINode")}
        self.assertEqual(present, set(subset))
        # quantify sees exactly the subset (no other ROIs in the scene)
        rows, _ = logic.quantifyRegions(vol, None, 870.0)
        self.assertEqual({r["name"] for r in rows}, set(subset))
        # adding more keeps the originals (create is additive)
        logic.createRegionRoiScaffold(vol, names=["L Insular Cortex"])
        self.assertEqual(len(slicer.util.getNodesByClass("vtkMRMLMarkupsROINode")), len(subset) + 1)
        # remove-unchecked reconciles to exactly the kept set, only touching taxonomy ROIs
        keep = ["L Thalamus", "R Thalamus"]
        removed = logic.removeUnselectedRegionRois(keep)
        self.assertEqual(removed, 3)  # 2 DLPFC + L Insular removed; L/R Thalamus kept
        present2 = {r.GetName() for r in slicer.util.getNodesByClass("vtkMRMLMarkupsROINode")}
        self.assertEqual(present2, set(keep))
        print(f"TEST selection created={n} after_remove={sorted(present2)}")
        self.delayDisplay("Region selection passed")

    def test_FocalHyperactive(self):
        self.delayDisplay("Focal hyperactivity (hot voxels vs diluted mean)")
        import numpy as np
        slicer.mrmlScene.Clear()
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()  # cortical hot spot ~760 (>85% of 870) at IJK (44,42,28)
        i2r = vtk.vtkMatrix4x4(); vol.GetIJKToRASMatrix(i2r)
        ras = i2r.MultiplyPoint([44.0, 42.0, 28.0, 1.0])[:3]
        roi = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", "L Posterior Frontal")
        roi.CreateDefaultDisplayNodes()
        roi.SetCenter(ras[0], ras[1], ras[2])
        roi.SetSize(45.0, 45.0, 45.0)  # loose box: hot core + surrounding ~600 tissue

        # focal ON: hyperactive from the hot cluster even though the mean is sub-threshold
        rows, meta = logic.quantifyRegions(vol, None, 870.0, focalHyper=True)
        r = next(x for x in rows if x["name"] == "L Posterior Frontal")
        self.assertEqual(r["cls"], "Hyperactive")
        self.assertTrue(r["focal"])
        self.assertGreaterEqual(r["hotN"], DEFAULT_MIN_HOT_VOXELS)
        self.assertLess(r["pct"], 85.0)               # the diluted mean alone would miss it
        self.assertIsNotNone(r["hotMeanPct"])
        self.assertGreaterEqual(r["hotMeanPct"], 85.0)  # hot-voxel mean is genuinely >= hyper %

        # COUNT criterion catches a small focus even when its FRACTION is tiny:
        # set the fraction trigger impossibly high so only the voxel count can fire.
        rows3, _ = logic.quantifyRegions(vol, None, 870.0, focalHyper=True,
                                         minHotVoxels=3, minHotFrac=100.0)
        r3 = next(x for x in rows3 if x["name"] == "L Posterior Frontal")
        self.assertEqual(r3["cls"], "Hyperactive")  # flagged via hotN, not fraction

        # focal OFF: classified by the mean -> Normal (reproduces the original complaint)
        rows2, _ = logic.quantifyRegions(vol, None, 870.0, focalHyper=False)
        r2 = next(x for x in rows2 if x["name"] == "L Posterior Frontal")
        self.assertEqual(r2["cls"], "Normal")
        self.assertFalse(r2["focal"])
        print(f"TEST focal mean%={r['pct']:.1f} hot%={r['hotFrac']:.1f} hotN={r['hotN']} "
              f"hotMean%={r['hotMeanPct']:.1f} focalON={r['cls']} countOnly={r3['cls']} focalOFF={r2['cls']}")
        self.delayDisplay("Focal hyperactivity passed")

    def test_QuantifyMasksExtracranial(self):
        self.delayDisplay("Brain mask excludes extracranial hot uptake")
        slicer.mrmlScene.Clear()
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        # add a DISCONNECTED extracranial hot blob (like salivary/sinus/scalp uptake,
        # never rendered red because the brain mask drops it) at a far corner
        arr = slicer.util.arrayFromVolume(vol)  # shape (k, j, i)
        arr[8:14, 8:14, 55:61] = 800.0  # ~216 voxels, 800 >= 85% of 870
        slicer.util.arrayFromVolumeModified(vol)
        i2r = vtk.vtkMatrix4x4(); vol.GetIJKToRASMatrix(i2r)
        ras = i2r.MultiplyPoint([58.0, 11.0, 11.0, 1.0])[:3]  # (i, j, k) center of the blob
        roi = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", "R Insular Cortex")
        roi.CreateDefaultDisplayNodes()
        roi.SetCenter(ras[0], ras[1], ras[2])
        roi.SetSize(30.0, 30.0, 30.0)

        # mask OFF: raw extracranial blob is counted -> false Hyperactive (the bug)
        rowsRaw, metaRaw = logic.quantifyRegions(vol, None, 870.0, maskToBrain=False)
        self.assertFalse(metaRaw["masked"])
        rRaw = next(x for x in rowsRaw if x["name"] == "R Insular Cortex")
        self.assertEqual(rRaw["cls"], "Hyperactive")
        self.assertGreater(rRaw["hotN"], DEFAULT_MIN_HOT_VOXELS)

        # mask ON: blob is not the largest island -> dropped; the all-extracranial
        # ROI then has no brain tissue and is excluded (no false hyperactive)
        rowsMask, meta = logic.quantifyRegions(vol, None, 870.0, maskToBrain=True,
                                               islandPct=30.0, islandErode=2)
        self.assertTrue(meta["masked"])
        self.assertNotIn("R Insular Cortex", {x["name"] for x in rowsMask})
        # the source volume must be untouched by masking (we worked on a copy)
        self.assertEqual(float(slicer.util.arrayFromVolume(vol)[11, 11, 58]), 800.0)
        print(f"TEST extracranial rawHotN={rRaw['hotN']} rawClass={rRaw['cls']} "
              f"maskedDropped={'R Insular Cortex' not in {x['name'] for x in rowsMask}}")
        self.delayDisplay("Extracranial masking passed")

    def test_OrdinalGrade(self):
        self.delayDisplay("Ordinal grade (-4..-1 hypo / +1..+4 hyper)")
        slicer.mrmlScene.Clear()
        logic = SPECTBrainRenderLogic()
        vol = self._makePhantom()
        arr = slicer.util.arrayFromVolume(vol)  # (k, j, i)
        arr[:] = 0.0  # blank canvas: each blob is isolated, surrounded by background
        # uniform blobs at known % of 870 -> known grades. (k, j, i) cubes.
        blobs = {  # name: (value, (k0, j0, i0))
            "+1 red": (770.0, (8, 8, 8)),       # 88.5% -> +1
            "+2 white": (835.0, (8, 8, 20)),    # 96.0% -> +2
            "+3 ceiling": (915.0, (8, 8, 32)),  # 105.2% -> +3
            "+4 top": (1010.0, (8, 8, 44)),     # 116.1% -> +4
            "-1 mild": (435.0, (8, 20, 8)),     # 50.0% -> -1
            "-4 severe": (174.0, (8, 20, 20)),  # 20.0% -> -4
            "normal": (600.0, (8, 20, 32)),     # 69.0% -> 0
        }
        i2r = vtk.vtkMatrix4x4(); vol.GetIJKToRASMatrix(i2r)
        for name, (value, (k0, j0, i0)) in blobs.items():
            arr[k0:k0 + 4, j0:j0 + 4, i0:i0 + 4] = value
            ras = i2r.MultiplyPoint([float(i0 + 1.5), float(j0 + 1.5), float(k0 + 1.5), 1.0])[:3]
            roi = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", name)
            roi.CreateDefaultDisplayNodes()
            roi.SetCenter(ras[0], ras[1], ras[2])
            roi.SetSize(12.0, 12.0, 12.0)
        slicer.util.arrayFromVolumeModified(vol)
        rows, meta = logic.quantifyRegions(vol, None, 870.0, maskToBrain=False)
        by = {r["name"]: r for r in rows}
        self.assertEqual(by["+1 red"]["grade"], 1)
        self.assertEqual(by["+2 white"]["grade"], 2)
        self.assertEqual(by["+3 ceiling"]["grade"], 3)
        self.assertEqual(by["+4 top"]["grade"], 4)
        self.assertEqual(by["-1 mild"]["grade"], -1)
        self.assertEqual(by["-4 severe"]["grade"], -4)
        self.assertEqual(by["normal"]["grade"], 0)
        self.assertEqual(by["normal"]["cls"], "Normal")
        gc = meta["gradeCounts"]
        self.assertEqual([gc[g] for g in (-4, -1, 1, 2, 3, 4)], [1, 1, 1, 1, 1, 1])
        self.assertEqual(gradeLabel(2), "+2")
        self.assertEqual(gradeLabel(-3), "-3")
        self.assertEqual(gradeLabel(0), "0")
        print("TEST grade +1=%d +2=%d +3=%d +4=%d -1=%d -4=%d normal=%d" % (
            by["+1 red"]["grade"], by["+2 white"]["grade"], by["+3 ceiling"]["grade"],
            by["+4 top"]["grade"], by["-1 mild"]["grade"], by["-4 severe"]["grade"],
            by["normal"]["grade"]))
        self.delayDisplay("Ordinal grade passed")
