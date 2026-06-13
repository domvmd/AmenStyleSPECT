# ============================================================
# AmenStyleSPECT - one-paste install helper for a new PC
# ------------------------------------------------------------
# Registers the SPECTBrainRender module path with Slicer and
# persists it, so the module loads automatically on restart.
#
# HOW TO USE (on the target PC):
#   1. Copy the AmenStyleSPECT folder onto this PC.
#   2. Edit MODULE_DIR below to point at the SPECTBrainRender
#      folder on THIS PC (the folder that contains
#      SPECTBrainRender.py).
#   3. In Slicer: View > Python Console, paste this whole file,
#      press Enter, and restart when prompted.
#
# If the module does not appear after restart, use the GUI method:
#   Edit > Application Settings > Modules > Additional module paths > Add.
# ============================================================

import os
import qt
import slicer

# ---- EDIT THIS to the SPECTBrainRender folder on this PC ----
# Windows example:  r"C:\Users\you\Documents\Slicer\AmenStyleSPECT\SPECTBrainRender"
# macOS/Linux:      "/Users/you/Slicer/AmenStyleSPECT/SPECTBrainRender"
MODULE_DIR = r"CHANGE_ME/AmenStyleSPECT/SPECTBrainRender"
# -------------------------------------------------------------


def install_amenstyle_spect(module_dir):
    # If MODULE_DIR wasn't edited (or doesn't exist), pop a folder picker so this
    # works as a true one-paste with no editing.
    if ("CHANGE_ME" in module_dir) or not os.path.isdir(os.path.expanduser(module_dir)):
        picked = qt.QFileDialog.getExistingDirectory(
            None, "Select the SPECTBrainRender folder (contains SPECTBrainRender.py)")
        if not picked:
            print("Install cancelled - no folder selected.")
            return
        module_dir = picked

    module_dir = os.path.abspath(os.path.expanduser(module_dir))
    script = os.path.join(module_dir, "SPECTBrainRender.py")
    if not os.path.isfile(script):
        raise RuntimeError(
            "SPECTBrainRender.py not found in:\n  %s\n"
            "Pick the folder that directly contains SPECTBrainRender.py." % module_dir)

    # Read the current additional module paths, append ours, de-duplicate.
    settings = slicer.app.revisionUserSettings()
    raw = settings.value("Modules/AdditionalPaths")
    if raw is None:
        paths = []
    elif isinstance(raw, str):
        paths = [raw]
    else:
        paths = list(raw)
    norm = os.path.normpath(module_dir)
    if norm not in [os.path.normpath(p) for p in paths]:
        paths.append(module_dir)
        settings.setValue("Modules/AdditionalPaths", paths)
        added = True
    else:
        added = True  # already present is also success

    print("AmenStyleSPECT path registered:\n  %s" % module_dir)
    print("Current additional module paths:")
    for p in paths:
        print("   - %s" % p)

    if added:
        if slicer.util.confirmYesNoDisplay(
                "AmenStyleSPECT registered.\nRestart Slicer now to load it?"):
            slicer.util.restart()
        else:
            print("Restart Slicer when ready; the module loads on the next start.")


install_amenstyle_spect(MODULE_DIR)
