cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(AmenStyleSPECT)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/your-org/AmenStyleSPECT")
set(EXTENSION_CATEGORY "Nuclear Medicine")
set(EXTENSION_CONTRIBUTORS "Dominic Velasco, MD (Brain SPECT Clinic)")
set(EXTENSION_DESCRIPTION "Amen-style perfusion SPECT visualization: generates a SURFACE rendering (hypoperfusion shown as holes/dents) and an ACTIVE rendering (hyperperfusion shown as red/white hot zones) from a reconstructed tomographic SPECT volume, using explicit, configurable percent-of-max thresholds. Research / visualization aid only; not a diagnostic device.")
set(EXTENSION_ICONURL "")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Build (only required when packaging through the Extensions Index / CMake).
# Developer-Mode install does NOT need this file to be built.
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SPECTBrainRender)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
