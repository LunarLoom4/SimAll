# Overlay triplet -- Release-only build of all vcpkg dependencies.
#
# Why this exists
# ---------------
# The default x64-windows triplet builds every port BOTH Debug AND Release.
# For our release pipeline we never ship Debug binaries, so the Debug
# variants of Qt6 + VTK + OpenCASCADE + HDF5 are pure waste -- they roughly
# double the dependency build time, pushing the first uncached release run
# past GitHub's hard 6-hour job timeout.
#
# Setting VCPKG_BUILD_TYPE=release tells every port to skip the Debug
# configuration entirely, halving the first-time build.
#
# Everything else mirrors the standard x64-windows triplet.

set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE dynamic)
set(VCPKG_BUILD_TYPE release)
