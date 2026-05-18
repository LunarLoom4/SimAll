#!/usr/bin/env bash
# =============================================================================
# SimAll Beta - scripts/package.sh
# Configure + build + cpack the project for the host platform.
# Produces installers under build/_packages/.
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

PRESET="${1:-linux-release}"
BUILD_DIR="${ROOT}/build"
OUT_DIR="${BUILD_DIR}/_packages"

echo "[package.sh] configure preset=${PRESET}"
cmake --preset "${PRESET}"

echo "[package.sh] build"
cmake --build "${BUILD_DIR}" --parallel

echo "[package.sh] cpack"
mkdir -p "${OUT_DIR}"
( cd "${BUILD_DIR}" && cpack -B "${OUT_DIR}" --config CPackConfig.cmake )

echo "[package.sh] artefacts:"
ls -lh "${OUT_DIR}"
