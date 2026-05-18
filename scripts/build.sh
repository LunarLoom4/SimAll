#!/usr/bin/env bash
# SimAll Beta — build driver (Linux / macOS).
# Usage: scripts/build.sh [preset] [target]
#   preset : cmake configure preset    (default: linux-release)
#   target : optional build target

set -euo pipefail

PRESET="${1:-linux-release}"
TARGET="${2:-}"
JOBS="${JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu)}"

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"

echo "==> Configuring preset '$PRESET'"
cmake --preset "$PRESET"

BUILD_ARGS=(--build build -j "$JOBS")
if [[ -n "$TARGET" ]]; then
    BUILD_ARGS+=(--target "$TARGET")
fi

echo "==> cmake ${BUILD_ARGS[*]}"
cmake "${BUILD_ARGS[@]}"
echo "==> OK"
