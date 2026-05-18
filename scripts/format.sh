#!/usr/bin/env bash
# =============================================================================
# SimAll Beta - scripts/format.sh
# Apply project clang-format style to every C/C++ source file under src/,
# tests/, apps/, plugins/. Idempotent; safe to run from CI or pre-commit.
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if ! command -v clang-format >/dev/null 2>&1 ; then
    echo "[format.sh] ERROR: clang-format not on PATH" >&2
    exit 1
fi

CF_VER="$(clang-format --version | head -n1)"
echo "[format.sh] Using $CF_VER"

mapfile -t FILES < <(find src tests apps plugins \
        -type f \( -name '*.hpp' -o -name '*.cpp' -o -name '*.h' -o -name '*.c' \) \
        2>/dev/null | sort)

if [ "${#FILES[@]}" -eq 0 ]; then
    echo "[format.sh] no sources found"; exit 0
fi

echo "[format.sh] formatting ${#FILES[@]} files…"
# Process in batches of 64 to keep argv short on macOS / WSL.
printf '%s\n' "${FILES[@]}" | xargs -n 64 clang-format -i --style=file
echo "[format.sh] done."
