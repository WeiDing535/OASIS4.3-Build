#!/bin/bash
# Compatibility entry point: PhaseBuild.command has been merged into OASIS4.3.command.
# This file is kept only for backward compatibility with existing user habits; the actual launch logic is unified in OASIS4.3.command.
HERE="$(cd -- "$(dirname -- "$0")" && pwd)"
exec "$HERE/OASIS4.3.command" "$@"
