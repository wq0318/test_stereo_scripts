#!/bin/bash
set +e

SCRIPT_DIR="/test_stereo_scripts/analysis"
export SCRIPT_DIR

exec "${SCRIPT_DIR}/downstream_analysis_impl.sh" "$@"
