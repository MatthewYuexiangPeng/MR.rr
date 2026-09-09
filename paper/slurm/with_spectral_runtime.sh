#!/usr/bin/env bash
# Source before R starts. A fixed numerical environment is required on every job.
export MRRR_R_MODULE="${MRRR_R_MODULE:-r/4.5.1-5zezbzn}"
if command -v module >/dev/null 2>&1; then module load "$MRRR_R_MODULE"; fi
export R_LIBS="${R_LIBS:-$HOME/R/mrrr-4.5:$HOME/R/library}"
export OPENBLAS_CORETYPE=HASWELL OPENBLAS_VERBOSE=2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [[ "$(uname -s)" == Linux ]]; then
  # HASWELL kernels require AVX2 and FMA. Reject an incompatible allocated node.
  if ! awk '/^flags[[:space:]]*:/ {
    avx2=0; fma=0; for (i=1; i<=NF; i++) {if ($i=="avx2") avx2=1; if ($i=="fma") fma=1}
    if (!avx2 || !fma) bad=1; seen=1
  } END {exit (!seen || bad)}' /proc/cpuinfo; then
    echo 'Allocated CPU lacks AVX2/FMA for HASWELL. Select compatible nodes using MRRR_REMAIN_CONSTRAINT or MRRR_REMAIN_NODE.' >&2
    return 2 2>/dev/null || exit 2
  fi
fi
command -v Rscript >/dev/null || { echo 'Rscript unavailable after module loading.' >&2; return 2 2>/dev/null || exit 2; }
