#!/bin/bash
# Source from the holstein-diagmc root on Kestrel.
# Use the installed executables directly: the RHEL 9 module hierarchy does not
# expose these RHEL 8 module names. Spack embeds their library search paths.
export HOLSTEIN_CXX="/nopt/nrel/apps/cpu_stack/compilers/06-24/linux-rhel8-sapphirerapids/gcc-12.3.0/gcc-12.3.0-3bukajtjn5ntkb5nm762hdjtnh36hnh4/bin/g++"
export PATH="$(dirname "$HOLSTEIN_CXX"):$PATH"
export HOLSTEIN_PYTHON="${HOLSTEIN_PYTHON:-$PWD/../venv/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONUNBUFFERED=1
test -x "$HOLSTEIN_PYTHON"
test -x "$HOLSTEIN_CXX"
