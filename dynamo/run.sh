#!/usr/bin/env bash
# run.sh [RUNDIR] -- one Surya run from init.dat in a fresh directory, then the butterfly plot.
#
# surya_2.f opens rad.dat, butbot.dat, series_b, ber.dat and psi.dat with access='append', so a second run in the
# same directory stacks its output onto the first one without any warning. Hence the empty-directory check.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
dir=${1:-run}

make -s -C "$here" surya
if [ -d "$dir" ] && [ -n "$(ls -A "$dir")" ]; then
    echo "$dir is not empty -- remove it or pass another directory" >&2
    exit 1
fi
mkdir -p "$dir"
cp "$here/init.dat" "$dir/"
(cd "$dir" && time "$here/surya")
python3 "$here/plot_run.py" "$dir" --out "$dir"
