#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/_env.sh"
R1=$(mktemp -d)
R2=$(mktemp -d)
trap 'rm -rf "$R1" "$R2"' EXIT
for D in "$R1" "$R2"; do
  ./galaxy_sim galaxy --output_dir="$D" --physics_package=Newtonian \
    --n_stars=40 --n_steps=2 --snapshot_every=1 --save_run_info=true \
    --galaxy_init_seed=13579 --galaxy_init_template=symmetric_disk
done
cmp -s "$R1/snapshot_00000.csv" "$R2/snapshot_00000.csv"
