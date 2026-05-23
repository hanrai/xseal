#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Deprecated name: forwards to genome-partition baseline (large per-thread slice).
# For 128 KiB block dispatch use repro_block128k_scaling.sh instead.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
echo "NOTE: repro_thread_scaling.sh -> repro_genome_partition_scaling.sh (legacy genome partition)" >&2
echo "      For 128 KiB blocks use: repro_block128k_scaling.sh" >&2
exec "$ROOT/experiments/rigorous_bench/scripts/repro_genome_partition_scaling.sh" "$@"
