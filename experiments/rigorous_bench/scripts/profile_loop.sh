#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Loop a command N times (for long enough AMD uProf / perf windows).
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "usage: $0 <count> <cmd> [args...]" >&2
  exit 2
fi
n="$1"
shift
for ((i = 1; i <= n; i++)); do
  "$@"
done
