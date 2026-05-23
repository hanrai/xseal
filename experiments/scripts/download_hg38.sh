#!/usr/bin/env bash

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 hanrai. All Rights Reserved.

# Download GRCh38/hg38 reference FASTA from UCSC into repo data/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${HG38_OUT:-$ROOT/data/hg38.fa}"
URL="${HG38_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz}"
GZ="${OUT}.gz"

mkdir -p "$(dirname "$OUT")"

if [[ -f "$OUT" ]]; then
  echo "hg38 already present: $OUT ($(wc -c <"$OUT") bytes)"
  exit 0
fi

echo "Downloading $URL ..."
wget -c -O "$GZ" "$URL"
echo "Decompressing ..."
gunzip -f "$GZ"
echo "Done: $OUT ($(wc -c <"$OUT") bytes)"
md5sum "$OUT" || true
