#!/bin/bash

set -euo pipefail

DATASET_DIR="$1"

if [ ! -d "$DATASET_DIR" ]; then
	mkdir "$DATASET_DIR"
fi

pushd "$DATASET_DIR"
echo "Downloading to $DATASET_DIR"
wget -O cmps_elegans_multi.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/Ef359x9kGYBGj2u9tqNT5I0BRiIUqiGTEyUfa8fAeZxr7w?e=vZLaly&download=1'
wget -O seqs_elegans_seqs.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/Ec5vuzv1PRdHimtQOSmnNHwBvQH_UYkvwpQ8s_anqkFidQ?e=lbuZeN&download=1'
wget -O cmps_ecoli_multi_single.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/Eb7k4SH0_3pAgWTdRch5PCcBP5qqorg-npgcfKzmLozyiw?e=3oUDO8&download=1'
wget -O seqs_ecoli_multi_single.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/EUP1_e2wRaNJvc8FAorMNsEBKJdVXZhG0q4gSJXghm854A?e=xEYtB9&download=1'
wget -O cmps_ecoli100_multi.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/EVtNR7Yaq_dMvyQH5JSIYZsBCerh4lM3rzjCjhnKtVXB_w?e=380yq9&download=1'
wget -O seqs_ecoli100_multi.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/EYARatFXtaRCmWVyvvG6OtoBC3lTWAh1rlG8trMIjM_h2g?e=6FSb88&download=1'
wget -O cmps_simulated85_multi.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/EQjeZHwGfblLhBx-mxjzKSMBwqH_CuHhMXNvn4GDAnWM8g?e=d53abr&download=1'
wget -O seqs_simulated85_multi.json.gz 'https://charitede-my.sharepoint.com/:u:/g/personal/max_zhao_charite_de/Ef0gzSn7BNZFsLiRdf3T2NkBj3UQzM-0ySLiHTlnNF4UbQ?e=obbv2l&download=1'

echo "Decompressing all downloaded files"
if command -v pigz &> /dev/null; then
	pigz -d *.json.gz
else
	gunzip *.json.gz
fi