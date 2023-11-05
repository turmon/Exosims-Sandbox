#!/bin/sh -f

SNRs="10 15 20 30 40"
Rss="35 105 140 200"

for snr in $SNRs; do 
  for rs in $Rss; do 
    echo SNR=$snr and Rs=$rs
    sed -e "s/VAL_SNR/$snr/" -e "s/VAL_RS/$rs/" < test/foo2.json > test/H6C_CODulzE_baseA_SNR_${snr}_Rs_${rs}.json
  done
done
