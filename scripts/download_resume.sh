#!/bin/bash
cd ~/rnaseq-meme-kanseri
while read srr; do
    if [ -f "raw_data/${srr}_1.fastq" ] || [ -f "raw_data/${srr}.fastq" ]; then
        echo "$srr zaten tamamlanmis, atlaniyor."
    else
        echo "$srr indiriliyor..."
        prefetch $srr -O raw_data/
        fasterq-dump raw_data/$srr -O raw_data/ -e 4
    fi
done < raw_data/srr_list.txt
