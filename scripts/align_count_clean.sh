#!/bin/bash
cd ~/rnaseq-meme-kanseri

GTF="ref/Homo_sapiens.GRCh38.110.gtf"
INDEX="ref/GRCh38_index"

while read srr; do
    if [ -f "results/${srr}_counts.txt" ]; then
        echo "$srr zaten islenmis, atlaniyor."
        continue
    fi

    # Eger FASTQ yoksa indir
    if [ ! -f "raw_data/${srr}_1.fastq" ]; then
        echo "$srr indiriliyor..."
        prefetch $srr -O raw_data/
        fasterq-dump raw_data/$srr -O raw_data/ -e 4
        rm -rf raw_data/$srr   # .sra onbellek klasorunu sil
    fi

    echo "$srr hizalaniyor..."
    hisat2 -p 4 -x $INDEX \
        -1 raw_data/${srr}_1.fastq -2 raw_data/${srr}_2.fastq \
        --summary-file results/${srr}_align_summary.txt \
        | samtools sort -@ 4 -o results/${srr}.bam -

    echo "$srr icin sayim yapiliyor..."
    featureCounts -T 4 -p -a $GTF -o results/${srr}_counts.txt results/${srr}.bam

    echo "$srr ham dosyalari siliniyor (yer acmak icin)..."
    rm -f raw_data/${srr}_1.fastq raw_data/${srr}_2.fastq
    rm -f results/${srr}.bam

    echo "$srr TAMAMLANDI."
done < raw_data/srr_list.txt

echo "TUM ORNEKLER ISLENDI."
