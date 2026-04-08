# reads mapping
## unzip the reference files
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.57.gtf.gz

## Build the STAR genome index
cd Analysis_replication
mkdir -p star_index

## mapping
ulimit -n 100000
mkdir -p Mapping
index="/Analysis_replication/star_index"
raw="raw_data"
out="Mapping"
threads=8          
max_jobs=8         
total=$(ls ${raw}/*_1.fastq.gz | wc -l)
i=0
run_star() {
    r1=$1
    base=$(basename $r1 _1.fastq.gz)
    r2=${raw}/${base}_2.fastq.gz
    echo "Starting $base"

    STAR \
    --runThreadN $threads \
    --genomeDir $index \
    --readFilesIn $r1 $r2 \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outTmpDir ${out}/${base}_STARtmp \
    --outFileNamePrefix ${out}/${base}_

    echo "Finished $base"
}

export -f run_star
export raw out index threads
for r1 in ${raw}/*_1.fastq.gz
do
    base=$(basename $r1 _1.fastq.gz)
    ((i++))
    echo "[$i/$total] Queueing $base"
    if [ -f ${out}/${base}_Aligned.sortedByCoord.out.bam ]; then
        echo "Skipping $base (already finished)"
        continue
    fi
    run_star $r1 &
    while [ $(jobs -r | wc -l) -ge $max_jobs ]
    do
        sleep 5
    done
done
wait
echo "All samples finished."

# Step 2: make mapping table
cd /stor/work/Chen/current_member/qiaoyunai/Analysis_replication/Mapping
echo -e "Sample Total_reads(M) Mapped_reads(M) Mapped_reads(%) Uniquely_mapped(%) Multi_mapped(%) Unmapped_too_short(%)" > mapping_table.txt

for f in *_Log.final.out; do
    sample=$(basename "$f" _Log.final.out)

    total=$(grep "Number of input reads" "$f" | awk -F'|' '{print $2}' | tr -d ' ')
    unique=$(grep "Uniquely mapped reads %" "$f" | awk -F'|' '{print $2}' | tr -d ' %')
    multi=$(grep "% of reads mapped to multiple loci" "$f" | awk -F'|' '{print $2}' | tr -d ' %')
    short=$(grep "% of reads unmapped: too short" "$f" | awk -F'|' '{print $2}' | tr -d ' %')

    mapped_pct=$(echo "$unique + $multi" | bc)
    mapped_reads=$(echo "$total * $mapped_pct / 100" | bc)

    total_m=$(echo "scale=2; $total/1000000" | bc)
    mapped_m=$(echo "scale=2; $mapped_reads/1000000" | bc)

    printf "%-12s %-14s %-15s %-15s %-18s %-15s %-20s\n" \
    "$sample" "$total_m" "$mapped_m" "$mapped_pct" "$unique" "$multi" "$short" >> mapping_table.txt

done

# Make gene_count_matrix
## extract gene IDs
cut -f1 ERR10359863_ReadsPerGene.out.tab | tail -n +5 > genes.txt

## extract counts
for f in *_ReadsPerGene.out.tab
do
    cut -f2 $f | tail -n +5 > ${f%.out.tab}.counts
done

## merge
paste genes.txt *.counts > gene_count_matrix.txt