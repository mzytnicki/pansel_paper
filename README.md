# Analysis for the `pansel` paper

## Architecture

    - Data
      - Graphs
      - Annotations
      - Genomes
      - Alignments
    - Results
      - 1000
      - 10000
      - 100000

## Download data

### Graphs

    cd Data/Graphs
    for i in `seq 2 22`
    do
      wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2022_03_11_minigraph_cactus/chrom-graphs-hprc-v1.1-mc-chm13-full/chr${i}.vg
      vg convert -fW chr${i}.vg | gzip -c > chr_${i}_hprc-v1.1-mc-chm13-full.gfa.gz
      rm chr${i}.vg
    done

### Annotations

HG38 gaps

    wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.agp.gz | gzip -dc | awk '($5 == "N") || ($5 == "U")' | cut -f 1-3 > Data/Annotations/Annotationshg38.gaps.bed

Structural variations

    wget -O - https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz | gzip -dc | sed '/^#/! s/^/chr/g' | grep -v "" > Data/Annotations/GRCh38.variant_call.all.vcf

Conservation, from PhyloP

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

Coding exons

    wget -O - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz | gzip -dc | awk '$3 == "exon"' | grep 'gene_type "protein_coding"' > Data/Annotations/gencode.v44.annotation_coding_exons.gtf

ChromHMM results

    wget https://public.hoffman2.idre.ucla.edu/ernst/2K9RS//full_stack/full_stack_annotation_public_release/hg38/hg38_genome_100_segments.bed.gz

Human alignments in HAL format

    wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2022_03_11_minigraph_cactus/hprc-v1.1-mc-grch38-full.hal

List of separate genomes

    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/y1v2_all/hprc.y1v2.fa.gz

Conserved regions (PSTs)

    wget --no-check-certificate https://dna-discovery.stanford.edu/publicmaterial/datasets/pangenome/pst-31mer.grch38.bed.gz


## Analysis

### Compute pansel results for different resolutions

    for resolution in 1000 10000 100000
    do
        for i in `seq 1 22`
        do
          sed "s/^/chr${i}\t/g" Data/Graphs/chr_${i}_hprc-v1.1-mc-chm13-full_GRCh38.0.chr${i}.tsv
        done > Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.tsv
    done

### Exclude gaps, and divide into different 8 bins (most conserved, to least conserved)

    for resolution in 1000 10000 100000
    do
        python3 divideRegions.py Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.tsv ${resolution} Data/Annotations/hg38.gaps.bed 8 Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.regions
    done

### Transform `pansel` results to BED format

    for binSize in 1000 10000 100000
    do
        for i in `seq 1 22`
        do
            sed "s/^/chr${i}\t/g" Results/${binSize}/chr_${i}_hprc-v1.1-mc-chm13-full_GRCh38.0.chr${i}.tsv
        done | awk '{print $1 "\t" ($3-1) "\t" $4 "\tbin_" NR "\t" ($5 / ($4 - $3)) "\t+"}' > Results/${binSize}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.bed
    done

### Fit to conserved/divergent distributions

    for binSize in 1000 10000 100000
    do
        Rscript getExtremes.R -i Results/${binSize}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.bed -p 0.05 -P 0.05 -t Results/${binSize}/fit.png -o Results/${binSize}/fitConserved.bed -O /Scratch/mazytnicki/PanSel/Results/${binSize}/fitDivergent.bed &> Results/${binSize}/fit.log
    done

### Compute number of overlaps with structural variations

    for resolution in 1000 10000 100000
    do
        out=Results/${resolution}/variants.tsv
        echo -e "category\tscore\ttype" > $out
        for i in $( seq 0 7 )
        do 
            echo -e -n "$i\t"
            bedtools intersect -a Data/Annotations/GRCh38.variant_call.all.vcf -b Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.regions_${i}.bed -u | wc -l | tr -d "\n"
            echo -e "\t#variants"
        done >> $out
    done

### Compute average conservation score

    for resolution in 1000 10000 100000
    do
        for i in Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.regions_?.bed
        do
          for j in Acet BivProm DNase EnhA EnhWk GapArtf HET PromF Quies ReprPC TSS Tx TxEnh TxEx TxWk znf
          do
            echo -n $i" "$j" "
            bedtools intersect -a Data/Annotations/hg38_genome_100_segments.${j}.bed -b $i | awk 'BEGIN{s = 0}{s += $3 - $2}END{print s}'
          done
        done > Results/${resolution}/chromHmmResults.txt
    done

### Compute number of overlaps with exons

    for resolution in 1000 10000 100000
    do
        out=Results/${resolution}/n_genes.tsv
        echo -e "category\tscore\ttype" > $out
        for i in $( seq 0 7 )
        do
            echo -e -n $i "\t"
            awk '$10 > 0' Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.regions_${i}.coverage.refGene.bed | wc -l | tr -d "\n"
            echo -e "\t#genes"
        done >> $out
    done

### Compute conservation scores

    for resolution in 1000 10000 100000
    do
        for i in Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.regions_?.bed
        do
          multiBigwigSummary BED-file --bwfiles Data/Annotations/hg38.phyloP100way.bw --BED $i -o tmp.npz --outRawCounts ${i%bed}phyloP.tsv
        done
    done

### Run PhastCons

    zgrep '^>' Data/Genomes/hprc.y1v2.fa.gz  | cut -c 2- > Data/Genomes/hprc.y1v2.chr.listset
    # Drop haplotypes
    rev Data/Genomes/hprc.y1v2.chr.listset | cut -f 2-3 -d "#" | rev | sort -u > Data/Genomes/hprc.y1v2.listset
    # Extract individual genomes
    while IFS= read -r line
        do
        grep $line Data/Genomes/hprc.y1v2.chr.listset > Data/Genomes/selected.listset
        seqtk subseq Data/Genomes/hprc.y1v2.fa.gz Data/Genomes/selected.listset | pigz -c -p 10 > Data/Genomes/hprc.y1v2.${line}.fa.gz
    done < Data/Genomes/hprc.y1v2.listset
    # Compute mash distance
    mashtree --numcpus 12 Data/Genomes/hprc.y1v2.*.fa.gz --outmatrix Data/Genomes/phylip.tsv > Data/Genomes/mashtree.dnd
    # Convert HAL 2 MAF
    cactus-hal2maf ./js Data/Alignments/hprc-v1.1-mc-grch38-full.hal Data/Alignments/hprc-v1.1-mc-grch38-full.maf.gz --noAncestors --refGenome GRCh38 --filterGapCausingDupes --chunkSize 100000 --batchCores 10 --batchCount 10 --noAncestors --batchParallelTaf 10 --batchSystem slurm --logFile Data/Alignments/hprc-v1.1-mc-grch38-full.maf.gz.log --caching=false
    # Extract alignment of chr1 
    grep -A 5 -B 5 -m 2 -n "maf version=1 scoring=N/A" Data/Alignments/hprc-v1.1-mc-grch38-full.maf
    head -11265399 Data/Alignments/hprc-v1.1-mc-grch38-full.maf > Data/Alignments/hprc-v1.1-mc-grch38-chr1.maf
    # Run a script to fix names
    fixNames.py
    # Run PhyloFit to get neutral model
    phyloFit --tree Data/Genomes/mashtree_fix.dnd --msa-format MAF Data/Alignments/hprc-v1.1-mc-grch38-chr1_fix.maf > Data/Alignments/neutralmodel.mod
    # Run a script to split MAF file to chromosomes
    


### Compute ChromHmm

    for resolution in 1000 10000 100000
    do
        for i in Results/${resolution}/chr_all_hprc-v1.1-mc-chm13-full_GRCh38.0.chrall.regions_?.bed
        do
          for j in Acet BivProm DNase EnhA EnhWk GapArtf HET PromF Quies ReprPC TSS Tx TxEnh TxEx TxWk znf
          do
            echo -n $i" "$j" "
            bedtools intersect -a Data/Annotations/hg38_genome_100_segments.${j}.bed -b $i | awk 'BEGIN{s = 0}{s += $3 - $2}END{print s}'
          done
        done > Results/${resolution}/chromHmmResults.txt
    done

### Plot figures

    for resolution in 1000 10000 100000
    do
        Rscript gatherAll.R Results/${resolution}/variants.tsv Results/${resolution}/phyloP.tsv Results/${resolution}/n_genes.tsv Results/${resolution}/chromHmm.tsv '#variants,conservation,#genes,HET,Tx' Results/${resolution}/scores.png
    done
