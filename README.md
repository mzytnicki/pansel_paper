# Analysis for the `pansel` paper

## Architecture

    - Data
      - Graphs
        - MGC # MiniGraph-Cactus
        - PGGB
        - M_xanthus
      - Annotations
      - Genomes
      - Alignments
    - Results
      - MGC
        - 1000
        - 10000
        - 100000
      - PGGB
        - 1000
        - 10000
        - 100000
      - M_xanthus
        
## Download data

### Graphs

#### MiniGraph-Cactus

    cd Data/Graphs/MGC
    for i in `seq 1 22`
    do
      wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2022_03_11_minigraph_cactus/chrom-graphs-hprc-v1.1-mc-chm13-full/chr${i}.vg
      vg convert -fW chr${i}.vg | gzip -c > chr${i}.gfa.gz
      rm chr${i}.vg
    done

#### PGGB

    cd Data/Graphs/PGGB
    for i in `seq 1 22`
    do
      wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_07_30_pggb/chroms/chr${i}.pan.fa.a2fb268.e820cd3.9ea71d8.smooth.gfa.gz && mv chr${i}.pan.fa.a2fb268.e820cd3.9ea71d8.smooth.gfa.gz chr${i}.gfa.gz
    done

### Create graph for the *M_xanthus*

    singularity exec docker://quay.io/comparative-genomics-toolkit/cactus:latest cactus-pangenome ./js ./seqfiles.txt --outDir ./Data/Graphs/M_xanthus --outName m_xanthus --reference GCA_000012685.1 --mgCores 12 --mapCores 12 --consCores 12 --gfa full

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

*M. xanthus* annotation

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/012/685/GCA_000012685.1_ASM1268v1/GCA_000012685.1_ASM1268v1_genomic.gtf.gz
    
## Analysis with MiniGraph-Cactus

### Compute pansel results for different resolutions

    for resolution in 1000 10000 100000
    do
        for i in `seq 1 22`
        do
          sed "s/^/chr${i}\t/g" Data/Graphs/chr${i}.tsv
        done > Results/MGC/${resolution}/chrall.tsv
    done

### Exclude gaps, and divide into different 8 bins (most conserved, to least conserved)

    for resolution in 1000 10000 100000
    do
        python3 divideRegions.py Results/MGC/${resolution}/chrall.tsv ${resolution} Data/Annotations/hg38.gaps.bed 8 Results/MGC/${resolution}/chrall.regions
    done

### Transform `pansel` results to BED format

    for binSize in 1000 10000 100000
    do
        for i in `seq 1 22`
        do
            sed "s/^/chr${i}\t/g" Results/MGC/${binSize}/chr${i}.tsv
        done | awk '{print $1 "\t" ($3-1) "\t" $4 "\tbin_" NR "\t" ($5 / ($4 - $3)) "\t+"}' > Results/MGC/${binSize}/chrall.bed
    done

### Fit to conserved/divergent distributions

    for binSize in 1000 10000 100000
    do
        Rscript getExtremes.R -i Results/MGC/${binSize}/chrall.bed -p 0.05 -P 0.05 -t Results/MGC/${binSize}/fit.png -o Results/MGC/${binSize}/fitConserved.bed -O Results/MGC/${binSize}/fitDivergent.bed &> Results/MGC/${binSize}/fit.log
    done

### Compute number of overlaps with structural variations

    for resolution in 1000 10000 100000
    do
        out=Results/MGC/${resolution}/variants.tsv
        echo -e "category\tscore\ttype" > $out
        for i in $( seq 0 7 )
        do 
            echo -e -n "$i\t"
            bedtools intersect -a Data/Annotations/GRCh38.variant_call.all.vcf -b Results/MGC/${resolution}/chrall.regions_${i}.bed -u | wc -l | tr -d "\n"
            echo -e "\t#variants"
        done >> $out
    done

### Compute average conservation score

    for resolution in 1000 10000 100000
    do
        for i in Results/MGC/${resolution}/chrall.regions_?.bed
        do
          for j in Acet BivProm DNase EnhA EnhWk GapArtf HET PromF Quies ReprPC TSS Tx TxEnh TxEx TxWk znf
          do
            echo -n $i" "$j" "
            bedtools intersect -a Data/Annotations/hg38_genome_100_segments.${j}.bed -b $i | awk 'BEGIN{s = 0}{s += $3 - $2}END{print s}'
          done
        done > Results/MGC/${resolution}/chromHmmResults.txt
    done

### Compute number of overlaps with exons

    for resolution in 1000 10000 100000
    do
        out=Results/MGC/${resolution}/n_genes.tsv
        echo -e "category\tscore\ttype" > $out
        for i in $( seq 0 7 )
        do
            echo -e -n $i "\t"
            awk '$10 > 0' Results/MGC/${resolution}/chrall.regions_${i}.coverage.refGene.bed | wc -l | tr -d "\n"
            echo -e "\t#genes"
        done >> $out
    done

### Compute conservation scores

    for resolution in 1000 10000 100000
    do
        for i in Results/MGC/${resolution}/chrall.regions_?.bed
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
    phyloFit --tree Data/Genomes/mashtree_fix.dnd --msa-format MAF Data/Alignments/hprc-v1.1-mc-grch38-chr1_fix.maf > Data/Alignments/phyloFit.mod
    # Run a script to split MAF file to chromosomes
    splitMaf.py
    # Run PhastCons (PhyloP did not work, see https://github.com/CshlSiepelLab/phast/issues/69)
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 MT X Y
    do
        phastCons Data/Alignments/hprc-v1.1-mc-grch38-chr${i}.maf Data/Alignments/phyloFit.mod > Data/Alignments/phyloPscores_${i}.wig
    done
    # Transform WIG to BigWig
    for i in $(seq 22); do echo sed "s/(null)/chr$i/g" Data/Alignments/phyloPscores_${i}.wig; done > Data/Alignments/phyloPscores.wig
    wigToBigWig Data/Alignments/phyloPscores.wig Data/Genomes/hs.size Data/Alignments/phyloPscores.bw # Requires 40GB RAM

    for resolution in 1000 10000 100000
    do
        for i in Results/MGC/${resolution}/chrall.regions_?.bed
        do
          out=${i%bed}phyloP_human.tsv
          rm -f $out
          multiBigwigSummary BED-file --bwfiles Data/Alignments/phyloPscores.bw --BED $i -o tmp.npz --outRawCounts ${out}
        done
    done

### PSTs

    for resolution in 1000 10000 100000
    do
        out=Results/MGC/${resolution}/conserved.tsv
        echo -e "category\tscore\ttype" > $out
        for i in $( seq 0 7 )
        do 
            echo -e -n "$i\t"
            bedtools intersect -a Data/Annotations/pst-31mer.grch38.bed.gz -b Results/MGC/${resolution}/chrall.regions_${i}.bed -u | wc -l | tr -d "\n"
            echo -e "\tPST"
        done >> $out
    done

### Use Corer (unsuccessful)

    ls -l Data/Genomes/hprc.y1v2.*.fa.gz | tr -s ' ' | cut -d' ' -f9 > Data/Genomes/list.txt
    Bifrost build -r Data/Genomes/list.txt -o Results/bifrost -c -t 20  # Requires 160GB RAM
    Corer -i Results/bifrost.gfa.gz -c Results/bifrost.color.bfg -o Results/core -q 90 -d 60 -t 80 # This takes more than 4 days
        
### Compute ChromHmm

    for resolution in 1000 10000 100000
    do
        for i in Results/MGC/${resolution}/chrall.regions_?.bed
        do
          for j in Acet BivProm DNase EnhA EnhWk GapArtf HET PromF Quies ReprPC TSS Tx TxEnh TxEx TxWk znf
          do
            echo -n $i" "$j" "
            bedtools intersect -a Data/Annotations/hg38_genome_100_segments.${j}.bed -b $i | awk 'BEGIN{s = 0}{s += $3 - $2}END{print s}'
          done
        done > Results/MGC/${resolution}/chromHmmResults.txt
    done

### Plot score figures

    for resolution in 1000 10000 100000
    do
        Rscript gatherAll.R Results/MGC/${resolution}/variants.tsv Result/MGC/${resolution}/phyloP.tsv Results/MGC/${resolution}/phyloP_human.tsv Results/MGC/${resolution}/conserved.tsv Results/MGC/${resolution}/n_genes.tsv Results/MGC/${resolution}/chromHmm.tsv '#variants,conservation,cons. HS,PST,#genes,HET,Tx' Results/MGC/${resolution}/scores.png
    done

### Extract subgraph for Bandage-NG (for NBPF20, aka ENSG00000162825.18)

    # Find genomic positions
    grep ENSG00000162825.18 Data/Annotations/gencode.v44.annotation.gtf | awk '$3 == "gene"'
    # Convert the right gfa to og format
    odgi build -g Data/Graphs/chr_1_hprc-v1.1-mc-chm13-full.gfa -o Data/Graphs/chr_1_hprc-v1.1-mc-chm13-full.og -P -t 10
    # Extract the sub-graph around the gene
    odgi extract -i Data/Graphs/chr_1_hprc-v1.1-mc-chm13-full.og -r GRCh38.0.chr1:145289900-145425603 -o tmp.og -t 10 -P
    # Remove spurious _MINIGRAPH_ paths
    odgi view -i tmp.og -g -t 10 -P | grep -v MINIGRAPH > tmp.gfa
    # Convert again to og
    odgi build -g tmp.gfa -o tmp.og -P -t 10 -s -O
    # Transform exon GFF lines to BED
    grep ENSG00000162825.18 Data/Annotations/gencode.v44.annotation.gtf | gff2bed | awk '$8 == "exon"' | cut -f 1-6 | sort -u | awk '{$1 = "GRCh38.0.chr1"; $4 = "ENSG00000162825.18_" NR; print}' | tr " " "\t" > tmp.bed
    # Adapt bed, so that it is on the right path
    odgi procbed -i tmp.og -b tmp.bed > tmp.proc.bed
    # Set the bed as a path
    odgi inject -i tmp.og -b tmp.proc.bed -o tmp2.og
    # Sort, and visualize
    odgi sort -i tmp2.og -o - -O | odgi viz -i - -o tmp.png -s '#' -P -t 10
    odgi view -i tmp2.og -g > Results/NBPF20.gfa
    
## Comparison between MiniGraph-Cactus and PGGB

    Rscript compareResults.R Results Results/comparisonMGC_PGGB.png &> Results/comparisonMGC_PGGB.log

