setwd("/sc/arion/projects/breen_lab/AD_editing/MSBB/BM_10/scripts")

file_list=read.table("samples.txt", stringsAsFactors=FALSE) # a text file with names of all unmapped bam files
files=basename(file_list$V1)
full_run_input <- NULL
full_run_input <- paste('#!/bin/bash\n', sep="")
for(item in 1:length(files)){
        file_name1=(files[item])
        sh_file <- paste(file_name1,"_config.sh", sep="")
        full_run_input <- paste(full_run_input,'bsub < ', sh_file, ';\n sleep 1;\n', sep="")
        sh_runfile <- "run_all_samples.sh"
        fileConn <- file(sh_file)

writeLines(paste('
#!/bin/bash
#BSUB -W 5:00
#BSUB -n 8
#BSUB -q express
#BSUB -P acc_breen_lab
#BSUB -cwd /sc/arion/projects/breen_lab/AD_editing/MSBB/BM_10/scripts
#BSUB -o /sc/arion/projects/breen_lab/AD_editing/MSBB/BM_10/scripts/logs/',file_name1,'.log.out
#BSUB -e /sc/arion/projects/breen_lab/AD_editing/MSBB/BM_10/scripts/logs/',file_name1,'.log.err
#BSUB -u alexander.chung@icahn.mssm.edu
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -L /bin/bash

#load all modules required...
#module load samtools/1.1            #this will query all recoding sites
#module load star/2.7.3a         #this does the mapping
#module load subread             #this does the counting
#module load rnaeditingindexer   #this generates the AEI
#module load picard              #load picard tools for QC                   

#move into local directory...
cd /sc/arion/projects/breen_lab/AD_editing/MSBB/BM_10/editing_sites

#Map FASTQ files to generate bam files
#STAR --genomeDir /sc/arion/projects/H_PBG/REFERENCES/GRCh38/star/2.7.1a/Gencode.v30.overhang100/chr_primary --sjdbGTFfile /sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.primary_assembly.annotation.gtf --runThreadN 12 --twopassMode None --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD --outFileNamePrefix ',file_name1,'. --readFilesCommand zcat  --readFilesIn ',file_name1,'_1.fastq.gz ',file_name1,'_2.fastq.gz --limitBAMsortRAM 71361722155

#Input bam file and generate RNA-seq counts for each gene (#if stranded use -s1 or -s2, depending on directionality. otherwise, if unstranded use -s0)
#featureCounts -T 5 -p -t exon -s 0 -g gene_id -a /sc/arion/projects/H_PBG/REFERENCES/GRCh38/Gencode/release_30/gencode.v30.primary_assembly.annotation.gtf ', file_name1,'.Aligned.sortedByCoord.out.bam -o ', file_name1,'.counts.txt

#getting halted before removal of temporary files 
#This removes temp files left over from STAR
#rm ', file_name1,'.Unmapped.out.mate1
#rm ', file_name1,'.Unmapped.out.mate2
#rm -r ',file_name1,'._STARgenome
#rm ',file_name1,'.Log.progress.out
#rm ',file_name1,'.Log.out
#rm ',file_name1,'.SJ.out.tab

#This generates PICARD QC
#java -jar $PICARD CollectAlignmentSummaryMetrics R=/sc/arion/projects/breen_lab/references/hg38.fa I= ', file_name1,'.Aligned.sortedByCoord.out.bam O= ', file_name1, '_picard.txt

#load all modules required...
module load samtools/1.1            #this will query all recoding sites
module load subread             #this does the counting
export PERL5LIB=$PERL5LIB:/sc/arion/projects/breen_lab/RecodingProject/BrainSpan/

perl /sc/arion/projects/breen_lab/RecodingProject/BrainSpan/query_known_sites.pl /sc/arion/projects/breen_lab/RecodingProject/BrainSpan/CNS_REDi_combined.txt /sc/arion/projects/adineto/mayo/bam_files/MayoRNAseq_Cerebellum_BAMs/',file_name1,'.aligned.sort.coord.bam  /sc/arion/projects/breen_lab/AD_editing/MAYO/CER/editing_sites/',file_name1,'_edits.txt

module load rnaeditingindexer   #this generates the AEI
mkdir /sc/arion/projects/breen_lab/AD_editing/MAYO/CER/AEI/',file_name1,'_AEI
cp /sc/arion/projects/adineto/mayo/bam_files/MayoRNAseq_Cerebellum_BAMs/',file_name1,'.aligned.sort.coord.bam /sc/arion/projects/breen_lab/AD_editing/MAYO/CER/AEI/',file_name1,'_AEI
cd /sc/arion/projects/breen_lab/AD_editing/MAYO/CER/AEI/',file_name1,'_AEI

RNAEditingIndex -d . -f ',file_name1,'.aligned.sort.coord.bam -o . --genes_expression /sc/arion/projects/breen_lab/AEI/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz --refseq /sc/arion/projects/breen_lab/AEI/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz --snps /sc/arion/projects/breen_lab/AEI/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150v2.bed.gz -gf /sc/arion/projects/breen_lab/AEI/Genomes/HomoSapiens/ucscHg38Genome.fa -rb /sc/arion/projects/breen_lab/AEI/Regions/HomoSapiens/ucscHg38Alu.bed.gz --genome UserProvided --paired_end

#This renames the original output file (EditingIndex.csv) and gives it a new title 
mv EditingIndex.csv  ',file_name1,'_AEI.csv
mv ',file_name1,'_AEI.csv /sc/arion/projects/breen_lab/AD_editing/MAYO/CER/AEI


echo Finished running ', file_name1, '
        ', sep=""), fileConn)
        close(fileConn)

}

fileConn <- file(sh_runfile)
writeLines(paste(full_run_input, sep=""), fileConn)
close(fileConn)                                                                                                                               
