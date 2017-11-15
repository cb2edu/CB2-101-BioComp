### NGS code for CB2-101
### Variant detection

# Download data file
wget https://github.com/cb2edu/CB2-101-BioComp/raw/master/08-Alignment_and_variant_calling/data/cb2-101-variant_detection_sample_data.tar.xz

# Download bwa
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbio-bwa%2Ffiles%2F&ts=1479318023&use_mirror=superb-dca2

#Download samtools
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2

# Install some needed software
su -c "yum install ncurses ncurses-devel"

# Donwload GATK
wget -O Gatk.tar.bz2 https://software.broadinstitute.org/gatk/download/auth?package=GATK

# Get picard
wget https://github.com/broadinstitute/picard/releases/download/2.15.0/picard.jar

# Put all the software in path


####################### Alignment #############################

# Step 1 create BWA index
bwa index chr20.fa

# Align the sequence agains the reference
bwa mem chr20.fa SRR765989_F1.fastq SRR765989_F2.fastq >SRR765989.sam

# Convert sam to bam files
samtools view -bS SRR765989.sam >SRR765989.bam

# Sort the bam files
samtools sort SRR765989.bam >SRR765989.sorted.bam



## Generate index chr20
/share/apps/bwa-0.7.10/bwa index chr20.fa -p chr20

## Convert sam to bam
samtools view -b SRR765989.sam >SRR765989.bam

## Sort the bam file based on coordinate
samtools sort SRR765989.bam >SRR765989.sorted.bam


################### Variant detection ################################
## Mark duplicates
java -jar picard.jar MarkDuplicates INPUT=SRR765989.sorted.bam OUTPUT=SRR765989.dup.bam\
    METRICS_FILE=picard_metrics.txt VALIDATION_STRINGENCY=LENIENT

# Create dictionary for the reference genome
java -jar ../../picard-tools-1.140/picard.jar CreateSequenceDictionary REFERENCE=chr20.fa OUTPUT=chr20.dict

# Create an index of the reference
samtools faidx chr20.fa

# Add readgroups to bam file
# For the description of readgroups look here:
# https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
java -jar ../picard-tools-1.140/picard.jar AddOrReplaceReadGroups I=SRR765989.dup.bam O=SRR765989.dup.rg.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 CREATE_INDEX=true

## create realignment target
java -jar ../GenomeAnalysisTK.jar -T RealignerTargetCreator -R gatk_ref/chr20.fa -o SRR765989.paired.bam.list -I SRR765989.dup.rg.bam

# Realign
java -jar ../GenomeAnalysisTK.jar -I SRR765989.dup.rg.bam -R gatk_ref/chr20.fa -T IndelRealigner -targetIntervals SRR765989.paired.bam.list -o SRR765989.realigned.bam

# Create vcf
# Use Haplotypecaller for all cases
# For options see:
# https://software.broadinstitute.org/gatk/documentation/article?id=2803
# -stand_emit_conf shows error
java -jar GenomeAnalysisTK.jar  -T HaplotypeCaller -R chr20.fa  -I SRR765989.realigned.bam --genotyping_mode DISCOVERY -stand_call_conf 30 -o raw_variants.vcf

# Variant annotation using SNPEFF
## Download SNPEFF
wget -O snpeff.zip https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download

## Find all the available databases for annotation
java -Xmx2g -Djava.io.tmpdir=. -jar ../snpEff/snpEff.jar databases | grep -i Homo_sapiens

## Annotate VCF file
java -Xmx2g -Djava.io.tmpdir=. -jar snpEff.jar -v hg19 SRR765989.vcf >SRR765989.snpeff.vcf


#---------------------Old code------------------------------
# Generate samtools faidx
/share/apps/samtools/samtools faidx human/chr20.fa
# Generate dic file with picard
java -jar picard-tools-1.119/CreateSequenceDictionary.jar  REFERENCE=human/chr20.fa  OUTPUT=human/chr20.dict
# Align F1 and F2 with bwa
/share/apps/bwa-0.7.10/bwa aln -t 1 -f SRR765989_F1.fastq.sai  human/chr20 SRR765989_F1.fastq
/share/apps/bwa-0.7.10/bwa aln -t 1 -f SRR765989_F2.fastq.sai  human/chr20 SRR765989_F2.fastq
# Generate SAM file
/share/apps/bwa-0.7.10/bwa sampe human/chr20 SRR765989_F1.fastq.sai SRR765989_F2.fastq.sai SRR765989_F1.fastq SRR765989_F2.fastq  -f SRR765989.sam -r "@RG\tID:. .\tLB:. .\tSM:. .\tPL:ILLUMINA"
# Convert SAM to BAM
/share/apps/samtools/samtools view -T human/chr20.fa -bS SRR765989.sam > SRR765989.bam
# Sort BAM and index
/share/apps/samtools/samtools sort SRR765989.bam -o SRR765989.bam.sorted
/share/apps/samtools/samtools index SRR765989.bam.sorted SRR765989.bai
# Filetering out unmapped reads
/share/apps/samtools/samtools view -h -F 4 SRR765989.bam.sorted -o SRR765989.mapped.bam
# Select properly paired reads
/share/apps/samtools/samtools view -h -f 0X0002  SRR765989.mapped.bam -o SRR765989.paired.bam
# Mark duplicate
java -Xmx4g -Djava.io.tmpdir=/tmp -jar picard-tools-1.119/MarkDuplicates.jar INPUT=SRR765989.paired.bam  OUTPUT=SRR765989.marked.bam METRICS_FILE=metrics CREATE_INDEX=true  VALIDATION_STRINGENCY=LENIENT
# Local Realignment
java -Xmx4g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R human/chr20.fa -o SRR765989.paired.bam.list   -I SRR765989.marked.bam
# Indel issue
java -Xmx4g -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar -I SRR765989.marked.bam -R human/chr20.fa -T IndelRealigner -targetIntervals SRR765989.paired.bam.list -o SRR765989.marked.realigned.bam
# Picard realignment
java -Djava.io.tmpdir=/tmp/flx-auswerter -jar FixMateInformation.jar INPUT=SRR765989.marked.realigned.bam OUTPUT=SRR765989.marked.realigned.fixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
# Recalibration
java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator  -R human/chr20.fa  -I SRR765989.marked.realigned.fixed.bam -knownSites ../dbsnp_138.hg19.vcf -o SRR765989.recal_data.csv
java -jar GenomeAnalysisTK.jar -T PrintReads -R human/chr20.fa  -I SRR765989.marked.realigned.fixed.bam -BQSR SRR765989.recal_data.csv -o SRR765989.marked.realigned.fixed.recal.bam
# Call
java -Xmx4g -jar GenomeAnalysisTK.jar -glm BOTH -R human/chr20.fa -T UnifiedGenotyper -I SRR765989.marked.realigned.fixed.recal.bam -D ../dbsnp_138.hg19.vcf -o SRR765989.vcf -metrics snps.metrics -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 1000 -A AlleleBalance
# Recalibration
# Filtering
#java -jar GenomeAnalysisTK.jar -T ApplyRecalibration -R human/chr20.fa -input SRR765989.vcf -mode SNP -recalFile SRR765989.snps.recal -tranchesFile raw.SNPs.tranches -o SRR765989.recalibrated.vcf -ts_filter_level 99.0
#java -Xmx4g -jar GenomeAnalysisTK.jar -R human/chr20.fa -T VariantFiltration -V SRR765989.recalibrated.vcf -o SRR765989.recalibrated.filtered.vcf --clusterWindowSize 10
