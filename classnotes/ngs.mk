trim: SRR765989_F1.trim.fastq SRR765989_F2.trim.fastq
SRR765989_F1.trim.fastq SRR765989_F2.trim.fastq: SRR765989_F1.fastq SRR765989_F2.fastq
	java -jar ../Trimmomatic-0.38/trimmomatic-0.38.jar PE \
SRR765989_F1.fastq  SRR765989_F2.fastq \
SRR765989_F1.trim.fastq SRR765989_F1.trim.u.fastq \
SRR765989_F2.trim.fastq SRR765989_F2.trim.u.fastq  \
ILLUMINACLIP:../Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:2 \
TRAILING:2 \
SLIDINGWINDOW:4:15 \
MINLEN:30

index_ref: chr20.fa.bwt
chr20.fa.bwt: chr20.fa
	bwa index chr20.fa

align: SRR765989.sam
SRR765989.sam: SRR765989_F1.fastq  SRR765989_F2.fastq chr20.fa.bwt
	bwa mem chr20.fa SRR765989_F1.fastq SRR765989_F2.fastq >SRR765989.sam

mapped_reads: SRR765989.mapped.sam SRR765989.unmapped.sam
SRR765989.mapped.sam SRR765989.unmapped.sam: SRR765989.sam
	samtools view -F 4  SRR765989.sam > SRR765989.mapped.sam
	samtools view -f 4 SRR765989.sam > SRR765989.unmapped.sam

bam: SRR765989.bam
SRR765989.bam: SRR765989.sam
	samtools view -bS SRR765989.sam >SRR765989.bam

sort: SRR765989.sorted.bam
SRR765989.sorted.bam: SRR765989.bam
	samtools sort SRR765989.bam >SRR765989.sorted.bam

rm_duplicate: SRR765989.dup.bam
SRR765989.dup.bam: SRR765989.sorted.bam
	java -jar ../picard.jar MarkDuplicates INPUT=SRR765989.sorted.bam OUTPUT=SRR765989.dup.bam\
    METRICS_FILE=picard_metrics.txt VALIDATION_STRINGENCY=LENIENT

add_rg: SRR765989.dup.rg.bam
SRR765989.dup.rg.bam: SRR765989.dup.bam
	java -jar ../picard.jar AddOrReplaceReadGroups I=SRR765989.dup.bam \
O=SRR765989.dup.rg.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 \
RGSM=20 CREATE_INDEX=true

index_fasta: chr20.fa.fai
chr20.fa.fai: chr20.fa
	samtools faidx chr20.fa

dict: chr20.dict
chr20.dict: chr20.fa
	java -jar ../picard.jar CreateSequenceDictionary REFERENCE=chr20.fa OUTPUT=chr20.dict


realign_target: SRR765989.paired.bam.list
SRR765989.paired.bam.list: SRR765989.dup.rg.bam
	java -jar gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R chr20.fa -o SRR765989.paired.bam.list -I SRR765989.dup.rg.bam

realign: SRR765989.realigned.bam
SRR765989.realigned.bam: SRR765989.paired.bam.list
	java -jar gatk/GenomeAnalysisTK.jar -I SRR765989.dup.rg.bam -R chr20.fa \
-T IndelRealigner -targetIntervals SRR765989.paired.bam.list -o SRR765989.realigned.bam

vcf: raw_variants.vcf
raw_variants.vcf: SRR765989.realigned.bam
	java -jar gatk/GenomeAnalysisTK.jar  -T HaplotypeCaller -R chr20.fa  \
	-I SRR765989.realigned.bam --genotyping_mode DISCOVERY -stand_call_conf 30 \
	 -o raw_variants.vcf

annotate: SRR765989.snpeff.vcf
SRR765989.snpeff.vcf: raw_variants.vcf
	java -Xmx2g -Djava.io.tmpdir=. -jar  snpEff/snpEff.jar -v hg19 raw_variants.vcf \
	>SRR765989.snpeff.vcf
