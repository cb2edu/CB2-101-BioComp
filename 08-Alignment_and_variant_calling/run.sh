# Generate index chr20
/share/apps/bwa-0.7.10/bwa index chr20.fa -p chr20
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
