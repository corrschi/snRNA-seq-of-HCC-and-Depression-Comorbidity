# Script: mapping_mouse_mm10_batch.sh
#The mouse reference genome (mm10, refdata-gex-mm10-2020-A) was used for alignment
# https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

### Step 1：softlink Create
mkdir -p /public/home/chidm/Workspace/rawdata_softlink

# Define raw data and link paths
RAW=/public/home/chidm/Workspace/rawData/mdd.liver
LINK=/public/home/chidm/Workspace/rawdata_softlink

# Sample names (with underscores)
samples=(Control_1 Control_2 Control_3 Stress_1 Stress_2 Stress_3)

# Create symlinks for Cell Ranger input
for s in "${samples[@]}"; do
ln -sf "$RAW/$s/${s}_R1.fq.gz" "$LINK/${s}_S1_L001_R1_001.fastq.gz"
ln -sf "$RAW/$s/${s}_R2.fq.gz" "$LINK/${s}_S1_L001_R2_001.fastq.gz"
done

cd /public/home/chidm/Workspace/mapping/cDNA

# step2. mapping.lsf template

BSUB -J AAA
BSUB -n 16
BSUB -o %J.our
BSUB -e %J.err
BSUB -R span[hosts=1]
BSUB -q normal

/public/home/chidm/software/bin/cellranger6 count \
--id=__SAMPLE__ \
--transcriptome=/public/home/chidm/reference_data/refdata-gex-mm10-2020-A \
--fastqs=/public/home/chidm/Workspace/rawdata_softlink \
--sample=__SAMPLE__ \
--localcores=__THREAD__ \
--force-cells=13000

# step3. Batch submission
THREAD=16
for s in Control_1 Control_2 Control_3 Stress_1 Stress_2 Stress_3; do
sed -e "s#__SAMPLE__#$s#g" -e "s#__THREAD__#$THREAD#g" mapping.lsf | bsub
done


