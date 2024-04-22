#############
# ARGUMENTS #
#############

# 1: input directory (required)

########################
# Read pre-processment #
########################

echo "Starting analysis"
echo "Pre-processing the reads"
echo "------------------------"
echo ""

###################
# Adapter removal #
###################

echo ""
echo "1. Adapter removal"
echo "------------------"

time (for i in $1*.fastq
do
	porechop -t 32 -i $i -o ${i::-6}-porechop.fastq
done)

#################################################
# Length filtering pb ###########################
# May need to change max lenght for plasmopara ##
#################################################

echo ""
echo "2. Length filtering"
echo "-------------------"

time (for i in $1*porechop*
do
	cat $i | NanoFilt -l 400 > ${i::-6}-nanofilt.fastq
done)

######
# QC #
######

echo ""
echo "3. Quality check"
echo "----------------"

time(for i in $1*porechop-nanofilt*
do 
	NanoStat -t 32 --fastq $i > ${i::-6}-NanoStat.txt
done)

###################
# Chimera removal #
###################

echo ""
echo "4. Chimera removal"
echo "------------------"

time(for i in $1*nanofilt.fastq
do 
	minimap2 -x ava-ont -g 500 -t 32 $i $i > ${i%%.*}.paf
	yacrd -i ${i%%.*}.paf -o ${i%%.*}.yacrd -c 4 -n 0.4 scrubb -i $i -o ${i%%.*}.scrubb.fastq

	rm ${i%%.*}.paf
done)

####################
# Quality filtering#
####################

echo ""
echo "5. Quality filtering"
echo "--------------------"


time (for i in $1*scrubb.fastq*
do
	cat $i | NanoFilt -q 12  > ${i::-6}-scrubb.refilt.fastq
done)

###########################
# Read mapping (minimap2) #
###########################

#################
# Map sequences #
#################

echo ""
echo "6. Mapping"
echo "----------"

# Silva 138
time (for i in $1*scrubb.refilt.fastq
do echo "Mapping" $i; minimap2 \
-x map-ont \
-t 32 \
--secondary=no \
-K 10M /mnt/DATA_1/scripts/Minimap_pipeline/general_ITS_DB_V02.mmi $i > $i.paf
done)

# For understanding PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md

#############
# Filtering #
#############

echo ""
echo "7. Filtering the PAF files"
echo "--------------------------"

# Although secondary alignments are turned off, some query reads have been mapped to multiple database seqs.
# Let's filter the output file to just keep one alignment per read.
# We will keep the largest alignment.
# We will remove alignments <250 pb

time (for i in $1*.paf
do /mnt/DATA_1/scripts/Minimap_pipeline/filterPAF.py -i $i > ${i::-4}-f.paf
done)

rm -r $1filteredPAFs
mkdir $1filteredPAFs
mv $1*-f.paf $1filteredPAFs

########################################
# Summarize & Merge filtered PAF files #
########################################

echo ""
echo "8. Summarizing the PAF files"
echo "----------------------------"

# It creates a table with the number of sequences assigned to each Database's ID for each sample
# It's something like a QIIME (v1) OTU table.

/mnt/DATA_1/scripts/Minimap_pipeline/merfePAF.py -i $1filteredPAFs/ > $1otu_table_P.viticola_2016.csv

####################################
# Create a phyloseq taxonomy table #
####################################

echo ""
echo "9. Creating the taxonomy table"
echo "------------------------------"

# This script creates the taxonomy table needed for loading the data into the phyloseq package.

/mnt/DATA_1/scripts/Minimap_pipeline/taxonomyTable.py \
-i $1otu_table_P.viticola_2016.csv \
-t /mnt/DATA_1/scripts/Minimap_pipeline/ITS_taxonomy_V02.tsv > $1phyloseq_taxonomy_P.viticola_2016.csv

