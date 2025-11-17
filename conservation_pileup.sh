
### CHANGE ME
WRK=/Path/to/Fox_NFIA_CTCF/
OUTDIR=/home/ubuntu/honors_research/conservation_analysis
###

# Dependencies
# - java
# - perl
# - python

set -exo
module load anaconda
#source activate bx

# Fill in placeholder constants with your directories
#MOTIF=$WRK/data/RefPT-Motif/1000bp
CDIR=/home/ubuntu/honors_research/conservation_analysis
#pip install pyBigWig   

# Script shortcuts
SCRIPTMANAGER=/home/ubuntu/ScriptManager-v0.14.jar
COMPOSITE=/home/ubuntu/honors_research/conservation_analysis/sum_Col_CDT.pl
#WIGTOBG=$WRK/bin/convert_wig_to_bedgraph.py
PILEUPBG=/home/ubuntu/honors_research/conservation_analysis/pileup_BedGraph_on_RefPT.py
PILEUPBW=/home/ubuntu/honors_research/conservation_analysis/pileup_BigWig_on_RefPT.py
#SPILEUPBW=$WRK/bin/pileup_BigWig_on_RefPT_stranded.py
SUMMAT=/home/ubuntu/honors_research/conservation_analysis/sum_each_CDT.py

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Define set of BED files to pileup on (this is nucleosome!!)
BEDFILE="/home/ubuntu/honors_research/peak_align/reference_files/edited_TSS_all_adj+1Nuc_4_250bp.bed"
BED=`basename $BEDFILE ".bed"`

[ -d $OUTDIR/$BED ] || mkdir $OUTDIR/$BED
[ -d $OUTDIR/$BED/CDT ] || mkdir $OUTDIR/$BED/CDT
[ -d $OUTDIR/$BED/Composites ] || mkdir $OUTDIR/$BED/Composites
[ -d $OUTDIR/$BED/PNG ] || mkdir $OUTDIR/$BED/PNG
[ -d $OUTDIR/$BED/SVG ] || mkdir $OUTDIR/$BED/SVG

# Pileup conservation scores
for CONSERVATION in "/home/ubuntu/honors_research/conservation_analysis/hg38.phyloP30way.bw";
do
	CONS=`basename $CONSERVATION ".bw"`

	echo $BED x $CONS

	# Pileup BigWig
	python3 $PILEUPBW -i $CONSERVATION -r $BEDFILE -o $OUTDIR/$BED/CDT/$CONS\_$BED.cdt
	# Make composite
	perl $COMPOSITE $OUTDIR/$BED/CDT/$CONS\_$BED.cdt $OUTDIR/$BED/Composites/$CONS\_$BED.out

	# Count sites
	#NSITES=`wc -l $BEDFILE | awk '{print $1-1}'`

	# Two-color heatmap
	#java -jar $SCRIPTMANAGER figure-generation heatmap -p 0.95 --black $OUTDIR/$BED/CDT/$CONS\_$BED.cdt -o $OUTDIR/$BED/PNG/$CONS\_$BED.png
	#java -jar $SCRIPTMANAGER figure-generation label-heatmap $OUTDIR/$BED/PNG/$CONS\_$BED.png \
		#-l "-250" -m "0" -r "+250" -w 1 -f 20 \
		#-x $BED -y "$BED occurences (${NSITES} sites)" \
		#-o $OUTDIR/$BED/SVG/$CONS\_$BED.svg
done