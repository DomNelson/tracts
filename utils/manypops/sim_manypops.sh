set -e

MIGFILE=$1
NUMINDS=$2
OUTDIR=$3
POPNAME=$4
NUMPOPS=$5

TIMESTAMP=$(date +%Y%m%d)_$(date +%H)_$(date +%M)_$(date +%S)
NEWDIR=$OUTDIR$TIMESTAMP
mkdir $NEWDIR

git describe HEAD > $NEWDIR/git_info.txt
git diff >> $NEWDIR/git_info.txt

cp $(basename $0) $NEWDIR/$(basename $0)

echo Forward simulations on a group of migration matrices as requested by Simon > $NEWDIR/notes.txt

PLOTOUT=$NEWDIR/${POPNAME}_plot.png
POPOUT=$NEWDIR/${POPNAME}_pop.pkl
BEDOUT=$NEWDIR/$POPNAME/

echo Parameters: >> $NEWDIR/notes.txt
echo $MIGFILE >> $NEWDIR/notes.txt
echo $NUMINDS >> $NEWDIR/notes.txt
echo $BEDOUT >> $NEWDIR/notes.txt
echo $POPOUT >> $NEWDIR/notes.txt
echo $PLOTOUT >> $NEWDIR/notes.txt

cd ~/project/tracts

for ((i=1; i<=$NUMPOPS; i++)); do
	echo Simulating $POPNAME population $i of $NUMPOPS
	BEDOUTNUM=$BEDOUT$i/
	POPOUTNUM=None
	PLOTOUTNUM=$NEWDIR/${POPNAME}/${POPNAME}${i}_plot.png

	python tracts_sim.py $MIGFILE None None $NUMINDS forward $BEDOUTNUM $POPOUTNUM $PLOTOUTNUM
done
