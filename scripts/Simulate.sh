set -e

MIGFILE=$1
NUMINDS=$2
METHOD=$3
OUTDIR=$4
POPNAME=$5

TIMESTAMP=$(date +%Y%m%d)_$(date +%H)_$(date +%M)_$(date +%S)
NEWDIR=$OUTDIR$TIMESTAMP
mkdir $NEWDIR

git describe HEAD > $NEWDIR/git_info.txt
git diff >> $NEWDIR/git_info.txt

cp Simulate.sh $NEWDIR/Simulate.sh

echo Forward simulations on a group of migration matrices as requested by Simon  > $NEWDIR/notes.txt

PLOTOUT=$NEWDIR/${POPNAME}_plot.png
POPOUT=$NEWDIR/${POPNAME}_pop.pkl
BEDOUT=$NEWDIR/$POPNAME/

echo Parameters: >> $NEWDIR/notes.txt
echo $MIGFILE >> $NEWDIR/notes.txt
echo $NUMINDS >> $NEWDIR/notes.txt
echo $BEDOUT >> $NEWDIR/notes.txt
echo $POPOUT >> $NEWDIR/notes.txt
echo $PLOTOUT >> $NEWDIR/notes.txt

python tracts_sim.py $MIGFILE $NUMINDS $METHOD $BEDOUT $POPOUT $PLOTOUT
