set -e

MIGDIR=$1
OUTDIR=$2
NUMINDS=$3

TIMESTAMP=$(date +%Y%m%d)_$(date +%H)_$(date +%M)_$(date +%S)
NEWDIR=$OUTDIR$TIMESTAMP

mkdir $NEWDIR

git rev-parse HEAD > $NEWDIR/git_info.txt
git diff >> $NEWDIR/git_info.txt

cp simulate.sh $NEWDIR/simulate.sh

echo Forward simulations on a group of migration matrices as requested by Simon  > $NEWDIR/notes.txt

for f in $MIGDIR*
do
	MIGFILE=$f
	PLOTOUT=$NEWDIR/$(basename "$f")_plot.png
	POPOUT=$NEWDIR/$(basename "$f")_pop.pkl

	echo Parameters: >> $NEWDIR/notes.txt
	echo $MIGFILE >> $NEWDIR/notes.txt
	echo $PLOTOUT >> $NEWDIR/notes.txt
	echo $NUMIND >> $NEWDIR/notes.txt
	echo $POPOUT >> $NEWDIR/notes.txt

	python tracts_sim.py $MIGFILE $PLOTOUT $NUMINDS $POPOUT &
	sleep 2
done

wait
