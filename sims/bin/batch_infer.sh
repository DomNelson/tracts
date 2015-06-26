set -e

POPDIR=$1
OUTDIR=$2

for f in $POPDIR*.pkl
do
	TIMESTAMP=$(date +%Y%m%d)_$(date +%H)_$(date +%M)_$(date +%S)
	NEWDIR=$OUTDIR$TIMESTAMP

	mkdir $NEWDIR

	git rev-parse HEAD >> $NEWDIR/git_info.txt
	git diff >> $NEWDIR/git_info.txt

	echo Inference on Simons migmats  >> $NEWDIR/notes.txt

	POPFILE=$f

	echo Parameters: >> $NEWDIR/notes.txt
	echo $POPFILE >> $NEWDIR/notes.txt
	echo $NEWDIR >> $NEWDIR/notes.txt

	python taino_ppx_xxp_validation.py $POPFILE $NEWDIR/ &
	sleep 2
done

wait
