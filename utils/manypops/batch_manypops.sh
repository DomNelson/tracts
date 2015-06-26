set -e

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Argentina_mig 154 ~/project/tracts/sims/results/simon_batch/newmatrices/ Argentina 100 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Chile_mig 27 ~/project/tracts/sims/results/simon_batch/newmatrices/ Chile 100 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Colombia_mig 95 ~/project/tracts/sims/results/simon_batch/newmatrices/ Colombia 100 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Ecuador_mig 19 ~/project/tracts/sims/results/simon_batch/newmatrices/ Ecuador 100 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Peru_mig 117 ~/project/tracts/sims/results/simon_batch/newmatrices/ Peru 100 &