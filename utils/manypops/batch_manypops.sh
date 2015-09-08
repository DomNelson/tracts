set -e

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Argentina_mig 154 ~/project/tracts/sims/results/simon_batch/newmatrices/1000pops/ Argentina 1000 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Chile_mig 27 ~/project/tracts/sims/results/simon_batch/newmatrices/1000pops/ Chile 1000 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Colombia_mig 95 ~/project/tracts/sims/results/simon_batch/newmatrices/1000pops/ Colombia1 500 &


sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Colombia_mig 95 ~/project/tracts/sims/results/simon_batch/newmatrices/1000pops/ Colombia2 500 &


sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Ecuador_mig 19 ~/project/tracts/sims/results/simon_batch/newmatrices/1000pops/ Ecuador 1000 &

sleep 2

./sim_manypops.sh ~/project/tracts/sims/data/simon_batch/Peru_mig 117 ~/project/tracts/sims/results/simon_batch/newmatrices/1000pops/ Peru 1000 &
