import sys, os
sys.path.append(os.path.realpath('..'))
import cPickle
import tracts_ped as ped

popdir = os.path.expanduser(sys.argv[1])
outdir = os.path.expanduser(sys.argv[2])

popfiles = [os.path.join(popdir, popfile) for popfile in os.listdir(popdir)
			if os.path.splitext(popfile)[1] == '.pkl']

for popfile in popfiles:
	with open(popfile, 'r') as f:
		pop = cPickle.load(f)
		for i in range(len(pop.indivs)):
			indfile = "IND" + str(i)
			inddir = os.path.splitext(popfile)[0].strip('_pop')
			if not os.path.exists(os.path.join(outdir, inddir)):
				os.makedirs(os.path.join(outdir, inddir))
			outfile = os.path.join(outdir, inddir, indfile)
			print "Outputting to", outfile
			ped.tracts_ind_to_bed(pop.indivs[i], outfile)