#!/bin/bash
case $# in
	0)
		echo "Usage: bigrunner_ss.sh iproc0 // runs from idefault*1000 to idefault*1000 +999";
  	exit 1 ;;
	1)
		iproc0=$1
		nproc=24
		iprocf=`expr ${iproc0} + ${nproc}`
		for((i=iproc0;i<iprocf;i++))
		do
			mkdir -p smash_output/AuAu19.6/run${i}
			ln -s -f /home/scott/git/best_afterburner/external_codes/best_sampler/software/resinfo/pdg-SMASH.dat smash_output/AuAu19.6/run${i}/
			./sampler_and_smash -c smash_output/AuAu19.6/run${i}/config.yaml > logfiles/ss_${i}.txt &
		done
		
esac
