for name in calc_* ; do
  cp jobfile $name/job_${name/calc_/}
  cd $name
  qsub -l num_proc=1,s_rt=24:00:00,s_vmem=2G,h_fsize=2G -pe mpi 16 job_${name/calc_/}
  cd ..
done
