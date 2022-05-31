#/bin/bash

mkdir -p ../GENERATED_SCRIPTS
mkdir -p /glade/scratch/$USER/stdout

#=============================================================
#==  Functions to write a script to be submitted to cheyenne =
#=============================================================
function write_cheyenne_script() {
    generated_script=$1
    jobname=$2
    walltime=$3
    nthreads=$4
    outputfile=$5
    errorfile=$6
    called_script=$7
    # if we don't need all arguments, just don't provide them
    arg1=$8
    arg2=$9
    arg3=${10}
    arg4=${11}
    arg5=${12}
    

    
if [ ${nthreads} -eq 1 ]; then
  cpurequest=select=1:ncpus=${nthreads}
  queue=share
elif [ ${nthreads} -gt 1 ]; then
  cpurequest=select=1:ncpus=${nthreads}:ompthreads=${nthreads}
  queue=economy
else
  echo wrong number of threads requested: ${nthreads}
  exit 1
fi

    echo "#!/bin/bash
#PBS -N $jobname
#PBS -A P19010000
#PBS -l walltime=$walltime
#PBS -l ${cpurequest}
#PBS -q ${queue}
#PBS -j oe
#PBS -o $scratch_dir/LOGS/$outputfile
#PBS -m abe
#PBS -l pmem=100G

export TMPDIR=/glade/scratch/$USER/stdout

module unload netcdf
module unload intel
module load gnu
module load netcdf

$called_script \"$arg1\" \"$arg2\" \"$arg3\" \"$arg4\" \"$arg5\"

    " > $home_dir/GENERATED_SCRIPTS/${generated_script}

    return
}

#==========================================================
#==  Functions to write a script for multiple simulations on same job==
#==========================================================
function write_cheyenne_multipsim_script() {
    generated_script=$1
    jobname=$2
    walltime=$3
    nthreads=$4
    outputfile=$5
    errorfile=$6
    called_script=$7
    # if we don't need all arguments, just don't provide them
    arg1=$8
    arg2=$9
    arg3=${10}
    arg4=${11}
    

    
if [ ${nthreads} -eq 1 ]; then
  cpurequest=select=1:ncpus=${nthreads}
  queue=share
elif [ ${nthreads} -gt 1 ]; then
  cpurequest=select=1:ncpus=${nthreads}:mpiprocs=${nthreads}
  queue=economy
else
  echo wrong number of threads requested: ${nthreads}
  exit 1
fi
# we want 300 cores
# ie 10 nodes with 30 cores is the closest
#cpurequest=select=10:ncpus=30:mpiprocs=30

    echo "#!/bin/bash
#PBS -N $jobname
#PBS -A P19010000
#PBS -l walltime=$walltime
#PBS -l ${cpurequest}
#PBS -q ${queue}
#PBS -j oe
#PBS -o $scratch_dir/LOGS/$outputfile
#PBS -m abe
#PBS -l pmem=5G

export TMPDIR=/glade/scratch/$USER/stdout

module unload netcdf
module unload intel
module load gnu
module load netcdf
export MPI_SHEPHERD=true

mpiexec_mpt ./launch_cf.sh $called_script

    " > $home_dir/GENERATED_SCRIPTS/${generated_script}

    return
}


function write_cheyenne_dependentscript() {
    generated_script=$1
    jobname=$2
    walltime=$3
    nthreads=$4
    outputfile=$5
    errorfile=$6
    dependency_expression=$7
    called_script=$8
    # if we don't need all arguments, just don't provide them
    arg1=$9
    arg2=${10}
    arg3=${11}
    arg4=${12}
    arg5=${13}
    
if [ ${nthreads} -eq 1 ]; then
  cpurequest=select=1:ncpus=${nthreads}
  queue=share
elif [ ${nthreads} -gt 1 ]; then
  cpurequest=select=1:ncpus=${nthreads}:ompthreads=${nthreads}
  queue=economy
else
  echo wrong number of threads requested: ${nthreads}
  exit 1
fi


    echo "#!/bin/bash

#PBS -N $jobname
#PBS -A P19010000
#PBS -l walltime=$walltime
#PBS -l ${cpurequest}
#PBS -q ${queue}
#PBS -j oe
#PBS -o $scratch_dir/LOGS/$outputfile
#PBS -m abe
#PBS -W depend=$dependency_expression
#PBS -l pmem=100G

export TMPDIR=/glade/scratch/$USER/stdout

module unload netcdf
module unload intel
module load gnu
module load netcdf

$called_script \"$arg1\" \"$arg2\" \"$arg3\" \"$arg4\" \"$arg5\"

    " > ${home_dir}/GENERATED_SCRIPTS/${generated_script}

    return
}

function write_cheyenne_monoproc_script() {
    generated_script=$1
    jobname=$2
    walltime=$3
    outputfile=$4
    errorfile=$5
    called_script=$6
    # if we don't need all arguments, just don't provide them
    arg1=$7
    arg2=$8
    arg3=$9
    arg4=${10}

    echo "#!/bin/bash

#PBS -N $jobname
#PBS -A P19010000
#PBS -l walltime=$walltime
#PBS -l select=1:ncpus=1
#PBS -q economy 
#PBS -j oe
#PBS -o $scratch_dir/LOGS/$outputfile
#PBS -m abe
#PBS -l pmem=100G

export TMPDIR=/glade/scratch/$USER/stdout

module unload netcdf
module unload intel
module load gnu
module load netcdf

$called_script \"$arg1\" \"$arg2\" \"$arg3\" \"$arg4\"

    " > ${home_dir}/GENERATED_SCRIPTS/${generated_script}

    return
}

function write_cheyenne_monoproc_dependentscript() {
    generated_script=$1
    jobname=$2
    walltime=$3
    outputfile=$4
    errorfile=$5
    dependency_expression=$6
    called_script=$7
    # if we don't need all arguments, just don't provide them
    arg1=$8
    arg2=$9
    arg3=${10}
    arg4=${11}
    arg5=${12}
    arg6=${13}

    echo "#!/bin/bash

#PBS -N $jobname
#PBS -A P19010000
#PBS -l walltime=$walltime
#PBS -l select=1:ncpus=1
#PBS -q economy 
#PBS -j oe
#PBS -o $scratch_dir/LOGS/$outputfile
#PBS -m abe
#PBS -W depend=$dependency_expression
#PBS -l pmem=100G

export TMPDIR=/glade/scratch/$USER/stdout

module unload netcdf
module unload intel
module load gnu
module load netcdf

$called_script \"$arg1\" \"$arg2\" \"$arg3\" \"$arg4\" \"$arg5\"\"$arg6\"

    " > ${home_dir}/GENERATED_SCRIPTS/${generated_script}

    return
}
