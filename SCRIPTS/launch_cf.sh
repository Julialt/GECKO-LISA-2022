#!/bin/bash
# Sidd Ghosh Feb 22, 2017
# for command file jobs
if [ ! -z "$PMI_RANK" ] ; then
   line=$(expr $PMI_RANK + 1)
elif [ ! -z "$OMPI_COMM_WORLD_RANK" ] ; then
   line=$(expr $OMPI_COMM_WORLD_RANK + 1)
else
   echo "Unsupported MPI exiting..."
   exit -1
fi
INSTANCE=$(sed -n ${line}p $1)
echo $INSTANCE
eval "$INSTANCE"
sleep 1
