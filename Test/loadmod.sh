part=${HOME:1:7}
if [[ "$part" == "public1" ]];then
    source /${part}/soft/modules/module.sh
    module purge
    module load mpi/intel/18.2
elif [[ "$part" == "public3" ]];then
    source /${part}/soft/modules/module.sh
    module purge
    module load mpi/intel/18.0.2-public3
else
    echo "local cluster"
fi
me=$(basename "$0")
if [[ "$me" == "loadmod.sh" ]];then
make realclean && make -j 48 
fi
