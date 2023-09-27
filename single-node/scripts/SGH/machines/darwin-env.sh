### Load environment modules here
### Assign names as relevant

mygcc="gcc/9.4.0"
myclang="clang/13.0.0"
mycuda="cuda/11.4.0"
myrocm="rocm"

if [ "$1" = "hpc" ]
then
    module purge
    if [ "$2" = "cuda" ]
    then
        module purge
        module load ${mygcc}
        module load ${mycuda}
    elif [ "$2" = "hip" ]
    then
        module purge
        module load ${myclang}
        module load ${myrocm}
    else
        module load ${mygcc}
    fi
    module load cmake
    module -t list
fi
