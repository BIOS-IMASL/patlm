

all_files="parametros.f95 proteina.f95 bulkimeter.f95 llamador.f95 printout.f95 adivinador.f95 fkfun.f95 freen.f95 mem_pore-pep.f95"



##### ##### #####
kinsol_lib="/usr/lib" #edit the path to your kinsol library


for i in libsundials_fkinsol libsundials_fnvecserial libsundials_kinsol libsundials_nvecserial;do
kinsol=$kinsol$kinsol_lib/$i.a" "
done

##### ##### #####

#flags="-g -fcheck=all -fbacktrace"
flags="-fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal"

#flags="-O3" #prodution flags
#flags="-g -fcheck=all -fbacktrace -ffpe-trap=overflow"


gfortran $flags $all_files $kinsol -o mem-pore-pep.x


