# OPP_lab4_mpi_Jakobi
Compile command:

    mpicxx main.c -o main -Wpedantic -Werror -Wall -O3 --std=c++11

MPE compile command:

    mpecxx -mpilog  main.c -o main -Wpedantic -Werror -Wall -O3 --std=c++11

Run command:

    mpirun -oversubscribe -np $proc_count ./main
