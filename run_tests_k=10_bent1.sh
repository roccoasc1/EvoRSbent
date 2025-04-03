#!/bin/bash

# Initialize Conda
eval "$(conda shell.bash hook)"
conda activate evoRS

# Array of models and dimensions
k_range=(10)
pop=100
era=2000
popelite=5
exe_length=1
init_size_min1=80
init_size_max1=160
fit_range=(5)
pc=0.8
prm=0.2
pin=0.2
pren=0.1
max_tests=10

run=(0 1 2 3 4)
# Run python scripts in parallel
for k in "${k_range[@]}"; do
    #echo "Running python for k $k"
    for r in "${run[@]}"; do
        for fit in "${fit_range[@]}"; do
            python EvoRS_run_fitnl_new_mutation_nobal.py --k=$k --pop_size=$pop --init_size_min=$init_size_min1 --init_size_max=$init_size_max1 --n_cycles=$era --popelite_size=$popelite --exe_length=$exe_length --run=$r --pc=$pc --pin=$pin --prm=$prm --pren=$pren --fit=$fit --max_tests=$max_tests&
        done
    done
done

# Wait for all background processes to finish
wait


run=(5 6 7 8 9)

# Run python scripts in parallel
for k in "${k_range[@]}"; do
    #echo "Running python for k $k"
    for r in "${run[@]}"; do
        for fit in "${fit_range[@]}"; do
            python EvoRS_run_fitnl_new_mutation_nobal.py --k=$k --pop_size=$pop --init_size_min=$init_size_min1 --init_size_max=$init_size_max1 --n_cycles=$era --popelite_size=$popelite --exe_length=$exe_length --run=$r --pc=$pc --pin=$pin --prm=$prm --pren=$pren --fit=$fit --max_tests=$max_tests&
        done
    done
done

# Wait for all background processes to finish
wait

run=(10 11 12 13 14)
# Run python scripts in parallel
for k in "${k_range[@]}"; do
    #echo "Running python for k $k"
    for r in "${run[@]}"; do
        for fit in "${fit_range[@]}"; do
            python EvoRS_run_fitnl_new_mutation_nobal.py --k=$k --pop_size=$pop --init_size_min=$init_size_min1 --init_size_max=$init_size_max1 --n_cycles=$era --popelite_size=$popelite --exe_length=$exe_length --run=$r --pc=$pc --pin=$pin --prm=$prm --pren=$pren --fit=$fit --max_tests=$max_tests&
        done
    done
done

# Wait for all background processes to finish
wait


run=(15 16 17 18 19)

# Run python scripts in parallel
for k in "${k_range[@]}"; do
    #echo "Running python for k $k"
    for r in "${run[@]}"; do
        for fit in "${fit_range[@]}"; do
            python EvoRS_run_fitnl_new_mutation_nobal.py --k=$k --pop_size=$pop --init_size_min=$init_size_min1 --init_size_max=$init_size_max1 --n_cycles=$era --popelite_size=$popelite --exe_length=$exe_length --run=$r --pc=$pc --pin=$pin --prm=$prm --pren=$pren --fit=$fit --max_tests=$max_tests&
        done
    done
done

# Wait for all background processes to finish
wait

run=(20 21 22 23 24)
# Run python scripts in parallel
for k in "${k_range[@]}"; do
    #echo "Running python for k $k"
    for r in "${run[@]}"; do
        for fit in "${fit_range[@]}"; do
            python EvoRS_run_fitnl_new_mutation_nobal.py --k=$k --pop_size=$pop --init_size_min=$init_size_min1 --init_size_max=$init_size_max1 --n_cycles=$era --popelite_size=$popelite --exe_length=$exe_length --run=$r --pc=$pc --pin=$pin --prm=$prm --pren=$pren --fit=$fit --max_tests=$max_tests&
        done
    done
done

# Wait for all background processes to finish
wait


run=(25 26 27 28 29)

# Run python scripts in parallel
for k in "${k_range[@]}"; do
    #echo "Running python for k $k"
    for r in "${run[@]}"; do
        for fit in "${fit_range[@]}"; do
            python EvoRS_run_fitnl_new_mutation_nobal.py --k=$k --pop_size=$pop --init_size_min=$init_size_min1 --init_size_max=$init_size_max1 --n_cycles=$era --popelite_size=$popelite --exe_length=$exe_length --run=$r --pc=$pc --pin=$pin --prm=$prm --pren=$pren --fit=$fit --max_tests=$max_tests&
        done
    done
done

# Wait for all background processes to finish
wait


echo "All tasks have been launched."