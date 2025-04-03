import numpy as np
import EvoRSmodule_newmuta_nobal as my #same module as EvoRSmodule_newmuta with unbalanceness constrain
import test_func
import matplotlib.pyplot as plt
import pickle
import threading
from collections import defaultdict
from sympy import symbols
from Bool_module import BoolSymbols
from fitness import fit_nl, fit_nl_unbal_normalize, fit_kroneker, fit_max_values, fit_2, fit_nl_max_values

import logging
import threading
import time
import argparse

import argparse
import random

import os
import csv
    
def create_cd(cd):
    if not os.path.exists(cd):
        os.makedirs(cd)
        print(f"new directory created: {cd}")

k_bound_bent = {6: 28, 8: 120, 10: 496}

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", default=8,type=int)
    parser.add_argument("--pop_size", default=50,type=int)
    parser.add_argument("--n_cycles", default=2000,type=int)
    parser.add_argument("--popelite_size", default=5,type=int)
    parser.add_argument("--init_size_min", default=40,type=int)
    parser.add_argument("--init_size_max", default=80,type=int)
    parser.add_argument("--exe_length", default=1,type=int)
    parser.add_argument("--pc", default=0.8,type=float)
    parser.add_argument("--pin", default=0.2,type=float)
    parser.add_argument("--prm", default=0.2,type=float)
    parser.add_argument("--pren", default=0.1,type=float)
    parser.add_argument("--run", default=0,type=int)
    parser.add_argument("--fit", default=5,type=int)
    parser.add_argument("--max_tests", default=1,type=int)

    args = parser.parse_args()
    
    # Setting parameters
    k = args.k # number of boolean inputs
    pop_size = args.pop_size
    n_cycles = args.n_cycles
    popelite_size = args.popelite_size
    exe_length = args.exe_length
    init_size_min = args.init_size_min
    init_size_max = args.init_size_max
    run = args.run
    max_tests = args.max_tests
    
    if args.fit == 0:
        fit = fit_nl
        cd_fit = f'2024.12.3_fit_nl'
    elif args.fit == 1:
        fit = fit_nl_unbal_normalize
        cd_fit = f'2024.12.3_fit_nl_unbal'
    elif args.fit == 2:
        fit = fit_kroneker
        cd_fit = f'2024.12.12_fit_nl_kroneker'
    elif args.fit == 3:
        fit = fit_max_values
        cd_fit = f'2024.12.20_fit_max_values'
    elif args.fit == 4:
        fit = fit_2
        cd_fit = f'2025.01.09_fit_2_bent'
    elif args.fit == 5:
        fit = fit_nl
        cd_fit = f'2025.01.09_fit_nl_bent'
    elif args.fit == 8:
        fit = fit_nl_max_values
        cd_fit = f'2025.01.17_fit_nl_bent_max_values'
    

    create_cd(f'./{cd_fit}')

    for cd in [f'./{cd_fit}/all_fitness/', f'./{cd_fit}/best_RS/', f'./{cd_fit}/best_RS_properties/']:
        create_cd(cd)

    pc, pin, prm, pren = args.pc,args.pin,args.prm,args.pren
    
    print(f"parameters: \n {k=} \n {init_size_min=} \n {init_size_max=} \n {pop_size=} \n {n_cycles=} \n {popelite_size=} \n {exe_length=} \n {pc=}, {pin=}, {prm=}, { pren=}")
    
    def get_file_name(ccycle):
        return f'{k=}_run{run}_currentcycle{ccycle}_maxevcycle{n_cycles}_{init_size_min=}_{init_size_max=}_popsize{pop_size}_popelite{popelite_size}_exelength{exe_length}_{pc=}_{pin=}_{prm=}_{pren=}'
    
    def get_file_name_final():
        return f'{k=}_run{run}_evcycle{n_cycles}_{init_size_min=}_{init_size_max=}_popsize{pop_size}_popelite{popelite_size}_exelength{exe_length}_{pc=}_{pin=}_{prm=}_{pren=}'
        
    random.seed(run)
    
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO,datefmt="%H:%M:%S") 
    logging.info(f"Thread {k=}, {run=}: starting")
    
    # Inizialization
    columns = ['fit_max', 'fit_median', 'bestRS_nl','bestRS_unbal','unique_pop']
    rows = []
    
    #create csv file
    cyclesaved=0
    #file_name = f'run{i}_currentcycle{cyclesaved}_maxevcycle{n_cycles}_popsize{pop_size}_{k=}_popelite{popelite_size}_exelength{exe_length}_{pc=}_{pin=}_{prm=}_{pren=}'
    file_name = get_file_name(cyclesaved)

    cd = f'./{cd_fit}/all_fitness/'
    with open(f'{cd}{file_name}.csv', "w",newline='') as f: 
        # Create a csv.writer object
        writer = csv.writer(f)
        # Write data to the CSV file
        writer.writerow(columns)
        
    
    #setting paremeters depending on k
    n_symb = (k+1)
    
    ERS = my.EvoRS(pop_size, popelite_size, n_symb, init_size_min, init_size_max, exe_length, pc, pin, prm, pren, fit, BoolSymbols(k), max_tests)

    # Evolution Cycle
    for j in range(n_cycles):
        #termination criteria
        #if ERS.fit(ERS.pop[0]) >= 0.99:
            #break
        ERS.evolution_cycle()

        bestf = ERS.bool_res(ERS.bestRS)
        
        rows += [[ERS.fit(ERS.bestRS),ERS.fit(ERS.medianRS),bestf.nonlinearity(),bestf.unbalancedness(),ERS.n_unique_elements]]

        maxnl = bestf.nonlinearity()

        if (j%25 == 0 and j > 0) or j == n_cycles-1 or maxnl == k_bound_bent[k]:
            logging.info(f"Thread {k=}, {run=}: current cycle {j}")
            file_name_old = get_file_name(cyclesaved)
            file_name_new = get_file_name(j)
        
            cd = f'./{cd_fit}/all_fitness/'
            with open(f'{cd}{file_name_old}.csv', "a",newline='') as f: 
                writer = csv.writer(f)
                writer.writerows(rows)
                rows = []
            os.replace(f'{cd}{file_name_old}.csv',f'{cd}{file_name_new}.csv')

            cdRS = f'./{cd_fit}/best_RS/'

            if j%100 == 0:
                with open(f'{cdRS}{file_name_new}.pkl', 'wb') as f:
                    pickle.dump(ERS.bestRS, f, protocol=pickle.HIGHEST_PROTOCOL)
            
            cyclesaved=j

        if j == n_cycles-1 or maxnl == k_bound_bent[k]:
            cd = f'./{cd_fit}/all_fitness/'            
            file_name_old = get_file_name(cyclesaved)
            file_name_final = get_file_name_final()
            
            os.replace(f'{cd}{file_name_old}.csv',f'{cd}{file_name_final}.csv')

        if maxnl == k_bound_bent[k]:
            break

    # rifaccio la truncated selection prima di salvare la final pop
    ERS.truncated_selection(flag=False)

    # c'è un problema quando salviamo la popolazione finale qua non abbiamo fatto la truncated selection 
    # quindi non abbiamo gli individui migliori, meglio aggiungere una truncated selection alla fine dell'ultimo ciclo        
    cdall = f'./{cd_fit}/final_pop/'
    file_name_final = get_file_name_final()
    create_cd(cdall)
    with open(f'{cdall}{file_name_new}.pkl', 'wb') as f:
        pickle.dump(ERS.pop, f, protocol=pickle.HIGHEST_PROTOCOL)


    cd = f'./{cd_fit}/best_RS/'
    file_name_final = get_file_name_final()

    with open(f'{cd}{file_name_final}.pkl', 'wb') as f:
        pickle.dump(ERS.bestRS, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    cd = f'./{cd_fit}/best_RS_properties/'
    file_best_RS_properties = f"{cd}run{run}_evcycle{n_cycles}_popsize{pop_size}_{k=}_popelite{popelite_size}_exelength{exe_length}_{pc=}_{pin=}_{prm=}_{pren=}.txt"
    
    with open(file_best_RS_properties, 'w') as f:
        bestf = ERS.bool_res(ERS.bestRS)
        print(f'nl(f): {bestf.nonlinearity()}',file= f)
        print(f'w_H(f): {bestf.hamming_weight()} , unbalenceness: {abs(2**(bestf.n_input-1)-bestf.hamming_weight())}',file= f)
        print(f'W_max(f): {bestf.spectral_radius_new()}',file= f)
    
    logging.info(f"Thread {k=}, {run=}: finished")  
    print(f"parameters: \n {k=} \n {init_size_min=} \n {init_size_max=} \n {pop_size=} \n {n_cycles=} \n {popelite_size=} \n {exe_length=} \n {pc=}, {pin=}, {prm=}, { pren=}")