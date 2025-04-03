import random
import itertools
from sympy import symbols
from RSmodule import Reaction,RS
from Bool_module import BoolFunction
import copy
import numpy as np
from math import floor
from collections import defaultdict


def random_P_outsymb(RI_symb,out_symb,min_n_symb=3,max_symb=3):
    #nsymb = len(RI_symb)
    k = random.randint(min_n_symb,max_symb) # at least min_n_symb symbols and at most max_symb symbols
    RI_symb_list = list(RI_symb)
    np.random.shuffle(RI_symb_list)
    RI_symb = set(RI_symb_list[:k])
    # Reactants
    R = {x for x in RI_symb if random.choice((True, False))}
    # Inhibitors
    I = RI_symb.difference(R)
    # Products
    P = set()
    while P == set(): #to avoid empty products
        P = {x for x in out_symb if random.choice((True, False))}
    return Reaction(R,I,P)

class EvoRS:
    def __init__(self,pop_size, pop_elite, n_symb, init_size_min, init_size_max, exe_length, pc, pin, prm, pren, fit, BoolSymbols, max_mutation_tests):
        self.BoolSymb = BoolSymbols
        new_symb = n_symb - self.BoolSymb.n_symbols()
        if new_symb >= 0:
            t_symb = set(symbols(f't:{new_symb}'))
        else:
            t_symb = set()
            print("Number of symbols must be as many as the input symbols plus the output symbols")
            print(f"Setting n_symb = {self.BoolSymb.n_symbols()}\n")
        
        self.non_out_symb = self.BoolSymb.input_symb | t_symb #non-output symbols
        self.out_symb = self.BoolSymb.output_symb
        self.S = self.out_symb | self.non_out_symb # all symbols
        self.n_symb = len(self.S)
        
        #self.init_size = init_size
        self.init_size_min = init_size_min
        self.init_size_max = init_size_max

        self.exe_length = exe_length
        
        self.pc, self.pin, self.prm, self.pren = pc, pin, prm, pren
        self.max_tests = max_mutation_tests

        self.pop_size = pop_size
        self.pop = []
        self.popelite_size = pop_elite
        self.popelite = []
        
        if self.exe_length == 1:
            k = len(self.BoolSymb.input_tuple)
            #lambda *args: Reaction.random_P_outsymb(*args,min_n_symb=3)
            self.R_random = Reaction.random_P_outsymb   # function to generate random reaction 
            self.RS_random = RS.random_P_outsymb        # function to generate random RS
            if k%2:
                self.min_symb_per_reaction = k//2
            else:
                self.min_symb_per_reaction = k//2 + 1
        else:
            self.R_random = Reaction.random   # function to generate random reaction
            self.RS_random = RS.random        # function to generate random RS

        for _ in range(pop_size):
            init_size = random.randint(self.init_size_min,self.init_size_max)
            self.pop.append(self.RS_random(self.non_out_symb,self.out_symb,init_size,min_n_symb=self.min_symb_per_reaction))
        
        self.fit_func = fit
        
    def mutuation_insert_remove_2_nobal(self):
        # a number of reaction between min_rem and max_rem substituted
        min_rem = 1

        pop_size = len(self.pop)
        for i in range(pop_size):
            if random.random() < self.prm:
                rs = copy.deepcopy(self.pop[i])
                max_rem = floor(rs.n_reactions()*0.33)
                n_substitute = random.randint(min_rem, max_rem)
                to_remove = random.choices(self.pop[i].reactions,k = n_substitute)
                
                sizes = defaultdict(int)
                for r in to_remove:
                    sizes[len(r)] += 1 
                
                rs.remove_reactions(to_remove)
                
                new_reactions = []
                for size,n_new in sizes.items():
                    new_reactions += [random_P_outsymb(self.non_out_symb,self.out_symb,min_n_symb=size,max_symb=size) for _ in range(n_new)]
                rs.add_reactions(new_reactions)
                
                self.pop += [rs]
    
    def crossover_rs_same_len_r(self,rs1,rs2):
        sons = []
        sons_balanced = []
        reactions1_len = defaultdict(list)
        for r in rs1.reactions:
            reactions1_len[len(r)] += [r]   

        reactions2_len = defaultdict(list)
        for r in rs2.reactions:
            reactions2_len[len(r)] += [r]

        len_intersection = list(reactions1_len.keys() & reactions2_len.keys())
        #randomizzare qua per k 
        for k in len_intersection:
            k = random.choice(list(reactions1_len.keys() & reactions2_len.keys()))

            son1_reactions = list(itertools.chain(*[reactions1_len[x] for x in reactions1_len.keys() if x != k])) + reactions2_len[k]
            son2_reactions = list(itertools.chain(*[reactions2_len[x] for x in reactions2_len.keys() if x != k])) + reactions1_len[k]

            son1 = RS(son1_reactions,self.S)
            son2 = RS(son2_reactions,self.S)
            sons += [son1,son2]
            #for son in [son1,son2]:
            #    if self.bool_res(son).unbalancedness() == 0:
            #        sons_balanced += [son]
                    
        return sons
    
    def crossover_balanced(self):
        pop_size = len(self.pop)
        sons = []
        for i in range(1, pop_size):
            if random.random() < self.pc:
                #tournament selection (4 scelgo 2)
                gen = list(np.random.choice(self.pop,size=4,replace=False))
                list.sort(gen, key=self.fit,reverse=True)
                gen1,gen2 = gen[:2]
                if gen1.n_reactions() == 0 and gen2.n_reactions() == 0:
                    break
                sons_balanced,_ = self.crossover_rs_same_len_r(gen1, gen2)
                sons += sons_balanced
        self.pop += sons
    
    def crossover_nobalanced(self):
        pop_size = len(self.pop)
        sons = []
        for i in range(1, pop_size):
            if random.random() < self.pc:
                #tournament selection (4 scelgo 2)
                gen = list(np.random.choice(self.pop,size=4,replace=False))
                list.sort(gen, key=self.fit,reverse=True)
                gen1,gen2 = gen[:2]
                if gen1.n_reactions() == 0 and gen2.n_reactions() == 0:
                    break
                new_sons = self.crossover_rs_same_len_r(gen1, gen2)
                sons += new_sons
        self.pop += sons
    
    def mutuation_renewal(self):
        pop_size = len(self.pop)
        new_mutated = []
        for i in range(pop_size):
            if random.random() < self.pren:
                init_size = random.randint(self.init_size_min,self.init_size_max)
                new_mutated += [self.RS_random(self.non_out_symb,self.out_symb,init_size,min_n_symb=self.min_symb_per_reaction)]
        self.pop += new_mutated
    
    def minimization(self):
        for i,_ in enumerate(self.pop):
            self.pop[i].minimize()
                    
    def execute_rs_le(self,rs,W):
        # f(W) is True if and only if True \in res^{i}(W) for i=1,...,exe_length
        for _ in range(self.exe_length):
            W = rs.res(W)
            if self.BoolSymb.output_symb & W != set():
                return True
        return False
    
    def execute_rs_eq(self,rs,W):
        # f(W) is True if and only if True \in res^{exe_length}(W)
        for _ in range(self.exe_length):
            W = rs.res(W)
        
        if self.BoolSymb.output_symb & W != set():
            return True
        return False
    
    def bool_res(self,rs):
        def bool_res_conversion(*argv):
            input_tuple = self.BoolSymb.input_tuple
            W = set([x for x in input_tuple if argv[input_tuple.index(x)] == 1])
            return self.execute_rs_le(rs,W) # self.execute_rs_eq(rs,W) # 
        n_input = len(self.BoolSymb.input_tuple)
        return BoolFunction(bool_res_conversion,n_input)
    
    def fit(self,rs):
        return self.fit_func(self.bool_res(rs))
    
    def truncated_selection(self,flag):
        self.pop += self.popelite

        self.pop,_ = unique_duplicate_list(self.pop)

        # discard rs that rappresent the same bool function
        #self.select_unique_func()

        list.sort(self.pop, key=self.fit,reverse=True)
        if len(self.pop) > self.pop_size:
            self.pop = self.pop[:self.pop_size]
        
        self.n_unique_elements =  copy.deepcopy(len(self.pop))
        
        self.bestRS = copy.deepcopy(self.pop[0])
        self.worstRS = copy.deepcopy(self.pop[-1])
        mid = self.pop_size // 2
        self.medianRS = copy.deepcopy(self.pop[mid])

        self.popelite = copy.deepcopy(self.pop[:self.popelite_size])
        
        if flag:
            k=self.popelite_size + 5
            print(f"  numbers of unique RS: {len(self.pop)}")
            print(f"  fitness first {k} elements: {[self.fit(self.pop[j]) for j in range(k)]}")
            print(f"  fitness max,median,min: {self.fit(self.bestRS)}, {self.fit(self.medianRS)}, {self.fit(self.worstRS)}")

        np.random.shuffle(self.pop)
    
    def ordering(self):
        list.sort(self.pop, key=self.fit,reverse=True)

    def evolution_cycle(self,flag_print=False):
        # Selection
        self.truncated_selection(flag=flag_print)
        # Renwal
        self.mutuation_renewal()
        # Crossover
        #self.crossover()
        # Crossover nobalanced
        self.crossover_nobalanced
        # Random removal
        #self.mutuation_removal()
        # Random insertion
        #self.mutuation_insetion()
        # Insert-remove 
        #self.mutuation_insert_remove()
        self.mutuation_insert_remove_2_nobal()
        # Minimization
        self.minimization()       
    
    def __str__(self):
        s = ""
        for i,rs in enumerate(self.pop):
            s += f"reaction system {i} (fittness {self.fit(rs)}): \n    {rs}"
        return s
    
    def select_unique_func(self):
        uniqueList = []
        uniqueList_func = []
        for rs in self.pop:
            f = self.bool_res(rs)
            if f not in uniqueList_func:
                uniqueList_func.append(f)
                uniqueList.append(rs)
        self.pop = uniqueList

def unique_duplicate_list(lista):
    uniqueList = []
    duplicateList = []
    for i in lista:
        if i not in uniqueList:
            uniqueList.append(i)
        elif i not in duplicateList:
            duplicateList.append(i)
    return uniqueList, duplicateList

