import random
import itertools
from sympy import symbols
import numpy as np
from collections import defaultdict
#import graphviz

from colored import fg, attr

def print_set(T):
    if T == set():
        return "\u2205"
    T_print = list(T)
    s = "{"
    reset = attr('reset')
    for a in T_print[:-1]:
        if symbols('True') == a:
            s+=fg('blue')+"True"+reset+", "  
        elif a in symbols(f't:{100}'):
            s+=fg('green')+f"{a}"+reset+", "
        else:
            s += f"{a}, "

    if symbols('True') == T_print[-1]:
        s+=fg('blue')+"True"+reset+"}"
    elif T_print[-1] in symbols(f't:{100}'):
         s+=fg('green')+f"{T_print[-1]}"+reset+"}"
    else:
        s += f"{T_print[-1]}"+"}"
    return s

def HTMLset(T):
    if T == set():
        return "<\u2205>"
    T_print = list(T)
    s = "<{"
    for a in T_print[:-1]:
        if symbols('True') == a:
            s+='<font color="blue">'+str(a)+'</font>'+", "
        elif a in symbols(f't:{100}'):
            s+='<font color="green">'+str(a)+'</font>'+", "
        else:
            s += str(a) +", "

    if symbols('True') == T_print[-1]:
        s+='<font color="blue">'+str(T_print[-1])+'</font>'
    elif T_print[-1] in symbols(f't:{100}'):
         s+='<font color="green">'+str(T_print[-1])+'</font>'
    else:
        s += str(T_print[-1])
    s += "}>"
    return s


def powerset(iterable):
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))


class Reaction:
    def __init__(self,R,I,P):
        self.R = frozenset(R)
        self.I = frozenset(I)
        self.P = frozenset(P)
        self.size = len(self.R | self.I)
    
    def is_enabled(self,T):
        if T >= self.R and T & self.I == set():
            return True
        return False
    
    def res(self,T):
        if self.is_enabled(T):
            return self.P
        return set()
    
    def __eq__(self, other):
        if self.R == other.R and self.I == other.I and self.P == other.P:
            return True
        return False
    
    def __le__(self, other):
        if self.R >= other.R and self.I >= other.I and self.P <= other.P:
            return True
        return False
    
    def __str__(self):
         return f"({print_set(self.R)}, {print_set(self.I)}, {print_set(self.P)})"
    
    def __hash__(self):
        # prova
        return (hash(self.R) ^ hash(self.I) ^ hash(self.P))
    
    def __len__(self) -> int:
        """
        Returns the number of reactants + the number of inhibitors
        """
        return len(self.R | self.I)

    def ordering(self):
        # order the sets in lexicographic order so that the hash becomes unique (?)
        pass #todo
    
    @staticmethod
    def random(RI_symb,other_symbols):
        # Reactants
        #R = set()
        R = {x for x in RI_symb if random.choice((True, False))}
        # Inhibitors
        #I = set()
        I = {x for x in RI_symb.difference(R) if random.choice((True, False))}
        # Products
        S_symb = other_symbols.union(RI_symb)
        P = set()
        while P == set(): #to avoid empty products
            P = {x for x in S_symb if random.choice((True, False))}
        return Reaction(R,I,P)
    
    @staticmethod
    def random_maximal(RI_symb,out_symb):
        # Reactants
        R = {x for x in RI_symb if random.choice((True, False))}
        # Inhibitors
        I = RI_symb.difference(R)
        # Products
        P = set()
        while P == set(): #to avoid empty products
            P = {x for x in out_symb if random.choice((True, False))}
        return Reaction(R,I,P)
    
    @staticmethod
    def random_P_outsymb(RI_symb,out_symb,min_n_symb=3):
        # fare in modo che con k pari ci siano almeno k//2 elementi in totale tra reagenti e inibitori
        nsymb = len(RI_symb)
        k = random.randint(min_n_symb,nsymb) # at least min_n_symb symbols
        RI_symb_list = list(RI_symb)
        np.random.shuffle(RI_symb_list)
        RI_symb = set(RI_symb_list[:k])
        # Reactants
        R = {x for x in RI_symb if random.choice((True, False))}
        # Inhibitors
        I = RI_symb.difference(R) #{x for x in RI_symb.difference(R) if random.choice((True, False))}
        # Products
        P = set()
        while P == set(): #to avoid empty products
            P = {x for x in out_symb if random.choice((True, False))}
        return Reaction(R,I,P)
    
    
class RS:
    def __init__(self,reactions,S):
        self.reactions = reactions
        self.S = frozenset(S)

    def ordering(self):
        list.sort(self.reactions, key=len)
        
    def n_reactions(self):
        return len(self.reactions)
    
    def remove_reaction(self,reaction):
        self.reactions.remove(reaction) # se ci sono più copie di reaction ne tira via solo la prima che incontra
    
    def remove_reactions(self,reactions):
        self.reactions = [r for r in self.reactions if r not in reactions]
    
    def add_reactions(self,reactions):
         self.reactions += reactions
    
    def minimize(self):
        rm_reactions = []
        # minimization
        for r1,r2 in itertools.permutations(self.reactions,2): 
            if r1<=r2:
                rm_reactions.append(r1)
        self.remove_reactions(rm_reactions)

    def sizes_reactions(self):
        sizes = set([len(r) for r in self.reactions])
        sizes_dict = {size: 0 for size in sizes}
        for r in self.reactions:
            sizes_dict[len(r)] += 1
        return sizes_dict
    
    def res(self,T):
        return set().union(*[r.res(T) for r in self.reactions])
    
    @staticmethod
    def random(RI_symb,other_symbols,max_reactions):
        n_reactions = random.randint(1, max_reactions)
        reactions = [Reaction.random(RI_symb,other_symbols) for _ in range(n_reactions)]
        S_symb = RI_symb.union(other_symbols)
        return RS(reactions, S_symb)     

    @staticmethod
    def random_maximal(RI_symb,out_symb,n_reactions):
        reactions = [Reaction.random_maximal(RI_symb,out_symb) for _ in range(n_reactions)]
        S_symb = set().union(*[RI_symb,out_symb])
        return RS(reactions, S_symb)

    @staticmethod
    def random_P_outsymb(RI_symb,out_symb,max_reactions,min_n_symb=3):
        #n_reactions = random.randint(1, max_reactions)
        n_reactions = max_reactions
        reactions = [Reaction.random_P_outsymb(RI_symb,out_symb,min_n_symb) for _ in range(n_reactions)]
        S_symb = set().union(*[RI_symb,out_symb])
        return RS(reactions, S_symb)   

    def graph(self):
        S = tuple(self.S)
        sort_S={x: i for i,x in enumerate(S)}

        G = graphviz.Digraph()
        G.attr('node',style='filled')
        for T in powerset(S):
            if symbols('True') in T:
                fillcolor = 'lightblue'
            else:
                fillcolor = 'white'
            G.node(name= str(T),
                       label=HTMLset(set(iter(T))),
                       color='lightgray',
                       fillcolor= 'white',
                       fontsize="12")

        for T_in in powerset(S):
            T_out = self.res(set(x for x in T_in))
            if T_out != set():
                T_out = tuple(sorted(list(x for x in T_out),key=lambda x: sort_S[x]))
            else:
                T_out = ()
            G.edge(str(T_in),str(T_out),color = 'gray')
        return G

    def __str__(self):
        self.ordering()
        s = f"number of reactions: {self.n_reactions()}\n"
        s += f"reactions sizes: {self.sizes_reactions()}\n"
        for r in self.reactions:
            s += f"\t{r}\n"
        return s
    
    def __eq__(self,other):
        if self.S == other.S:
            if set(self.reactions) == set(other.reactions):
                return True
        return False


