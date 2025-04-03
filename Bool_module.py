from prettytable import PrettyTable
from sympy import symbols
import numpy as np
import math

class BoolSymbols:
    def __init__(self, n_input_symbols,set_true=True):
        self.input_symb = set(symbols(f'x:{n_input_symbols}'))
        self.input_tuple = tuple(symbols(f'x:{n_input_symbols}'))
        if set_true:
            self.output_symb = set([symbols('True')])
        else:
            self.output_symb = set()
        
    def n_symbols(self):
        return len(self.input_symb)+len(self.output_symb)
    
    def test_input(self):
        k = len(self.input_symb)
        return [tuple(map(int, bin(j)[2:].zfill(k))) for j in range(2**k)]
    

'''
    * Computes the Moebius Transform of a Boolean function using the Fast
    * Moebius Transform (FMT) algorithm, which requires O(NlogN) operations
    * (N=2^n is the length of the truth table). The method directly computes
    * the spectrum in the original vector, and it must be called with the
    * initial parameters (vector, 0, vector.length). The return value of the
    * method is the algebraic degree of the function.
    * 
    * The reference for the algorithm is Carlet, "Cryptography and
    * Error-Correcting Codes", chapter 8 in Crama, Hammer, "Boolean Models and
    * Methods in Mathematics, Computer Science and Engineering", p. 263.
    * 
    * NOTICE: the Moebius transform is an involution, meaning that this method
    * can also be used to recover the truth table of a function starting from
    * its ANF vector.
    * 
    * @param vector    An array boolean representing the boolean function.
    * @param start     The index of the truth table where to start computations.
    * @param length    The length of the subvector where to perform computations
    *                  starting from start.
    * @return          The algebraic degree

    Credits: Luca Mariot GitHub (https://github.com/rymoah/CryptoBoolFun/tree/master)
'''
def computeFMT(vector, start, length):
    half = length//2
    
    #Main cycle: split vector in two parts (v0 e v1),
    #update v1 as v1=v0 XOR v1.
    for i in range(start,start+half):
        vector[i+half] = vector[i] ^ vector[i+half]
    
    #Recursive call on v0 and v1.
    if half>1:
        val1 = computeFMT(vector,start,half);
        val2 = computeFMT(vector,start+half,half);
        
        #At the end of the recursive calls, compare val1 and val2 to decide
        #what is the algebraic degree in the portion of ANF included
        #between start and start+length
        if val1 > val2:
            return val1    
        else:
            return val2
    else:
        #If we have reached half=1 (function of 2 variables),
        #return the Hamming weight of the largest input with nonzero
        #coefficient in the current subvector. This is the algebraic degree
        #of the restriction of the function of 2 variables.
            
        #If both coefficient are zero, then the degree of the subfunction
        #is zero.
        if (vector[start] == False) and (vector[start+half] == False): 
            return 0
        else:
            #Compute the length of the input vectors.
            inplen = int(math.log(vector.size)/math.log(2))

            #If the coefficient of the higher vector is null,
            #then return the Hamming weight of the lower vector.
            if vector[start+half] == False:
                input_tuple =  tuple(np.binary_repr(start, width=inplen))
                input = np.array(input_tuple, dtype=bool)
                subdeg = sum(input)

                return subdeg
            
            else:
                #In all other cases, return the Hamming weight of the
                #higher vector.
                input_tuple =  tuple(np.binary_repr(start+half, width=inplen))
                input = np.array(input_tuple, dtype=bool)
                input = np.array(input_tuple, dtype=bool)
                subdeg = sum(input)

                return subdeg



'''
     * Computes the Walsh Transform of a Boolean function using the Fast Walsh
     * Transform (FWT) algorithm, which requires O(NlogN) operations (N=2^n is
     * the length of the truth table). The method directly computes the spectrum
     * in the original vector, and it must be called with the initial parameters
     * (vector, 0, vector.length). The return value is the spectral radius of
     * the function (=maximum absolute value of the Walsh transform).
     * The reference for the algorithm is C. Carlet, "Cryptography and
     * Error-Correcting Codes", chapter 8 in Crama, Hammer, "Boolean Models and
     * Methods in Mathematics, Computer Science and Engineering", p. 272.
     * 
     * @param vector    An array of int representing the truth table of a
     *                  Boolean function in polar form. (polar form means -1 False and 1 True)
     * @param start     The index of the truth table where to start computations.
     * @param length    Length of the portion of truth table where to perform
     *                  computations, starting from start
     * @return          The spectral radius of the function, computed as the
     *                  maximum absolute value of the Walsh transform
     
     Credits: Luca Mariot GitHub (https://github.com/rymoah/CryptoBoolFun/tree/master)
'''
def computeFWT(vector, start, length):
    half = length//2
    #Main cycle: split vector in two parts (v0 e v1), 
    #update v0 as v0=v0+v1, and v1 as v1=v0-v1.
    for i in range(start,start+half):            
        temp = vector[i]
        vector[i] += vector[i+half]
        vector[i+half] = temp - vector[i+half]

    #Recursive call on v0 and v1.
    if(half>1):
        val1 = computeFWT(vector,start,half)
        val2 = computeFWT(vector,start+half,half)

        #At the end of the recursive calls, compare val1 and val2 to decide
        #what is the spectral radius in the portion of truth table included
        #between start and start+length
        if(val1 > val2):                    
            return val1
        else:
            return val2
    else:
        #If we have reached half=1 (function of 2 variables),
        #return the highest coefficient in absolute value.
        if(abs(vector[start]) > abs(vector[start+half])):
            return abs(vector[start])
        else:
            return abs(vector[start+half])


def dot_XOR(v,w):
    dot = 0
    for vi,wi in zip(v,w):
        dot = dot ^ (vi*wi)
    return dot

def hamming(f, g):
    t1 = f.get_truth_table()
    t2 = g.get_truth_table()
    return sum(t1i ^ t2i for t1i, t2i in zip(t1, t2))
    
class BoolFunction:
    def __init__(self,f,n_input):
        self.f = f
        self.n_input = n_input
        #k = self.n_input
        #self.input = tuple(tuple(map(int, bin(j)[2:].zfill(k))) for j in range(2**k))
        #self.truth_table = tuple(self.f(*values) for values in self.input)
        #self.truth_table_polar = tuple(2*self.f(*values)-1 for values in self.input)
    
    def get_input(self):
        k = self.n_input
        return tuple(tuple(map(int, bin(j)[2:].zfill(k))) for j in range(2**k))
    
    def get_truth_table(self):
        return tuple(self.f(*values) for values in self.get_input())
    
    def get_truth_table_polar(self):
        return tuple(2*self.f(*values)-1 for values in self.get_input())
    
    def hamming_weight(self):
        return sum(self.get_truth_table())
    
    def Walsh_transform(self,w):
        return sum((-1)**(self.f(*x)^dot_XOR(w,x)) for x in self.get_input())
    
    def n_max_values_Walsh(self):
        Walsh_values = tuple(abs(self.Walsh_transform(w)) for w in self.get_input())
        max_value = max(Walsh_values)
        n_max_value = sum(tuple(1 for value in Walsh_values if value == max_value))
        return n_max_value
    
    def n_max_values_Walsh_2(self):
        current_max = 0
        count_max = 0
        for w in self.get_input():
            w_value = abs(self.Walsh_transform(w))
            if w_value > current_max:
                current_max = w_value
                count_max = 1
            elif w_value == current_max:
                count_max += 1
        return count_max

    def spectral_radius_old(self):
        return max(abs(self.Walsh_transform(w)) for w in self.get_input())
    
    def spectral_radius_new(self):
        v = np.asarray(self.get_truth_table_polar(),dtype=int)
        r = computeFWT(v, 0, v.size)
        return r
    
    def nonlinearity(self):
        return int(2**(self.n_input-1)- self.spectral_radius_new()/2)
    
    def upperbound_nonlinearity(self):
        n = self.n_input
        if n==5:
            return 12
        elif n==7:
            return 56
        else:
            #Covering Radius Bound
            return 2**(n-1)-2**(n/2-1)
    
    def unbalancedness(self):
        return abs(2**(self.n_input-1)-self.hamming_weight())
                
    def algebraic_degree(self): 
        #c'è qualcosa di sbagliato nel codice! RICONTROLLARE
        v = np.asarray(self.get_truth_table(),dtype=bool)
        d = computeFMT(v, 0, v.size)
        return int(d)
    
    def print_properties(self):
        print('Properties of f (non linearity, hamming weight, spectral radius, unbalancedness)')
        t = PrettyTable()
        t.field_names = ['nl(f)','w_H(f)','W_max(f)','unbal']
        t.add_row([self.nonlinearity(),self.hamming_weight(),self.spectral_radius_new(),self.unbalancedness()])
        print(t)
        

    def print_truth_table(self):
        t = PrettyTable()
        t.add_column(f'(x_1,..,x_{self.n_input})',self.get_input())
        t.add_column('Omega_f',self.get_truth_table())
        t.add_column('W_f',[self.Walsh_transform(w) for w in self.get_input()])
        print(t)

    def __eq__(self, other):
        if hamming(self,other) == 0:
            return True
        return False
    
    def __hash__(self):
        # prova 
        return (hash(self.f) ^ hash(self.n_input))
        
   