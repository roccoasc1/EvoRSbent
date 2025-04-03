def fit_w1nl_w2unbal(bfunc,w1,w2,normalize=True):
    # nonlinearity - unbalanceness # k dispari
    k = bfunc.n_input
    if normalize:
        max_nl = bfunc.upperbound_nonlinearity()
        max_unbal = 2**(k-1)
        return w1*bfunc.nonlinearity()/max_nl - w2*bfunc.unbalancedness()/max_unbal
    else:
        return w1*bfunc.nonlinearity() - w2*bfunc.unbalancedness()
    
def fit_nl_unbal_normalize(bfunc):
    return fit_w1nl_w2unbal(bfunc,w1=1.0,w2=1.0,normalize=True)

def fit_nl(bfunc):
    # nonlinearity --> the functions with max non linearity for n even are called bent functions
    # same as fit_w1nl_w2unbal with w1=1, w2=0, normalize=False
    return bfunc.nonlinearity()

def fit_nl_unbal(bfunc):
    # nonlinearity - unbalanceness # k dispari
    return int(bfunc.nonlinearity()-bfunc.unbalancedness())

#def fit_nl_unbal_kr(bfunc):
    # k dispari
#    if bfunc.unbalancedness() <= 2**(bfunc.n_input-2):
#        s = 1 - bfunc.unbalancedness()/2**(bfunc.n_input-2) # s should give the percentange of how far we are from being balanced
#        return s*bfunc.nonlinearity()
#    return 0

def fit_2(bfunc):
    # what does the part after - mean? # non linearity - (deviazione della hamming distance dall'essere bent)
    # solo k pari # ban
    return bfunc.nonlinearity()-abs(2**(bfunc.n_input-1)-2**(bfunc.n_input/2-1)-bfunc.hamming_weight())

def fit_kroneker(bfunc):
    # con questa fitness escludiamo tutti gli rs non bilanciati
    if bfunc.unbalancedness()==0:
        return bfunc.nonlinearity()
    return 0

def fit_max_values(bfunc):
    unbal = bfunc.unbalancedness()
    if unbal==0:
        return bfunc.nonlinearity()+(2**bfunc.n_input-bfunc.n_max_values_Walsh())/2**bfunc.n_input
    return -unbal

def fit_nl_max_values(bfunc):
    return bfunc.nonlinearity()+(2**bfunc.n_input-bfunc.n_max_values_Walsh())/2**bfunc.n_input

def fit_nl_unbal_kroneker(bfunc):
    unbal = bfunc.unbalancedness()
    if unbal==0:
        return bfunc.nonlinearity()
    return -unbal
