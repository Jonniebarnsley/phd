import numpy as np

def round_sig_figs(x, sig=2):
    return round(x, sig-int(np.floor(np.log10(abs(x))))-1)