print("Importing jnu.parser")

import numpy as np
import astropy.units as u
import astropy.constants as const

def Fnu_unit():
  return(1*u.erg/u.s/u.cm**2/u.Hz)

def Flambda_unit():
  return(1*u.erg/u.s/u.cm**2/u.AA)

class spectrum:
    '''
    Holds a spectrum
    '''

    def __init__(self, wave, H_lambda):
        '''Initializing the spectrum with the Eddington flux in per AA units'''
        
        self.wave = wave
        self.H_lambda = H_lambda
        
        
    def to_F_lambda(self):
        return(self)
        
    def to_Fa_lambda(self):
        return(self)
        
    def to_H_nu(self):
        return(self)
        
    def to_F_nu(self):
        return(self)
        
    def to_Fa_nu(self):
        return(self)
        
    def get_nu(self):
        '''Return the wavelength array converted to frequency (Hz)'''
        return((const.c/self.wave).to(u.Hz))
        
    # Overload __getitem and __setitem
    
def read_tlusty(file):

  tlusty = np.genfromtxt(file, unpack=True) # unnormalized TLUSTY spectrum
  ## The file gives wavelength in AA
  wave = tlusty[0]* u.AA
  ## The file gives the Eddington flux in wavelength (AA)
  h = tlusty[1]* Flambda_unit()
  ## The spectrum class is initialized with the eddington flux, so OK.
  return(spectrum(wave, h))
