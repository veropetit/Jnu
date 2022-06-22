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
        
        
    def get_F_lambda(self):
        return(self.H_lambda*4*np.pi)
        
    def get_Fa_lambda(self):
        return(self.H_lambda*4)
        
    def get_H_nu(self):
        return((self.H_lambda*self.wave**2/const.c).to(Fnu_unit()))
        
    def get_F_nu(self):
        return(self.get_H_nu()*4*np.pi)
        
    def get_Fa_nu(self):
        return(self.get_H_nu()*4)
        
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

def read_fastwind(file):
    
    fw = np.genfromtxt(file, unpack=True)
    # first column is wavelength in AA
    wave = fw[0]*u.AA
    # Second column is related to log(f_nu)
    # and to go to Hnu, we need to do
    # (according to the idlproc/readmod.pro):
    correction = np.log10(120.0**2/(4*np.pi))
    H_nu =  10**(fw[1]+correction) * Fnu_unit()
    # This correcpond to the 10**3rd column of the .txt files.

    ## The spectrum class is initialized with the eddington flux
    ## in per wavelength unit.
    ## Doing the conversion here
    ## wave dwave = nu dnu
    ## wave = nu dnu/dwave
    ## where dnu/dwave = c / wave**2

    H_lambda = (H_nu * const.c / wave**2).to(Flambda_unit())
    return(spectrum(wave,H_lambda))
