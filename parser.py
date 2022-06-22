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
    # Second column is f_nu
    # Not clear if this is Flux or Astro Flux.
    # Maybe best to use the H_nu column.
    #f_n = fw[1]
    # 3rd column is H_nu
    h_nu = fw[2] * Fnu_unit()
    ## The spectrum class is initialized with the eddington flux
    ## in per wavelength unit.
    ## Doing the conversion here
    ## wave dwave = nu dnu
    ## wave = nu dnu/dwave
    ## where dnu/dwave = c / wave**2
    h = (h_nu * const.c / wave**2).to(Flambda_unit())
    return(spectrum(wave,h))
    
