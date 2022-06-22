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


class FASTWIND_model :
    
    def __init__ (self, spfile, outfile, jfile) :
        '''Initialize a FASTWIND_model from a set of 3 files
        
            :param spfile: a 3 column file containing wavelength (AA), log(Fnu), log(Hnu)
            :param outfile: file containing a line profile
            :param jfile: file containing Jnu(r) for a set of wavelengths
            :return: a FASTWIND_model object.
        '''
    
        self.init_spfile (spfile)
        self.init_outfile (outfile)
        self.init_jfile (jfile)
        return

    def init_spfile (self, spfile) :
        '''Init the FASTWIND_model.spectrum
        
            :param file: the spectrum file
            :return: a spectrum object
            '''
        spectrum = read_fastwind(spfile)
        self.spectrum = spectrum
        return

    def init_outfile (self, outfile) :
        self.wavelength_OVI, self.spectrum_OVI = np.genfromtxt (outfile, unpack=True, skip_footer=1, usecols=(2,4))
        return

    def init_jfile (self, jfile) :
        fh = open (jfile, 'rt')
        lines = fh.readlines ()
        fh.close ()
        blocks = []
        imin = 0
        for i in range (len(lines)) :
            if not lines[i].split () :
                blocks.append (lines[imin:i])
                imin = i+1
        self.radius = np.array ([])
        self.jwavelengths = np.array ([])*u.AA
        self.jnu = np.array ([])
        for i in range (len (blocks)) :
            wl, r, j = self.parse_block (blocks[i])
            if i == 0 :
                self.jwavelengths = wl*u.AA
                self.radius = r
                ## Need to add the proper Jnu units here.
                self.jnu = j
            else :
                self.jwavelengths = np.append (self.jwavelengths, wl*u.AA)
                self.jnu = np.append (self.jnu, j, axis=0) # check correct axis
        
        
    def parse_block (self, block) :
        wlstring = block[2].replace ('|', ' ')
        wl = np.fromstring (wlstring, sep=' ')
        datalines = []
        for i in range (len(block)) :
            if i > 2 :
                temp = block[i].replace ('|', ' ')
                datalines.append (np.fromstring (temp, sep=' '))
        data = np.transpose (np.array (datalines))
        r = data[1,:]
        jdata = data[2:,:]
        return wl, r, jdata
            
