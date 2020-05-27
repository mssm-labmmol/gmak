import numpy as np

def init_property_from_string(property_string, value, err):
    if (property_string == 'density'):
        return Density(value, err)
    elif (property_string == 'dhvap'):
        return dHvap(0, value, 0, 0, err, 0, 0, 1, 0)
    elif (property_string == 'gamma'):
        return Gamma(value, err)
    else:
        raise ValueError("Property {} is not recognized.".format(property_string))
    
class PropertyBase:

    def __init__ (self):
        self.name  = "property"
        self.units = "u.a."
        self.symbol = "$P$"
        self.value = 0.0
        self.err   = 0.0

    def set_textual_elements (self, name, units, symbol):
        self.name = name
        self.units = units
        self.symbol = symbol

    def set_value (self, value, err):
        self.value = value
        self.err   = err

    def get_label (self):
        return "%s [%s]" % (self.symbol, self.units)

    def get_label_err (self):
        return "$\\delta$%s [%s]" % (self.symbol, self.units)

class Density (PropertyBase):

    def __init__ (self, value, err):
        self.set_textual_elements ("density", "kg m$^{-3}$", "$\\rho_\\mathrm{liq}$")
        self.value = value
        self.err = err

class Gamma (PropertyBase):

    def __init__ (self, value, err):
        self.set_textual_elements ("gamma", "mN m$^{-1}$", "$\\gamma$")
        self.value = value
        self.err = err

class dHvap (PropertyBase):

    # corrs are taken into account as dHvap_corr = dHvap - value_pol + corrs
    def __init__ (self, value_liq, value_gas, value_pol, err_liq, err_gas, err_pol,\
            corrs, nmols, temperature):
        R = 8.3144598e-3
        self.set_textual_elements ("dhvap", "kJ mol$^{-1}$", "$\\Delta H_\\mathrm{vap}$")
        self.value = value_gas - value_liq/nmols - value_pol + corrs + R*temperature
        self.err   = np.sqrt(err_gas**2 + (err_liq/nmols)**2 + (err_pol)**2)
        print( "************************************")
        print( "Initialized dHvap = %.2f +/- %.2f " % (self.value, self.err))
        print( "Components:")
        print( "U_liq = %.2f (%.2f)" % (value_liq, value_liq/nmols))
        print( "U_gas = %.2f " % (value_gas))
        print( "Polcorr = %.2f " % (value_pol))
        print( "Other corrections = %.2f " % (corrs))
        print( "RT = %.2f " % (R*temperature))
