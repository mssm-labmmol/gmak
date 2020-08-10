import numpy as np

def init_property_from_string(property_string, value, err):
    if (property_string == 'density'):
        return Density(value, err)
    elif (property_string == 'dhvap'):
        return dHvap(0, value, 0, 0, err, 0, 0, 1, 0)
    elif (property_string == 'gamma'):
        return Gamma(value, err)
    elif (property_string == 'ced'):
        # cohesive energy dummy
        ced = CohesiveEnergyDensity(0, 1.0, 1.0, 0, 0, 1.0, 0, 0, 0, 1)
        # modify
        ced.value = value
        ced.err = err
        return ced
    elif (property_string == 'gced'):
        # gamma via cohesive energy
        gamma_ced = GammaViaCohesiveEnergyDensity(0, 1.0, 1.0, 0, 0, 1.0, 1.0, 0, 0, 1)
        # modify
        gamma_ced.value = value
        gamma_ced.err   = err
        return gamma_ced
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

    def __init__ (self, value, err, nsurf=2):
        self.set_textual_elements ("gamma", "mN m$^{-1}$", "$\\gamma$")
        self.value = 0.1 * value/nsurf
        self.err = 0.1 * err/nsurf

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
        print( "U_gas = %.4f " % (value_gas))
        print( "Polcorr = %.2f " % (value_pol))
        print( "Other corrections = %.2f " % (corrs))
        print( "RT = %.2f " % (R*temperature))

class CohesiveEnergyDensity(PropertyBase):

    avogadroConstant = 6.02214076e+23

    def __init__(self, u_liq, u_gas, v, pol, err_liq, err_gas, err_v, err_pol, corrs, nmols):
        v_m = self.avogadroConstant * (v/nmols) * 1e-27 # in m3
        err_v_m = self.avogadroConstant * (err_v/nmols) * 1e-27 # in m3
        self.set_textual_elements("ced", "MPa", "$\\delta^2$")
        self.value = 1e-3 * (u_gas - u_liq/nmols - pol + corrs)/(v_m)
        self.err = self.value * np.sqrt((err_v_m/v_m)**2 + ((np.sqrt(err_gas**2 + (err_liq/nmols)**2 + err_pol**2))/(u_gas - u_liq/nmols - pol + corrs))**2)
        #print("**** COHESIVE ENERGY DENSITY ****")
        #print("molar volume (cm^3/mol) = {}".format(v_m * 1e6))
        #print("error (cm^3/mol)        = {}".format(err_v_m*1e6))
        #print("delta2 (MPa) = {}".format(self.value))
        #print("error (MPa) = {}".format(self.err))

class GammaViaCohesiveEnergyDensity(PropertyBase):

    # see https://www.sciencedirect.com/science/article/abs/pii/0021979772902457
    # we use the "quasi-thermodynamic" value derived in the last section of this paper, converted to our units
    # (in the paper, delta2 is in cal/cm3, gamma is is dynes/cm and vm is in cm3/mol)
    conversionConstant = 0.01 * 14.041 * (2.045**2)

    def __init__(self, u_liq, u_gas, v, pol, err_liq, err_gas, err_v, err_pol, corrs, nmols):
        ced = CohesiveEnergyDensity(u_liq, u_gas, v, pol, err_liq, err_gas, err_v, err_pol, corrs, nmols)
        v_m = ced.avogadroConstant * (v/nmols) * 1e-27 # in m3
        err_v_m = ced.avogadroConstant * (err_v/nmols) * 1e-27 # in m3
        self.set_textual_elements("gced", "mN m$^{-1}$", "$\\gamma_{\\delta^2}$")
        self.value = (v_m**(1.0/3)) * ced.value / self.conversionConstant
        #self.err   = np.abs( ced.value * (1.0/3)*(v_m)**(-2.0/3)*err_v_m + ced.err * (v_m**(1.0/3)) ) / self.conversionConstant
        self.err   = self.value * np.sqrt( (((1.0/3) * v_m**(-2.0/3) * err_v_m)/(v_m**(1./3)))**2 + (ced.err/ced.value)**2)
        
        
