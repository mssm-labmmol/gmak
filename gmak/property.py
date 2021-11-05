import numpy as np
import gmak.runcmd as runcmd
import os
from shutil import copyfile
import gmak.atomic_properties as atomic_properties
from abc import ABC, abstractmethod
from gmak.custom_attributes import CustomizableAttributesMixin


class PropertyDriver(ABC,
                     CustomizableAttributesMixin):

    def __init__(self, name, type, surrmodel, prots, prot_objs,
                 property_class):
        self.name = name
        self.type = type
        self.surrmodel = surrmodel
        self.protocols = prots
        self.protocol_objs = prot_objs
        self.property_cls = property_class

    @abstractmethod
    def create_atomic_properties(self):
        pass

    # Basic implementation that works for simple properties like
    # density, gamma, dg
    def compute(self, gridpoint, protocols):
        # note: protocols : list of BaseProtocol
        protocolObjs = self.protocol_objs
        mu, sigma = protocolObjs[0].get_avg_err_estimate_of_property_for_gridpoint(gridpoint, 
            self.type, self.surrmodel)
        prop_obj = self.property_cls(mu, sigma)
        gridpoint.add_property_estimate(self.name, self.type, prop_obj) 

    def _mask_list_with_none(self, input_list, none_str="none"):
        outlist = []
        for i,p in zip(input_list, self.protocols):
           if p == none_str:
               outlist.append(None)
           else:
               outlist.append(i)
        return outlist


class DensityDriver(PropertyDriver):
    def create_atomic_properties(self):
        return self._mask_list_with_none([
            atomic_properties.create_atomic_property("density")
        ])

class DHvapDriver(PropertyDriver):
    def create_atomic_properties(self):
        try:
            mu = self.get_custom_attribute('mu')
        except AttributeError:
            mu = None
        try:
            alpha = self.get_custom_attribute('alpha')
        except AttributeError:
            alpha = None
        return self._mask_list_with_none([
            atomic_properties.create_atomic_property("potential"),
            atomic_properties.create_atomic_property("potential"),
            atomic_properties.create_atomic_property("polcorr",
                                                     mu=mu,
                                                     alpha=alpha)
        ])


    def compute(self, gridpoint, protocols):
        protocolObjs = self.protocol_objs
        propertyTypes = [ap.name if ap is not None else None for ap in
                         self.create_atomic_properties()]

        if propertyTypes[0] is not None:
            mu_u_liq, sigma_u_liq = protocolObjs[0].get_avg_err_estimate_of_property_for_gridpoint(gridpoint, 
                propertyTypes[0], self.surrmodel)
        else:
            mu_u_liq, sigma_u_liq = 0.0, 0.0

        if propertyTypes[1] is not None:
            mu_u_gas, sigma_u_gas = protocolObjs[1].get_avg_err_estimate_of_property_for_gridpoint(gridpoint, 
                propertyTypes[1], self.surrmodel)
        else:
            mu_u_gas, sigma_u_gas = 0.0, 0.0

        if propertyTypes[2] is not None:
            mu_u_polcorr, sigma_u_polcorr = protocolObjs[2].get_avg_err_estimate_of_property_for_gridpoint(gridpoint, 
                propertyTypes[2], self.surrmodel)
        else:
            mu_u_polcorr, sigma_u_polcorr = 0.0, 0.0

        # get temperature from $compute block if it was supplied;
        # otherwise, try to get it from the liquid protocol
        try:
            temperature = self.get_custom_attribute('temperature')
        except AttributeError:
            try:
                temperature = protocolObjs[0].get_temperature()
            except:
                raise Exception("Couldn't figure out temperature for calculating DHvap!")

        # get corrections from $compute block; assume zero if not
        # supplied.
        try:
            corrs = self.get_custom_attribute('C')
        except AttributeError:
            corrs = 0.0

        try:
            nmols = self.get_custom_attribute('nmols')
        except AttributeError:
            raise Exception("Couldn't figure out nmols for calculating DHvap!")
            
        prop_obj = self.property_cls(
            mu_u_liq,
            mu_u_gas,
            mu_u_polcorr,
            sigma_u_liq,
            sigma_u_gas,
            sigma_u_polcorr,
            corrs,
            nmols,
            temperature)
        gridpoint.add_property_estimate(
            self.name,
            self.type,
            prop_obj)


class GammaDriver(PropertyDriver):
    def create_atomic_properties(self):
        return self._mask_list_with_none([
            atomic_properties.create_atomic_property("gamma")
        ])


class DGDriver(PropertyDriver):
    def create_atomic_properties(self):
        # get temperature from $compute block if it was supplied;
        # otherwise, try to get it from the liquid protocol
        try:
            temperature = self.get_custom_attribute('temperature')
        except AttributeError:
            try:
                temperature = self.protocol_objs[0].get_temperature()
            except:
                raise Exception("Couldn't figure out temperature for calculating DG!")
        return self._mask_list_with_none([
            atomic_properties.create_atomic_property("dg", temperature)
        ])

class CustomDriver(PropertyDriver):
    def create_atomic_properties(self):
        return self._mask_list_with_none([
            # NOTE: The use of `vars` below means that any custom attribute
            # that is not necessary for the property will raise an error!
            atomic_properties.create_atomic_property(self.type,
                                                     **vars(self.get_custom_attributes()))
        ])

    def compute(self, gridpoint, protocols):
        # note: protocols : list of BaseProtocol
        protocolObjs = self.protocol_objs
        mu, sigma = protocolObjs[0].get_avg_err_estimate_of_property_for_gridpoint(gridpoint,
            self.type, self.surrmodel)
        prop_obj = self.property_cls(self.type, mu, sigma)
        gridpoint.add_property_estimate(self.name, self.type, prop_obj)


class PropertyDriverFactory:

    @classmethod
    def create(cls, name, type, surrmodel, prots, prot_objs):
        if type == 'density':
            return DensityDriver(name, type, surrmodel, prots, prot_objs,
                                 Density)
        elif type == 'dhvap':
            return DHvapDriver(name, type, surrmodel, prots, prot_objs, dHvap)
        elif type == 'gamma':
            return GammaDriver(name, type, surrmodel, prots, prot_objs, Gamma)
        elif type == 'dg':
            return DGDriver(name, type, surrmodel, prots, prot_objs,
                            DGAlchemicalAnalysis)
        else:
            return CustomDriver(name, type, surrmodel, prots, prot_objs,
                                CustomProperty)

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
        # custom property
        return CustomProperty(property_string, value, err)

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

    id_string = "density"

    def __init__ (self, value, err):
        self.set_textual_elements ("density", "kg m$^{-3}$", "$\\rho_\\mathrm{liq}$")
        self.value = value
        self.err = err

class Gamma (PropertyBase):

    id_string = "gamma"

    def __init__ (self, value, err, nsurf=2):
        self.set_textual_elements ("gamma", "mN m$^{-1}$", "$\\gamma$")
        self.value = 0.1 * value/nsurf
        self.err = 0.1 * err/nsurf

class dHvap (PropertyBase):

    id_string = "dhvap"

    # corrs are taken into account as dHvap_corr = dHvap - value_pol + corrs
    def __init__ (self, value_liq, value_gas, value_pol, err_liq, err_gas, err_pol,\
            corrs, nmols, temperature):
        R = 8.3144598e-3
        self.set_textual_elements ("dhvap", "kJ mol$^{-1}$", "$\\Delta H_\\mathrm{vap}$")
        self.value = value_gas - value_liq/nmols - value_pol + corrs + R*temperature
        self.err   = np.sqrt(err_gas**2 + (err_liq/nmols)**2 + (err_pol)**2)
        #print( "************************************")
        #print( "Initialized dHvap = %.2f +/- %.2f " % (self.value, self.err))
        #print( "Components:")
        #print( "U_liq = %.2f (%.2f)" % (value_liq, value_liq/nmols))
        #print( "U_gas = %.4f " % (value_gas))
        #print( "Polcorr = %.2f " % (value_pol))
        #print( "Other corrections = %.2f " % (corrs))
        #print( "RT = %.2f " % (R*temperature))

class CohesiveEnergyDensity(PropertyBase):

    id_string = "ced"
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

    id_string = "gced"

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

class DGAlchemicalAnalysis(PropertyBase):

    id_string = "dg"

    def __init__(self, value, err):
        # Set attributes
        self.value = value
        self.err = err
        self.set_textual_elements("dg",
                                  "kJ mol$^{-1}$",
                                  "$\\Delta G_{\\mathrm{solv}}$")



class CustomProperty(PropertyBase):
    # NOTE: A CustomProperty is always associated with a CustomAtomicProperty.
    def __init__ (self, name, value, err):
        # Just to check if it works.
        _ = atomic_properties.\
            CustomAtomicPropertyFactory.\
            create_custom_atomic_property(name)
        # Now do stuff.
        self.set_textual_elements(name, f"{name} units", name)
        self.value = value
        self.err = err

