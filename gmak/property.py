import numpy as np
import gmak.runcmd as runcmd
import os
from shutil import copyfile
import gmak.component_properties as component_properties
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
    def create_component_properties(self):
        pass

    # Basic implementation that works for simple properties like
    # density, gamma, dg
    def compute(self, gridpoint):
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
    def create_component_properties(self):
        return self._mask_list_with_none([
            component_properties.create_component_property("density")
        ])

class DHvapDriver(PropertyDriver):
    def create_component_properties(self):
        try:
            mu = self.get_custom_attribute('mu')
        except AttributeError:
            mu = None
        try:
            alpha = self.get_custom_attribute('alpha')
        except AttributeError:
            alpha = None
        if (mu is None) or (alpha is None):
            return self._mask_list_with_none([
                component_properties.create_component_property("potential"),
                component_properties.create_component_property("potential"),
                None,
            ])
        else:
            return self._mask_list_with_none([
                component_properties.create_component_property("potential"),
                component_properties.create_component_property("potential"),
                component_properties.create_component_property("polcorr",
                                                         mu=mu,
                                                         alpha=alpha)
            ])


    def compute(self, gridpoint):
        protocolObjs = self.protocol_objs
        propertyTypes = [ap.name if ap is not None else None for ap in
                         self.create_component_properties()]

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
    def create_component_properties(self):
        return self._mask_list_with_none([
            component_properties.create_component_property("gamma")
        ])


class DGDriver(PropertyDriver):

    def create_component_properties(self):
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
            component_properties.create_component_property("dg", temperature)
        ])


class CustomDriver(PropertyDriver):

    def __init__(self, name, type, surrmodel, prots, prot_objs,
                 component_props, system_objs, property_factory):

        super().__init__(name, type, surrmodel, prots, prot_objs,
                         property_factory)
        # find desired systems
        # this list is filled in the same order of the protocols in prot_objs
        self.system_objs = []
        for prot in prot_objs:
            for system in system_objs:
                if system.name == prot.system:
                    self.system_objs.append(system)
                    break
        if component_props is None:
            raise ValueError(f"Component properties not specified for composite"
                             f"property {name}.")
        else:
            self.component_props = component_props

    def create_component_properties(self):
        out = []
        for prop_type, system_obj in zip(self.component_props,
                                         self.system_objs):
            ap = component_properties.create_component_property(prop_type,
                                                          system=system_obj)
            self.clone_custom_attributes(ap)
            out.append(ap)
        return self._mask_list_with_none(out)

    def compute(self, gridpoint):
        protocolObjs = self.protocol_objs
        mus = []
        sigmas = []
        for protocolObj, propType in zip(self.protocol_objs,
                                         self.component_props):
            mu, sigma = protocolObj.get_avg_err_estimate_of_property_for_gridpoint(
                gridpoint, propType, self.surrmodel)
            mus.append((propType, mu))
            sigmas.append((propType, sigma))

        prop_obj = self.property_cls.create(self.type, mus, sigmas)
        self.clone_custom_attributes(prop_obj)
        prop_obj.compute_values()
        gridpoint.add_property_estimate(self.name, self.type, prop_obj)


class PropertyDriverFactory:

    @classmethod
    def create(cls, name, type, surrmodel, prots, prot_objs, system_objs,
               component_props=None):
        # component_props is None for all drivers except the custom
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
                                component_props, system_objs, CustomPropertyFactory)

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

class DGAlchemicalAnalysis(PropertyBase):

    id_string = "dg"

    def __init__(self, value, err):
        # Set attributes
        self.value = value
        self.err = err
        self.set_textual_elements("dg",
                                  "kJ mol$^{-1}$",
                                  "$\\Delta G_{\\mathrm{solv}}$")


class CustomProperty(PropertyBase,
                     CustomizableAttributesMixin):
    def __init__ (self, name, values, errs, calculator=None):
        """
        values : list of (protype, value) pairs for each component.
        errs : list of (protype, value) pairs for each component.
        """
        self.set_textual_elements(name, f"{name} units", name)
        if calculator is None:
            self.calculator = self._default_calculator
        else:
            self.calculator = calculator
        self._values = values
        self._errs = errs

    def compute_values(self):
        self.value, self.err = self.calculator(self._values,
                                               self._errs,
                                               self.get_custom_attributes())

    @staticmethod
    def _default_calculator(values, errs, propertry_attrs):
        if len(values) > 1:
            raise ValueError(f"If the composite property has more than one"
                             " component property, you most supply a proper"
                             " calculator function.")
        return float(values[0][1]), float(errs[0][1])


class CustomPropertyFactory:

    ptable = {}

    @classmethod
    def add_custom_composite_property(cls, type_name, calculator=None):
        cls.ptable[type_name] = calculator

    @classmethod
    def create(cls, type_name, values, errs):
        return CustomProperty(type_name, values, errs, cls.ptable[type_name])


