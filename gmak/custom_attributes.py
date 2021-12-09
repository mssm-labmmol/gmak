class CustomizableAttributesMixin:

    class InputParameters:
        """
        An instance of this class stores the input parameters supplied in a
        block of the input file as attributes. Input parameters that correspond
        to a single value token are stored as objects with type appropriate to
        the token value. Input parameters that correspond to a list of multiple
        value tokens are stored as a list of objects, each one processed as above.

        For example, the block

        .. code-block:: gmi

            $coords
            name liquid
            type gmx_liquid
            coords molecule.gro
            nmols 1500
            box 3.0 3.0 3.0
            $end

        corresponds to an object ``coords_attr`` for which

            >>> coords_attr.name
            'liquid'
            >>> coords_attr.type
            'gmx_liquid'
            >>> coords_attr.coords
            'molecule.gro'
            >>> coords_attr.nmols
            1500
            >>> coords_attr.box
            [3.0, 3.0, 3.0]
        
        Note the automatic type conversion for the parameters ``nmols`` and
        ``box``, from :py:class:`string` and :py:class:`list` of
        :py:class:`string` to :py:class:`int` and :py:class:`list` of
        :py:class:`float`, respectively.

        .. note:: Since the parameter is used as the name of an instance
            attribute, the use of only lowercase characters and underscores for
            custom parameters is highly recommended to avoid run-time errors.
       """

        def __init__(self):
            pass

    # "private"
    def _set_custom_attribute(self, attr_name, attr_value):
        # custom attributes are attributes of
        # self._custom_attributes_data
        if not hasattr(self, "_custom_attributes_data"):
            self._custom_attributes_data = self.InputParameters()
        try:
            if len(attr_value) == 1:
                setattr(self._custom_attributes_data, attr_name, attr_value[0])
            else:
                setattr(self._custom_attributes_data, attr_name, attr_value)
        except TypeError:
            # then it has no length
            setattr(self._custom_attributes_data, attr_name, attr_value)

    def set_custom_attributes_from_blockdict(self, blockdict, systems,
                                             coordinates, protocols,
                                             exclude=None):
        self._custom_attributes_data = self.InputParameters()
        for k, v in blockdict.items():
            if (exclude is None) or (k not in exclude):
                try:
                    if (len(v) == 3) and (v[0] == "from"):
                    # syntax for recycling from other objects
                        if v[1] == "system":
                            origin = systems[v[2]]
                        elif v[1] == "coordinates":
                            origin = coordinates[v[2]]
                        elif v[1] == "protocol":
                            origin = protocols[v[2]]
                        else:
                            raise ValueError(f"Can't merge attribute from {v[1]}.")
                        self.merge_attr(origin, k)
                    else:
                        self._set_custom_attribute(k, v)
                except TypeError:
                    # v is not an iterable
                    self._set_custom_attribute(k, v)
                except KeyError:
                    # attribute could not be found
                    raise KeyError(f"Can't find attribute \"{k}\" in {v[1]} "
                                   f"\"{v[2]}\".")


    def clone_custom_attributes(self, dest):
        from copy import deepcopy
        try:
            dest._custom_attributes_data = deepcopy(
                self._custom_attributes_data)
        except AttributeError:
            pass

    def get_custom_attributes(self):
        try:
            return self._custom_attributes_data
        except AttributeError:
            return None

    def get_custom_attribute(self, name):
        # this will raise a proper error
        return getattr(self._custom_attributes_data, name)

    def merge_attr(self, other, attr_name):
        self._set_custom_attribute(attr_name,
                                   other.get_custom_attribute(attr_name))

