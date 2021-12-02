class CustomizableAttributesMixin:

    class CustomizableAttributesData:
        """
        This class stores the data that is passed to the program's customizable
        objects (e.g., protocols) in the input file. Every such object contains
        an instance of this class that is passed to the functions implemented
        and supplied by the user *via* the customization API.

        In the context of an input-file block defining a given customizable
        object, every line follows the syntax::

            ATTRNAME    VALUES

        ``ATTRNAME`` is a string containing only lowercase characters and
        underscores. It is separated by one or more white spaces from
        ``VALUES``, which is either a string not containing white spaces or a
        sequence of those (henceforth referred to as tokens) separated by any
        number of spaces. The lines below show examples with proper syntax::

            extend z 5.0
            total_mass 2.3e+03
            index_group /path/to/index.ndx OW

        The data defined in the input-file line is made available as an
        instance attribute with name ``ATTRNAME``. The value of the attribute
        depends on the syntax of ``VALUES``. If there is a single numerical
        token, the value is the numerical value of the corresponding string
        (type conversion is done automatically). If there is a single
        non-numerical token, it is the corresponding string itself. If there is
        a sequence of tokens, it is a list of values resulting from the
        processing of each token as described above, except if the first token
        equals ``from`` (more on that exception below).  For the examples
        above::

            >>> attr.extend
            ['z', 5.0]
            >>> attr.total_mass
            2.3e+03
            >>> attr.index_group
            ['/path/to/index.ndx', 'OW']


        where ``attr`` stands for the defined instance of type
        :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.CustomizableAttributesData`.

        .. note:: Since ``ATTRNAME`` is used as the name of an instance attribute,
           the use of only lowercase characters and underscores is highly
           recommended to avoid run-time errors.

        It is possible to recycle custom attributes from other objects using
        input lines with the following syntax:

        .. code-block:: text

           ATTRNAME from OBJECT OBJNAME

        ``OBJNAME`` and ``OBJECT`` are the name and type (either ``system``,
        ``coordinates`` or ``protocol``) of the object from which
        the attribute is recycled, and ``ATTRNAME`` is the name of the recycled
        attribute. For example, the line

        .. code-block:: text

           nmols from coordinates water_liq

        sets in the object being defined a custom attribute ``nmols`` with the
        same value as that for an existing coordinates object named
        ``water_liq``.
        """
        def __init__(self):
            pass

    # "private"
    def _set_custom_attribute(self, attr_name, attr_value):
        # custom attributes are attributes of
        # self._custom_attributes_data
        if not hasattr(self, "_custom_attributes_data"):
            self._custom_attributes_data = self.CustomizableAttributesData()
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
        self._custom_attributes_data = self.CustomizableAttributesData()
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

