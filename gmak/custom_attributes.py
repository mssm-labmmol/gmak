class CustomizableAttributesMixin:

    class CustomizableAttributesData:
        def __init__(self):
            pass

    # "private"
    def _set_custom_attribute(self, attr_name, attr_value):
        # custom attributes are attributes of
        # self._custom_attributes_data
        if not hasattr(self, "_custom_attributes_data"):
            self._custom_attributes_data = self.CustomizableAttributesData()
        setattr(self._custom_attributes_data, attr_name, attr_value)

    def set_custom_attributes_from_blockdict(self, blockdict):
        for k, v in blockdict.items():
            self._set_custom_attribute(k, v)

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
