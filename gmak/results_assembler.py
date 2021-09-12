from numpy import all, array
from numpy import mean
from numpy.linalg import norm

class ResultsAssembler:

    def __init__(self, parameter_names):
        self.parameters = []
        self.property_names = []
        self.parameter_names = parameter_names
        self.data = []

    def _checkParameters(self, parameter_values):
        out = -1
        for i, x in enumerate(self.parameters):
            if all(x == parameter_values):
                out = i
                break
        return out
            
    def addData(self, parameter_values, property_name, property_estimate, property_error):
        if (property_name not in self.property_names):
            self.property_names.append(property_name)

        p = self._checkParameters(parameter_values)
        if (p == -1):
            self.parameters.append(parameter_values)
            p = len(self.parameters) - 1
            self.data.append({})

        try:
            # updating data
            self.data[p][property_name]['estimate'] = property_estimate
            self.data[p][property_name]['error'] = property_error
        except KeyError:
            # setting new data
            self.data[p][property_name] = {'estimate': property_estimate, 'error': property_error}

    def writeToFile(self, prefix):
        for prop in self.property_names:
            fn = prefix + '_{}.csv'.format(prop)
            fp = open(fn, 'w')
            # Write title
            fp.write(" ".join(self.parameter_names) + " {} {}_err\n".format(prop, prop))
            for i in range(len(self.parameters)):
                # Write parameter values
                for j in range(len(self.parameter_names)):
                    fp.write("%14.6e" % self.parameters[i][j])
                # Write property value and error
                fp.write("%14.6f%14.6f\n" % (self.data[i][prop]['estimate'], self.data[i][prop]['error']))
            fp.close()

    def print(self, stream):
        stream.write(str(self.parameters) + '\n')
        stream.write(str(self.data) + '\n')
