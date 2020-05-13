
def initialize_from_input (input_file):
    output_grid = ParameterGrid ()
    output_protocols = []
    output_properties = []
    output_protocolsHash = {}
    output_workdir = ""
    output_optimizer = gridOptimizer ()
    output_gridshifter = GridShifter ()
    output_paramLoop = ParameterLoop ()
    output_gaCoverInterface = coverInterface ()
    output_reweightHash = {}
    fp = open(input_file, "r")
    for line in fp:
        if line[0] == '#':
            continue
        # grid
        if len(line.split()) < 1:
            continue
        if (line.split()[0] == "workdir"):
            output_workdir = os.path.abspath(line.split()[1])
        if (line.rstrip() == "$grid"):
            output_grid.read_from_stream (fp)
        if (line.rstrip() == "$protocol"):
            line = next(fp)
            if (line.split()[0] == 'type'):
                typeRead = line.split()[1]
                if (typeRead == 'liquid'):
                    new_protocol = LiquidProtocol("",0,"",[""],[],[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'gas'):
                    new_protocol = GasProtocol("","",0.0,0.0,[],[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                elif (typeRead == 'slab'):
                    new_protocol = SlabProtocol("",[],5.0,[])
                    new_protocol.read_from_stream (fp)
                    output_protocols.append(new_protocol)
                else:
                    print \
"""ERROR: Type \"%s\" is not supported.\n""" % typeRead
                    exit()
            else:
                print \
"""ERROR: First line after a $protocol flag MUST assign its type.\n"""
                exit()
        if (line.rstrip() == '$compute'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    break
                propRead = line.split()[0]
                nameRead = line.split()[1]
                output_properties.append(propRead)
                if (propRead == 'density'):
                    output_protocolsHash[propRead] = nameRead
                    # find protocol with name given
                    protocols = filter (lambda x: x.name == nameRead, \
                            output_protocols)
                    for protocol in protocols:
                        protocol.properties.append('density') 
                        # potential ALWAYS
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')
                elif (propRead == 'dhvap'):
                    nameLiq = line.split()[1]
                    nameGas = line.split()[2]
                    corr    = float(line.split()[3])

                    output_protocolsHash[propRead] = [nameLiq,nameGas]
                    # find protocol with name given - Liq
                    protocols = filter (lambda x: x.name == nameLiq, \
                            output_protocols)
                    for protocol in protocols:
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')

                    # find protocol with name given - Gas
                    protocols = filter (lambda x: x.name == nameGas, \
                            output_protocols)
                    for protocol in protocols:
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')
                        protocol.properties.append('polcorr')
                        protocol.set_other_corrections(corr)

                elif (propRead == 'gamma'):
                    # find protocol with name given
                    output_protocolsHash[propRead] = nameRead
                    protocols = filter (lambda x: x.name == nameRead, \
                            output_protocols)
                    for protocol in protocols:
                        protocol.properties.append('gamma')
                        # potential ALWAYS
                        if 'potential' not in protocol.properties:
                            protocol.properties.append('potential')
                else:
                    print \
"""ERROR: Property \"%s\" is not supported.\n""" % typeRead
                    exit()
        if (line.rstrip() == '$optimize'):
            output_optimizer.readFromStream (fp)
        if (line.rstrip() == '$parameters'):
            output_paramLoop.readFromStream (fp)
        if (line.rstrip() == '$gridshift'):
            output_gridshifter.readFromStream (fp)
        if (line.rstrip() == '$gacover'):
            output_gaCoverInterface.readFromStream (fp)
        if (line.rstrip() == '$reweight'):
            for line in fp:
                if line[0] == '#':
                    continue
                if line.rstrip() == '$end':
                    break
                option = line.split()[0]
                value  = line.split()[1]
                output_reweightHash[option] = value
    fp.close()
    return (output_workdir, output_grid, output_protocols, output_properties, \
            output_protocolsHash, output_optimizer, output_gridshifter,\
            output_paramLoop, output_gaCoverInterface, output_reweightHash)
