from gmak.systems import TopologyOutput, System
from gmak.custom_attributes import CustomizableAttributesMixin
from gmak.interaction_parameter import (InteractionParameter, InteractionPair,
                                        InteractionAtom,
                                        InteractionParameterType,
                                        CombinationRule)
from typing import List
from io import StringIO
from glob import glob
from subprocess import check_output
import re
import os
import platform


class GmxTopologyOutput(TopologyOutput):
    """
    Automatically sets a topology-file path with extension ``.top`` based on
    the directory path ``workdir``, the system name ``name``, the grid-shift
    iteration ``grid`` and the grid-point linear index ``state``.
    """
    def __init__(self,
                 workdir,
                 name,
                 grid,
                 state):
        self.fn = os.path.join(workdir, f"{name}_{grid}_{state}.top")

    @property
    def path(self):
        """
        The path of the topology file.
        """
        return os.path.abspath(self.fn)



def gmx_topo_out_creator(workdir: str,
                         name: str,
                         grid: int,
                         state: int,
                         attrs: CustomizableAttributesMixin.InputParameters):
    return GmxTopologyOutput(workdir, name, grid, state)


# ----------------------------------------------------------------------
# Writer
# ----------------------------------------------------------------------

def gmx_find_recursive_in_dir(dirpath, fn):
    in_directory = glob(os.path.join(dirpath, "**", fn), recursive=True)
    if len(in_directory) == 1:
        sub_fn = os.path.abspath(in_directory[0])
        print(f"Include: {fn} matched {sub_fn}.")
        return sub_fn
    elif len(in_directory) > 1:
        raise RuntimeError(f"Include: several matches for {fn} in current"
                           f" directory: {in_directory}.")
    else:
        return None


def gmx_find_include_file(fn, original_fn):
    # try to build path of the included file based on the location of the file
    # that includes it
    real_fn = os.path.join(
        os.path.dirname(os.path.abspath(original_fn)),
        fn)
    # the path is specified directly
    if os.path.isfile(real_fn):
        print(f"Include: {fn} matched {real_fn}.")
        return real_fn
    # search in subdirectories of the current directory
    sub_fn = gmx_find_recursive_in_dir(os.getcwd(), fn)
    if sub_fn:
        return sub_fn
    # search in GROMACS library as inferred from gmx
    if platform.system() == "Windows":
        which_command = "where"
    else:
        which_command = "which"
    bin_dir = os.path.dirname(
        check_output([which_command, "gmx"]).decode('utf-8'))
    lib_dir = os.path.join(bin_dir, "..")
    sub_fn = gmx_find_recursive_in_dir(lib_dir, fn)
    if sub_fn:
        return sub_fn
    # search in path inferred by GMXLIB environment variable
    try:
        lib_dir = os.environ["GMXLIB"]
        sub_fn = gmx_find_recursive_in_dir(lib_dir, fn)
        if sub_fn:
            return sub_fn
        raise RuntimeError(f"Could not find included file: {fn} in"
                           f" {original_fn}.")
    except KeyError:
        raise RuntimeError(f"Could not find included file: {fn} in"
                           f" {original_fn}.")


def gmx_top_ff_to_stream(fn: str, stream: StringIO):
    fp = open(fn, 'r')
    for line in fp:
        # If the line contains an include directive, recursively apply this
        # function to the file specified by the include.
        m = re.match(r'^#include *"{,1}(\S+?)"{,1}\s*$', line)
        if m:
            i_fn = gmx_find_include_file(m.group(1), fn)
            gmx_top_ff_to_stream(i_fn, stream)
        else:
            stream.write(line)
    fp.close()


def gmx_read_top_ff(fn: str) -> StringIO:
    stream = StringIO()
    gmx_top_ff_to_stream(fn, stream)
    return stream


def gmx_macro_replace(param: InteractionParameter,
                      istream: StringIO,
                      ostream: StringIO):
    for line in istream:
        ostream.write(line.replace(param.name, str(param.value)))


def gmx_replace_V_W(param,
                    V_type,
                    W_type,
                    old_value_V,
                    old_value_W,
                    ostream):
    at_i = param.particles.name_i
    at_j = param.particles.name_j
    if param.type == V_type:
        ostream.write(
            f"{at_i:6}{at_j:6} 1 {param.value:14.6e}{old_value_W:14.6e} "
            f"; GMAK\n")
    elif param.type == W_type:
        ostream.write(
            f"{at_i:6}{at_j:6} 1 {old_value_V:14.6e}{param.value:14.6e} "
            f"; GMAK\n")


def gmx_replace_pairblock(param: InteractionParameter,
                          istream: StringIO,
                          ostream: StringIO,
                          blockname: str,
                          V_type,
                          W_type):
    in_block = False
    has_pair = False
    for line in istream:
        # if comment, write and continue
        if re.match(r"^\s*;.*$", line):
            ostream.write(line)
            continue
        # if inside block, process and continue
        if in_block:
            tokens = line.split()
            if len(tokens) > 0:
                at_i = tokens[0]
                at_j = tokens[1]
                old_V = float(tokens[3])
                old_W = float(tokens[4])
                pair = InteractionPair(at_i, at_j)
                if param.particles == pair:
                    has_pair = True
                    gmx_replace_V_W(param, V_type, W_type, old_V, old_W, ostream)
                    in_block = False
                    continue
        # if not in block, test for start of block
        m = re.match(rf"^\[ {blockname} \]", line)
        if m:
            in_block = True
        ostream.write(line)
    if not has_pair:
        raise ValueError(f"We couldn't find the parameter {param.name} in the"
                         " force-field files. Please, make sure that they"
                         " include a [ {blockname} ] directive under which the "
                         f"desired pair ({param.particles.name_i},"
                         f"{param.particles.name_j}) is listed.")


def gmx_14_lj_replace(param: InteractionParameter,
                      istream: StringIO,
                      ostream: StringIO):
    gmx_replace_pairblock(param, istream, ostream, "pairtypes",
                          InteractionParameterType.LJ_14_V,
                          InteractionParameterType.LJ_14_W)


def gmx_lj_replace_pair(param: InteractionParameter,
                        istream: StringIO,
                        ostream: StringIO):
    gmx_replace_pairblock(param, istream, ostream, "nonbond_params",
                          InteractionParameterType.LJ_V,
                          InteractionParameterType.LJ_W)


def gmx_lj_replace_atom(param: InteractionParameter,
                        istream: StringIO,
                        ostream: StringIO,
                        combrule: CombinationRule):
    in_block = False
    old_atomic_param_value = None
    for line in istream:
        # if comment, write and continue
        if re.match(r'^\s*;.*$', line):
            ostream.write(line)
            continue
        if in_block:
            m = re.match(r'^\[ (\S+) \]', line)
            if m:
                ostream.write(line)
                last_block = m.group(1)
                in_block = False
                break
            line = line.split(';')[0]
            tokens = line.split()
            if len(tokens) > 0:
                at = tokens[0]
                V_loc = -2
                W_loc = -1
                V = float(tokens[V_loc])
                W = float(tokens[W_loc])
                if InteractionAtom(at) == param.particles:
                    if param.type == InteractionParameterType.LJ_V:
                        replace_loc = V_loc
                        old_atomic_param_value = V
                    elif param.type == InteractionParameterType.LJ_W:
                        replace_loc = W_loc
                        old_atomic_param_value = W
                    tokens[replace_loc] = f"{float(param.value):14.6e}"
                    tokens.append("; GMAK\n")
                    ostream.write("    ".join(tokens))
                    continue
        m = re.match(r'^\[ atomtypes \]', line)
        if m:
            in_block = True
        ostream.write(line)

    # read until [ nonbond_params ]
    if last_block != 'nonbond_params':
        for line in istream:
            ostream.write(line)
            if re.match(r'^\[ nonbond_params \]'):
                break

    # now, replace desired pairs
    for line in istream:
        # if comment, write and continue
        if re.match(r'^\s*;.*$', line):
            ostream.write(line)
            continue
        # if a new block, write and break
        if re.match(r'^\[ \S+ \]', line):
            ostream.write(line)
            break
        # in other cases, process the line
        tokens = line.split()
        if len(tokens) > 0:
            at_i = tokens[0]
            at_j = tokens[1]
            line_pair = InteractionPair(at_i, at_j)
            if line_pair.derives_from_atom(param.particles):
                # redefine parameter value
                if param.type == InteractionParameterType.LJ_V:
                    replace_loc = 3
                elif param.type == InteractionParameterType.LJ_W:
                    replace_loc = 4
                new_value = combrule.rescale_param_value(param,
                                                         float(tokens[replace_loc]),
                                                         old_atomic_param_value)
                tokens[replace_loc] = f"{new_value:14.6e}"
                tokens.append("; GMAK\n")
                # write line
                ostream.write("    ".join(tokens))
            else:
                # write line without modifications
                ostream.write(line)
        else:
            ostream.write(line)

    # write the rest
    for line in istream:
        ostream.write(line)


def gmx_lj_replace(param: InteractionParameter,
                   istream: StringIO,
                   ostream: StringIO,
                   combrule: CombinationRule):
    if isinstance(param.particles, InteractionPair):
        gmx_lj_replace_pair(param, istream, ostream)
    elif isinstance(param.particles, InteractionAtom):
        gmx_lj_replace_atom(param, istream, ostream, combrule)
    else:
        raise ValueError("LJ particles should be an atom or pair.")


class GmxCustomReplace:

    writer = None

    @classmethod
    def add_gmx_custom_parameter_writer(cls, writer):
        cls.writer = writer

    def __call__(self, p, s, o):
        if self.writer is None:
            raise NotImplementedError(f"Gmx custom type was not implemented and {p.name} is a custom parameter.")
        else:
            return self.writer(p, s, o)

def gmx_set_param(param: InteractionParameter, stream: StringIO, combrule:
                  CombinationRule) -> StringIO:
    ostream = StringIO()
    stream.seek(0)
    if param.type == InteractionParameterType.MacroParameter:
        gmx_macro_replace(param, stream, ostream)
    if param.type == InteractionParameterType.LJ_14_V:
        gmx_14_lj_replace(param, stream, ostream)
    if param.type == InteractionParameterType.LJ_14_W:
        gmx_14_lj_replace(param, stream, ostream)
    if param.type == InteractionParameterType.LJ_V:
        gmx_lj_replace(param, stream, ostream, combrule)
    if param.type == InteractionParameterType.LJ_W:
        gmx_lj_replace(param, stream, ostream, combrule)
    if param.type == InteractionParameterType.CustomParameter:
        GmxCustomReplace()(param, stream, ostream)
    return ostream


def gmx_read_combrule(stream: StringIO) -> CombinationRule:
    """Reads the combination rule from stream and seeks it to zero."""
    stream.seek(0)
    for line in stream:
        if re.match(r'^\[ defaults \]', line):
            break
    for line in stream:
        if re.match(r'^\s*;', line):
            continue
        else:
            combrule_code = int(line.split()[1])
            break
    stream.seek(0)
    return CombinationRule(combrule_code)


def gmx_topo_out_writer(params: List[InteractionParameter],
                        topo_out: GmxTopologyOutput,
                        attrs: CustomizableAttributesMixin.InputParameters):
    stream = gmx_read_top_ff(attrs.template)
    combrule = gmx_read_combrule(stream)
    for i, param in enumerate(params):
        replaced_stream = gmx_set_param(param, stream, combrule)
        stream.close()
        stream = StringIO(replaced_stream.getvalue())
        if i < len(params) - 1:
            replaced_stream.close()
    outstream = open(topo_out.path, 'w')
    outstream.write(replaced_stream.getvalue())
    outstream.close()
    replaced_stream.close()


class GmxSystem(System):

    def __init__(self, name: str):
        super().__init__(name, "gmx", gmx_topo_out_creator,
                         gmx_topo_out_writer)
