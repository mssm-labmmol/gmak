#!/usr/bin/env python3

from copy import deepcopy
import re
import os
from glob import glob

######################################################################
# Example
######################################################################

class BlockExampleData:
    def __init__(self, block_contents):
        self._data = block_contents

    def write(self, stream, formatter):
        stream.write(self.print(formatter))

    def print(self, formatter=None):
        if formatter == "code":
            _commentless = []
            for line in self._data:
                if not re.match('^\s*#.*$', line):
                    _commentless.append(line)
            return "".join(_commentless)
        elif formatter == "comments":
            _comment = []
            count = 0
            for line in self._data:
                m = re.match('^\s*#\s*(.*)$', line)
                if m:
                    _comment.append(m.group(1) + "\n")
                    count += 1
            return "".join(_comment)
        else:
            return formatter.print_example_block(self)

class ExampleData:
    def __init__(self, data):
        self._data = {}
        for k, v in data.items():
            self._data[k] = BlockExampleData(v)

    def __getitem__(self, blockname):
        return self._data[blockname]

######################################################################
# Input
######################################################################

class InputData:
    pass

    def write(self, stream, formatter):
        stream.write(self.print(formatter))


######################################################################
# Option Types
######################################################################

class OptionType:
    pass

class NumericalOption(OptionType):
    def print(self, formatter):
        return formatter.print_numerical(self)

class ExprOption(OptionType):
    def print(self, formatter):
        return formatter.print_expr(self)

class PathOption(OptionType):
    def print(self, formatter):
        return formatter.print_path(self)

class ObjectNameOption(OptionType):
    def print(self, formatter):
        return formatter.print_objname(self)

class StringOption(OptionType):
    def print(self, formatter):
        return formatter.print_string(self)

class ListOption(OptionType):
    def __init__(self, list_of_otype):
        self.subtypes = list_of_otype

    def print(self, formatter):
        return formatter.print_list(self)

class UnionOption(OptionType):
    def __init__(self, possibilities):
        self.subtypes = possibilities

    def print(self, formatter):
        return formatter.print_union(self)

######################################################################
# Usage Types
######################################################################

class UsageType:
    pass

class RequiredUsage(UsageType):
    def __init__(self):
        self.msg = " "

class OptionalUsage(UsageType):
    def __init__(self):
        self.msg = "*(Optional)* "

class CustomUsage(UsageType):
    def __init__(self, msg):
        self.msg = f"*({msg})* "

######################################################################
# Options read from file
######################################################################

# e.g. workdir
class RawOption(InputData):
    def __init__(self, name, descr, otype, utype, remarks, mergeable):
        self.name = name
        self.descr = descr
        self.otype = otype
        self.usage = utype
        self.remarks = remarks
        self.mergeable = mergeable

    def print(self, formatter):
        return formatter.print_raw(self)

    def is_optional(self):
        return self.usage.is_optional()


class BlockData(InputData):
    def __init__(self, name, descr, depends, options, uniqueness):
        self.name = name
        self.descr = descr
        self.depends = depends
        self.uniqueness = uniqueness
        # options is a list of RawOption
        self.options = options

    def is_unique(self):
        return self.uniqueness

    def print(self, formatter):
        return formatter.print_block(self)

######################################################################
# Parsers
######################################################################

class TypeFactory:
    @classmethod
    def create(cls, string):
        if string == 'expr':
            return ExprOption()
        if string == 'numerical':
            return NumericalOption()
        elif string == 'path':
            return PathOption()
        elif string == 'string':
            return StringOption()
        elif string == 'objname':
            return ObjectNameOption()
        elif string.startswith("list"):
            subtypes_s = string.split()[1:]
            subtypes = [cls.create(s) for s in subtypes_s]
            return ListOption(subtypes)
        elif string.startswith("union"):
            subtypes_s = string.replace("union", "").split('|')
            subtypes = [cls.create(s.strip()) for s in subtypes_s]
            return UnionOption(subtypes)

class UsageFactory:
    @classmethod
    def create(cls, string):
        if string.startswith('custom'):
            m = re.match('custom\s+(.*)$', string)
            msg = m.group(1)
            return CustomUsage(msg)
        if string == 'required':
            return RequiredUsage()
        elif string == 'optional':
            return OptionalUsage()
        else:
            raise ValueError

class RawParser:
    @classmethod
    def parse(cls, stream):
        name = stream.readline().strip()
        if name == '':
            # this signals the end
            return None
        usage = UsageFactory.create(stream.readline().strip())
        type = TypeFactory.create(stream.readline().strip())
        description = stream.readline().strip()
        remarks = []
        for line in stream:
            if line == '\n':
                break
            else:
                remarks.append(line.strip())
        if "mergeable" in remarks:
            mergeable = True
            remarks.remove("mergeable")
        else:
            mergeable = False
        return RawOption(name, description, type, usage, remarks, mergeable)


class BlockParser:
    @classmethod
    def parse(cls, stream, name):
        descr = stream.readline().strip()
        unique = stream.readline().strip()
        if unique == "unique":
            uniqueness = True
        else:
            uniqueness = False
        depends = stream.readline().strip().split(',')
        depends = [d.strip() for d in depends]
        options = []
        while True:
            opt = RawParser.parse(stream)
            if opt is None:
                break
            else:
                options.append(opt)
        return BlockData(name, descr, depends, options, uniqueness)


class FileParser:
    @classmethod
    def parse(cls, fn, isblock=True):
        block = os.path.basename(fn)
        fp = open(fn, 'r')
        if isblock:
            out = BlockParser.parse(fp, block)
        else:
            out = RawParser.parse(fp)
        fp.close()
        return out


class ExampleParser:
    @classmethod
    def parse(cls, fn):
        data = {}
        block = ""
        in_block = False
        fp = open(fn, 'r')
        for line in fp:
            m = re.match(r'^\s*\$(\S+)$', line)
            if m:
                in_block = True
                _block = m.group(1)
                if _block != "end":
                    block = _block
            if m and _block != "end" and block != "":
                data[block] = []
            if in_block:
                data[block].append(line)
            if line.strip() == '$end':
                in_block = False
        fp.close()
        print(data)
        return ExampleData(data)

######################################################################
# Formatters
######################################################################

class TypeFormatter:

    @classmethod
    def print_expr(cls, input):
        return "expression"

    @classmethod
    def print_objname(cls, input):
        return "object reference"

    @classmethod
    def print_numerical(cls, input):
        return "numerical"

    @classmethod
    def print_path(cls, input):
        return "path"

    @classmethod
    def print_string(cls, input):
        return "string"

    @classmethod
    def print_list(cls, input):
        out = "List["
        subs = [x.print(cls) for x in input.subtypes]
        out += ", ".join(subs)
        out += "]"
        return out

    @classmethod
    def print_union(cls, input):
        out = " OR ".join([x.print(cls) for x in input.subtypes])
        return out


class TOCFormatter:

    @classmethod
    def print_raw(cls, input):
        if input.mergeable:
            merge = "*(Mergeable)* "
        else:
            merge = ""
        descr = f"{merge}{input.usage.msg}{input.descr}"
        out = f"""
   * - {input.name}
     - {input.otype.print(TypeFormatter)}
     - {descr}
     - """
        for rm in input.remarks:
            out += f"""{rm} """
        return out

    @classmethod
    def print_block(cls, input):
        if input.depends[0].lower() == "none":
            depends_string = ""
        else:
            depends_string = ", ".join([":doc:`$%s <blocks/%s>`" % (x, x) for x in input.depends])
            depends_string = "Depends on " + depends_string + "."
        return f"""
   * - :doc:`{input.name} </usage/blocks/{input.name}>`
     - {input.descr}
     - {depends_string}
"""

class TOCSectionizer:
    @staticmethod
    def get_section(opt):
        # section is determined by remarks of the form
        # "For type = XXX."
        section = "Block parameters"
        for remark in opt.remarks:
            m = re.match(r'For type = (\S+).', remark)
            if m:
                section = f"Type ``{m.group(1)}``"
        return section

    @classmethod
    def get_sections(cls, opts):
        sections = [cls.get_section(opt) for opt in opts]
        # make sections unique
        unique_sections = []
        for section in sections:
            if section not in unique_sections:
                unique_sections.append(section)
        return unique_sections

    @classmethod
    def get_section_opts(cls, section, opts):
        section_opts = []
        for opt in opts:
            if cls.get_section(opt) == section:
                section_opts.append(opt)
        return section_opts

    @classmethod
    def get_sections_opts(cls, opts):
        sections_opts = []
        sections = cls.get_sections(opts)
        for section in sections:
            sections_opts.append((section, cls.get_section_opts(section, opts)))
        # remove remarks
        for section, opts in sections_opts:
            for opt in opts:
                for remark in deepcopy(opt.remarks):
                    m = re.match(r'For type = (\S+).', remark)
                    if m:
                       opt.remarks.remove(remark)
        return sections_opts


class BlockFileFormatter:

    @staticmethod
    def print_section(title, level):
        title_len = len(title)
        level_markers = ["#", "=", "-", "~"]
        output = "\n"
        if level == 0:
            output += level_markers[0] * title_len
            output += "\n"
        output += title
        output += "\n"
        output += level_markers[level] * title_len
        output += "\n"
        return output

    @staticmethod
    def print_custom_attributes_note():
        output = ""
        output += f"""
.. note:: Parameters that are not listed above can also be supplied.
   They are not recognized by the program in any special way, but are
   parsed and made available in the :doc:`/usage/customization_api`,
   together with all the other block parameters, as an
   :py:class:`~gmak.custom_attributes.CustomizableAttributesMixin.InputParameters`
   object.
"""
        return output

    @staticmethod
    def print_descr(input_data):
        return input_data.descr + "\n"

    @staticmethod
    def print_uniqueness(input_data):
        if input_data.is_unique():
            return "This block can appear only once in the input file.\n"
        else:
            return "This block can appear multiple times in the input file.\n"

    @staticmethod
    def print_dependencies(input_data):
        ndep = len(input_data.depends)
        if ndep > 0:
            if input_data.depends[0] == 'none':
                return ""
            deps = ", ".join([
                ":doc:`$%s <%s>`" % (x, x)
                for x in input_data.depends])
            if ndep == 1:
                return f"It intrinsically depends on the {deps} block."
            else:
                return f"It intrinsically depends on the {deps} blocks."
        else:
            return ""

    @staticmethod
    def print_toc(options):
        output = """
 .. list-table::
   :header-rows: 1
   :widths: 10 10 10 10
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks
"""
        for opt in options:
            output += opt.print(TOCFormatter)
        return output

    @classmethod
    def print_example_block(cls, block_example):
        output = cls.print_section("Example", 1)
        output += "\n"
        output += "See :doc:`/examples/tutorial` for a commented example."
        # output += ".. code-block:: gmi\n\n"
        # output += block_example.print(formatter="code")
        # output += "\n\n"
        # output += block_example.print(formatter="comments")
        return output

    @classmethod
    def print_block(cls, input_data):
        output = ""
        output += cls.print_section(f"The ``${input_data.name}`` block", 0)
        output += "\n"
        output += cls.print_descr(input_data)
        output += cls.print_uniqueness(input_data)
        output += cls.print_dependencies(input_data)
        output += "\n\nThe input parameters that are set in this block are listed below.\n"
        secs_opts = TOCSectionizer.get_sections_opts(input_data.options)
        for sec, opt in secs_opts:
            output += cls.print_section(sec, 1)
            output += cls.print_toc(opt)
            output += "\n"
        output += cls.print_custom_attributes_note()
        return output


class MergeableLineFormatter:
    '''Formats lines for the Mergeables table. Assumes it is mergeable.'''
    def __init__(self, blockname):
        self.blockname = blockname

    def print_raw(self, input):
        remark = " "
        for rem in input.remarks:
            m = re.match(r'For type = (\S+).', rem)
            if m:
                remark = f"For type ``{m.group(1)}``"
        return f"""
   * - :doc:`{self.blockname} <blocks/{self.blockname}>`
     - {input.name}
     - {remark}
"""


class MergeableBlockFormatter:
    '''Formats all lines of a block for the Mergeables table.'''
    @classmethod
    def print_block(cls, input):
        out = ""
        for opt in input.options:
            if opt.mergeable:
                out += opt.print(MergeableLineFormatter(input.name))
        return out


class MergeableTableFormatter:
    '''Formats the entire table.'''
    @classmethod
    def write(cls, fn, blocks):
        fp = open(fn, 'w')
        fp.write(f"""
.. list-table::
   :header-rows: 1
   :widths: 10 10 10
   :align: center

   * - Block
     - Parameter
     - Remarks""")
        for block in blocks:
            block.write(fp, MergeableBlockFormatter)
        fp.close()


if __name__ == '__main__':
    blocks = []
    raws = []
    for fn in glob("blocks/*"):
        blocks.append(FileParser.parse(fn))
    for fn in glob("raw/*"):
        raws.append(FileParser.parse(fn, isblock=False))

    examples = ExampleParser.parse("example")

    # now write mergeables
    MergeableTableFormatter.write("mergeables.rst", blocks)

    fp = open('toc.rst', 'w')
    fp.write("""
 .. list-table::
   :header-rows: 1
   :widths: auto
   :align: center

   * - Parameter
     - Value Type
     - Description
     - Remarks
""")
    for raw in raws:
        raw.write(fp, TOCFormatter)
    fp.close()

    fp = open(f"blocks.rst", 'w')
    fp.write("""
 .. list-table::
   :header-rows: 1
   :align: center

   * - Block name
     - Description
     - Remarks
""")
    for block in blocks:
        block.write(fp, TOCFormatter)
    fp.close()

    # create the toctree
    fp = open("toctree.rst", "w")
    fp.write(f"""
.. toctree::
   :hidden:
""")
    # only blocks have pages
    for block in blocks:
        fp.write(f"""
   blocks/{block.name}
""")
    fp.close()

    # now create files for each block
    for block in blocks:
        fp = open(f"{block.name}.rst", "w")
        block.write(fp, BlockFileFormatter)
        examples[block.name].write(fp, BlockFileFormatter)
        fp.close()
