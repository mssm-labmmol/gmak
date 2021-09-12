from datetime import datetime
import gmak.runcmd as runcmd
import os

class Logger (object):

    def __init__ (self, fn, indentString='  '):
        self.fn = fn
        self.fp = open(fn, 'w')
        self.indentation = indentString
        self._currentIndentation = 0
        print(f"Starting log file {fn}.")

    def __del__(self):
        self.fp.close()

    @staticmethod
    def getTimestamp():
        return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def indent(self, number=1):
        self._currentIndentation += number

    def unindent(self, number=1):
        self._currentIndentation -= number

    def putMessage(self, msg, dated=False):
        self.fp.write(self._currentIndentation * self.indentation)
        self.fp.write(msg)
        self.fp.write(' @ ' + Logger.getTimestamp())
        self.fp.write('\n')
        self.fp.flush()

class EmptyLogger (Logger):

    def __init__ (self):
        pass

    def putMessage(self, msg, dated=False):
        pass

    def __del__(self):
        pass

globalLogger = EmptyLogger() # this will be altered during program execution
