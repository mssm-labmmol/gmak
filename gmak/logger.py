from datetime import datetime
import gmak.runcmd as runcmd
import os

class Logger (object):
    def __init__ (self, fn, indentString='  '):
        self.fn = fn
        self.fp = open(fn, 'w')
        self.indentation = indentString
        self._currentIndentation = 0
        #print(f"Starting log file {fn}.")

    def __del__(self):
        self.fp.close()

    @staticmethod
    def getTimestamp():
        return datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def indent(self, number=1):
        """
        Indents the stream by a certain number of indentation strings.

        :param number: The number of indentation strings (default is 1).
        :type  number: int
        """
        self._currentIndentation += number

    def unindent(self, number=1):
        """
        Unindents the stream by a certain number of indentation strings.

        :param number: The number of indentation strings (default is 1).
        :type  number: int
        """
        self._currentIndentation -= number

    def putMessage(self, msg, dated=True):
        """
        Writes a message to the log-file stream with proper indentation.

        :param msg: The message.
        :type msg: str
        :param dated: If :py:obj:`True`, write a timestamp at the end of the message.
        :type dated: bool
        """
        self.fp.write(self._currentIndentation * self.indentation)
        self.fp.write(msg)
        if dated:
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

#: The global :py:class:`~gmak.logger.Logger` instance that has access
#: to the log-file stream.
globalLogger = EmptyLogger()
