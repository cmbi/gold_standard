#!/usr/bin/python
import inspect


class CustomException(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

        caller = inspect.stack()[1][3]
        self.message = "{} raised an exception: {}".format(caller, message)
