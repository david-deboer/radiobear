# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.
import six
import os


def setup(log):
    if isinstance(log, LogIt):
        return log
    return LogIt(log)


class LogIt:
    def __init__(self, log):
        if isinstance(log, six.string_types):
            self.logfile = log
            log_directory = os.path.dirname(log)
            self.log_directory = log_directory
            if not os.path.isdir(log_directory):
                os.mkdir(log_directory)
            self.fp = open(log, 'a')
        else:
            self.fp = log

    def add(self, msg, printOut=True):
        if self.fp is not None:
            self.fp.write(msg + '\n')
            if printOut:
                print(msg)

    def close(self):
        if self.fp is not None:
            self.fp.close()