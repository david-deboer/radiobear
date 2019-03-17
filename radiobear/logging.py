# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.


class Log:
    def __init__(log):
        if isinstance(log, six.string_types):
            self.logfile = log
            log_directory = os.path.dirname(log)
            if not os.path.isdir(log_directory):
                os.mkdir(log_directory)
            self.fp = open(log, 'a')
        else:
            self.fp = log

    def log(msg, printOut=True):
        if self.fp is not None:
            self.fp.write(msg + '\n')
            if printOut:
                print(msg)

    def close():
        if self.fp is not None:
            self.fp.close()
