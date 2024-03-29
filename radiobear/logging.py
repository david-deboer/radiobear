# -*- mode: python; coding: utf-8 -*-
"""Log function."""
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.
import os


def setup(log):
    """Set the log up, if not already."""
    if isinstance(log, LogIt):
        return log
    return LogIt(log)


class LogIt:
    """Log class."""

    def __init__(self, log):
        """Initialize instance."""
        if isinstance(log, str):
            self.logfile = log
            log_directory = os.path.dirname(log)
            self.log_directory = log_directory
            if not os.path.isdir(log_directory):
                os.mkdir(log_directory)
            self.fp = open(log, 'a')
        else:
            self.fp = log

    def add(self, msg, printOut=True):
        """Add to log."""
        if self.fp is not None:
            self.fp.write(msg + '\n')
            if printOut:
                print(msg)

    def close(self):
        """Close log."""
        if self.fp is not None:
            self.fp.close()
