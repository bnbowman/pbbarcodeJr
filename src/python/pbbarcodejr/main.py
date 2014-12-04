#!/usr/bin/env python
#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$
import os
import sys
import argparse
import logging
import tempfile
import shutil
import pkg_resources
import re
import subprocess
import random
import shutil

from multiprocessing import Pool

import h5py as h5
import numpy as n

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbcore.io import BaxH5Reader, BasH5Reader
from pbcore.io.BarcodeH5Reader import *
from pbcore.io import FastaReader, FastaRecord

from pbbarcodejr.BarcodeScorer import BarcodeScorer
from pbbarcodejr._version import __version__

SCORE_MODES    = ['symmetric', 'paired']

def makeBarcodeH5FromBasH5(basH5):
    """The workhorse function for creating a barcode H5 file from a
    base H5 file."""
    labeler = BarcodeScorer(basH5, FastaReader(runner.args.barcodeFile),
                            runner.args.adapterSidePad, runner.args.insertSidePad,
                            scoreMode = runner.args.scoreMode,
                            maxHits = runner.args.maxAdapters)
    if runner.args.nZmws < 0:
        zmws = basH5.sequencingZmws
    else:
        zmws = basH5.sequencingZmws[0:runner.args.nZmws]

    logging.debug("Labeling %d ZMWs from: %s" % (len(zmws), basH5.filename))
    labels = list(labeler.labelZmws(zmws))
    logging.debug("Labelled %d ZMWs from: %s" % (len(labels), basH5.filename))

    return labels

def mpWrapper(f):
    return makeBarcodeH5FromBasH5(BasH5Reader(f))

def makeBarcodeFofnFromBasFofn():
    inputFofn = runner.args.inputFile
    inFiles = open(inputFofn).read().splitlines()

    if not all(map(os.path.exists, inFiles)):
        raise IOError("All files in input.fofn must exist.")

    logging.debug("Using %d processes." % runner.args.nProcs)
    if runner.args.nProcs <= 1:
        lines = map(mpWrapper, inFiles)
    else:
        pool = Pool(runner.args.nProcs)
        labelLists = pool.map(mpWrapper, inFiles)

    for labelList in labelLists:
        for line in labelList:
            print line

class Pbbarcode(PBMultiToolRunner):
    def __init__(self):
        desc = ['Utilities for labeling and annoting reads with barcode information.']
        super(Pbbarcode, self).__init__('\n'.join(desc))
        subparsers = self.subParsers

        desc = ['Creates a barcode.h5 file from base h5 files.']
        parser_m = subparsers.add_parser('labelZmws', description = "\n".join(desc),
                                         help = 'Label zmws with barcode annotation',
                                         formatter_class = \
                                             argparse.ArgumentDefaultsHelpFormatter)
        parser_m.add_argument('--outDir',
                              help = 'Where to write the newly created barcode.h5 files.',
                              default = os.getcwd())
        parser_m.add_argument('--outFofn', help = 'Write to outFofn',
                              default = 'barcode.fofn')
        parser_m.add_argument('--adapterSidePad', help = 'Pad with adapterSidePad bases',
                              default = 0, type = int)
        parser_m.add_argument('--insertSidePad', help = 'Pad with insertSidePad bases',
                              default = 25, type = int)
        parser_m.add_argument('--scoreMode',
                              help = 'The mode in which the barcodes should be scored.',
                              choices = SCORE_MODES, default = 'symmetric', type = str)
        parser_m.add_argument('--maxAdapters', type = int, default = 20,
                              help = 'Only score the first maxAdapters')
        parser_m.add_argument('--nZmws', type = int, default = -1,
                              help = 'Use the first n ZMWs for testing')
        parser_m.add_argument('--nProcs', type = int, default = 8,
                              help = 'How many processes to use')
        parser_m.add_argument('barcodeFile', metavar = 'barcode.fasta',
                              help = 'Input barcode fasta file')
        parser_m.add_argument('inputFile', metavar = 'input.fofn',
                              help = 'Input base fofn')

    def getVersion(self):
        return  __version__

    def run(self):
        logging.debug("Arguments" + str(self.args))

        if self.args.subCommand == 'labelZmws':
            makeBarcodeFofnFromBasFofn()
        else:
            sys.exit(1)

runner = Pbbarcode()

def main():
    """The entry point for pbbarcode"""
    sys.exit(runner.start())
