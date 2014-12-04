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
import logging

from pbcore.io import BasH5Reader, BaxH5Reader
from pbcore.io.FastaIO import *
import pbbarcodejr.SWaligner as Aligner
import numpy as n

from pbcore.io.BarcodeH5Reader import BARCODE_DELIMITER

__RC_MAP__ = dict(zip('ACGTacgt-N','TGCAtgca-N'))

class BarcodeScorer(object):
    def __init__(self, basH5, barcodeFasta,
                 adapterSidePad = 0,
                 insertSidePad = 4,
                 scoreMode = 'symmetric',
                 maxHits = 10):
        """A BarcodeScorer object scores ZMWs and produces summaries
        of the scores. Various parameters control the behavior of the
        object, specifically the padding allows the user to add a
        little extra on each side of the adapter find for safety. The
        most relevant parameter is the scoreMode which dictates how
        the barcodes are scored, either paired or symmetric."""

        self.basH5 = basH5
        self.barcodeFasta = list(barcodeFasta)
        self.aligner = Aligner.SWaligner()
        self.adapterSidePad = adapterSidePad
        self.insertSidePad = insertSidePad
        self.maxHits = maxHits

        if scoreMode not in ['symmetric', 'paired']:
            raise Exception("scoreMode must either be symmetric or paired")
        self._scoreMode = scoreMode

        # Test the barcodes to ensure they all have the same length barcodes
        self.barcodeLength = n.unique(map(lambda x : len(x.sequence), self.barcodeFasta))
        if len(self.barcodeLength) > 1:
            raise Exception("Currently, all barcodes must be the same length.")
        else:
            self.barcodeLength = int(self.barcodeLength)

        # Next we build the barcode scoring tool, using the forward sequence
        #   of the forward barcodes and the reverse-complement of the reverse
        #   barcodes
        self.barcodeSeqs = [bc.sequence.upper() if (i%2)==0 else
                            self._rc(bc.sequence.upper())
                            for i, bc in enumerate(self.barcodeFasta)]
        self.barcodeScorer = self.aligner.makeScorer(self.barcodeSeqs)

        logging.debug(("Constructed BarcodeScorer with scoreMode: %s," + \
                           "adapterSidePad: %d, insertSidePad: %d") \
                          % (scoreMode, adapterSidePad, insertSidePad))

    @property
    def movieName(self):
        return self.basH5.movieName

    def makeBCLabel(self, s1, s2):
        return BARCODE_DELIMITER.join((s1, s2))

    @property
    def barcodeLabels(self):
        """The barcode labels are function of the barcodeNames and the
        scoreMode, they represent the user-visible names."""
        if self.scoreMode == 'paired':
            return n.array([self.makeBCLabel(self.barcodeFasta[i].name,
                                             self.barcodeFasta[i+1].name) for i
                            in xrange(0, len(self.barcodeSeqs), 2)])
        else:
            return n.array([self.makeBCLabel(x.name, x.name) for x in self.barcodeFasta])

    @property
    def barcodeNames(self):
        """The barcode names are the FASTA names"""
        return n.array([x.name for x in self.barcodeFasta])

    @property
    def scoreMode(self):
        return self._scoreMode

    def _rc(self, s):
        """
        Reverse-Complement a DNA Sequence
        """
        return "".join([__RC_MAP__[c] for c in s[::-1]])

    def fromRange(self, zmw, rStart, rEnd):
        """
        Given a ZMW and the start/end positions of an adapter region in it,
        parse the left and right flanking regions if they exist
        """
        # Parse and reverse complement the left-side flanking region
        try:
            qSeqLeftRaw = zmw.read(rStart - (self.barcodeLength + self.insertSidePad),
                                   rStart + self.adapterSidePad).basecalls()
            qSeqLeft = self._rc(qSeqLeftRaw)
        except IndexError:
            qSeqLeft = None

        # Parse right-side flanking region
        try:
            qSeqRight = zmw.read(rEnd - self.adapterSidePad,
                                 rEnd + self.barcodeLength +
                                 self.insertSidePad).basecalls()
        except IndexError:
            qSeqRight = None

        # Return the results
        return (qSeqLeft, qSeqRight)

    def _flankingSeqs(self, zmw):
        """
        Extract flanking sequences around each adapter for the first X
        adapters for scoring
        """
        adapterRegions = zmw.adapterRegions
        if len(adapterRegions) > self.maxHits:
            adapterRegions = adapterRegions[0:self.maxHits]

        return [self.fromRange(zmw, start, end) for (start, end) in adapterRegions]

    def scoreZmw(self, zmw):
        adapters = self._flankingSeqs(zmw)
        adapterScores = [[]]*len(adapters)
        barcodeScores = n.zeros(len(self.barcodeSeqs))

        for i,adapter in enumerate(adapters):
            fscores  = self.barcodeScorer(adapter[0])
            rscores  = self.barcodeScorer(adapter[1])

            scored = 2.0 if adapter[0] and adapter[1] \
                else 1.0 if adapter[0] or  adapter[1] \
                else 0

            # An adapter score is the average barcode score for
            # each barcode -- that way, you can compare across
            # adapters even if the different adapters have
            # different numbers of flanking sequence.
            if scored == 0:
                adapterScores[i] = barcodeScores
            else:
                adapterScores[i] = (fscores + rscores)/scored

        barcodeScores = reduce(lambda x, y: x + y, adapterScores) if adapterScores \
            else n.zeros(len(self.barcodeSeqs))

        return (zmw, len(adapters), barcodeScores, adapterScores)

    def displaySymmetric(self, o):
        zmw       = o[0]
        numAdp    = o[1]
        bcScores  = o[2]
        sortedBc  = n.argsort(-bcScores)
        bestIdx, secondIdx = sortedBc[0], sortedBc[1]
        bestLabel, secondLabel = self.barcodeLabels[bestIdx], self.barcodeLabels[secondIdx]
        bestScore, secondScore = bcScores[bestIdx], bcScores[secondIdx]
        return "{0},{1},{2},{3},{4},{5},{6}".format(zmw.zmwName, numAdp, bestIdx, bestLabel, bestScore, secondIdx, secondLabel, secondScore)

    def displayPaired(self, o):
        zmw       = o[0]
        numAdp    = o[1]
        bcScores  = o[2]
        adpScores = o[3]
        if o[1] == 1:
            results = n.array([max(bcScores[i], bcScores[i + 1]) for i in \
                               xrange(0, len(self.barcodeSeqs), 2)])
            sortedBc = n.argsort(-results)
            sortedScores = results[sortedBc]
        else:
            # score the pairs by scoring the two alternate
            # ways they could have been put on the molecule. A
            # missed adapter will confuse this computation.
            scores  = o[3]
            results = n.zeros(len(self.barcodeSeqs)/2)
            for i in xrange(0, len(self.barcodeSeqs), 2):
                pths = [0,0]
                for j in xrange(0, len(scores)):
                    pths[j % 2] += scores[j][i]
                    pths[1 - j % 2] += scores[j][i + 1]
                results[i/2] = max(pths)
            sortedBc     = n.argsort(-results)
            sortedScores = results[sortedBc]
        bestIdx, secondIdx = sortedBc[0], sortedBc[1]
        bestLabel, secondLabel = self.barcodeLabels[bestIdx], self.barcodeLabels[secondIdx]
        bestScore, secondScore = sortedScores[0], sortedScores[1]
        return "{0},{1},{2},{3},{4},{5},{6}".format(zmw.zmwName, numAdp, bestIdx, bestLabel, bestScore, secondIdx, secondLabel, secondScore)

    def labelZmws(self, holeNumbers):
        """Return a list of LabeledZmws for input holeNumbers"""

        # Select the function for choosing the best barcode by
        if self.scoreMode == 'symmetric':
            display = self.displaySymmetric
        elif self.scoreMode == 'paired':
            display = self.displayPaired
        else:
            raise Exception("Unsupported scoring mode in BarcodeLabeler.py")

        for zmw in holeNumbers:
            scored = self.scoreZmw(self.basH5[zmw])
            if len(scored) > 2 and scored[1]:
                if self.scoreMode == 'symmetric':
                    yield self.displaySymmetric(scored)
                elif self.scoreMode == 'paired':
                    yield self.displayPaired(scored)
