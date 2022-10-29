#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement
# scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        seq1 = seq1[:align_length]
        seq2 = seq2[:align_length]

        table = {}
        backtrace = []

        finali = 0
        finalj = 0

        if banded:
            for i in range(MAXINDELS + 1):
                table[(i, 0)] = i * INDEL
            for j in range(MAXINDELS + 1):
                table[(0, j)] = j * INDEL
            for i in range(1, len(seq1) + 1):
                for j in range(1, len(seq2) + 1):
                    if abs(i - j) <= MAXINDELS:
                        diag = table[(i - 1, j - 1)] + SUB
                        if table.get((i - 1, j)) is None:
                            up = float('inf')
                        else:
                            up = table[(i - 1, j)] + INDEL
                        if table.get((i, j - 1)) is None:
                            left = float('inf')
                        else:
                            left = table[(i, j - 1)] + INDEL

                        if seq1[i - 1] == seq2[j - 1]:
                            table[(i, j)] = table[(i - 1, j - 1)] + MATCH
                            backtrace.append('diag')
                        elif diag <= up and diag <= left:
                            table[(i, j)] = diag
                            backtrace.append('diag')
                        elif up <= diag and up <= left:
                            table[(i, j)] = up
                            backtrace.append('up')
                        else:
                            table[(i, j)] = left
                            backtrace.append('left')
                        finali = i
                        finalj = j
                    else:
                        continue

        else:
            for i in range(len(seq1) + 1):
                table[(i, 0)] = i * INDEL
            for j in range(len(seq2) + 1):
                table[(0, j)] = j * INDEL
            for i in range(1, len(seq1) + 1):
                for j in range(1, len(seq2) + 1):
                    diag = table[(i - 1, j - 1)] + SUB
                    up = table[(i - 1, j)] + INDEL
                    left = table[(i, j - 1)] + INDEL
                    if seq1[i - 1] == seq2[j - 1]:
                        table[(i, j)] = table[(i - 1, j - 1)] + MATCH
                        backtrace.append('diag')
                    elif diag <= up and diag <= left:
                        table[(i, j)] = diag
                        backtrace.append('diag')
                    elif up <= diag and up <= left:
                        table[(i, j)] = up
                        backtrace.append('up')
                    else:
                        table[(i, j)] = left
                        backtrace.append('left')
                    finali = i
                    finalj = j

        ##################################################################################################
        # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
        score = table[finali, finalj]
        alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
            len(seq1), align_length, ',BANDED' if banded else '')
        alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
            len(seq2), align_length, ',BANDED' if banded else '')
        ###################################################################################################

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
