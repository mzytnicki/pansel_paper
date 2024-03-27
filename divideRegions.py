#! /usr/bin/env python3

import sys
from collections import namedtuple, defaultdict


Bin = namedtuple("Bin", "chromosome start end f")


if __name__ == "__main__":
  inputFileName = sys.argv[1]
  binSize       = int(sys.argv[2])
  gapsFileName  = sys.argv[3]
  nGroups       = int(sys.argv[4])
  outputPrefix  = sys.argv[5]

  gaps = defaultdict(dict)
  with open(gapsFileName) as gapsFile:
    for line in gapsFile:
      line       = line.strip().split()
      chromosome = line[0]
      start      = int(line[1]) // binSize
      end        = (int(line[2]) - 1) // binSize
      for b in range(start, end + 1):
        gaps[chromosome][b] = True

  nBins         = 0
  nSelectedBins = 0
  bins          = []
  with open(inputFileName) as inputFile:
    for line in inputFile:
      line       = line.strip().split()
      chromosome = line[0]
      start      = int(line[2])
      end        = int(line[3])
      f          = float(line[4]) / float(line[5])
      inGap      = False
      nBins     += 1
      if end - start <= 1.5 * binSize:
        for b in range(start // binSize, end // binSize):
          if chromosome in gaps and b in gaps[chromosome]:
            inGap = True
        if not inGap:
          nSelectedBins += 1
          bins.append(Bin(chromosome, start, end, f))
  print(f"{nBins} seen, {nSelectedBins} selected.")

  nBins      = len(bins)
  bins       = sorted(bins, key = lambda b: b.f)
  thresholds = [int(round(nBins * (float(i) / nGroups))) for i in range(1, nGroups)]
  breakId    = 0
  idStart    = 0
  fStart     = bins[0].f
  outputFile = open(f"{outputPrefix}_0.bed", 'w')
  for i in range(nBins):
    if (breakId < len(thresholds)) and (i == thresholds[breakId]):
      print(f"Break #{breakId}: ids: {idStart} - {i-1}, f: {fStart} - {bins[i-1].f}")
      breakId   += 1
      nElements = 0
      fStart    = bins[i-1].f
      idStart   = i-1
      outputFile.close()
      outputFile = open(f"{outputPrefix}_{breakId}.bed", 'w')
    outputFile.write(f"{bins[i].chromosome}\t{bins[i].start-1}\t{bins[i].end}\tbin_{i}_{breakId}\t{bins[i].f}\t+\n")                                                                                               
  print(f"Break #{breakId}: ids: {idStart} - {i-1}, f: {fStart} - {bins[i-1].f}")
