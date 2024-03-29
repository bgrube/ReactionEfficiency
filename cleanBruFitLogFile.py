#!/usr/bin/env python3


import argparse
import functools
import glob
import sys


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def filterLinesInFile(
  file,
  markers,
  removeTrailingEmptyLines = True
):
  """
  generator that returns all lines from the given file handle, except
  those text blocks that start with markers given by the keys of the
  `markers` dict; the dict values give the block length in number of
  lines (including the marker line)
  """
  # see https://codereview.stackexchange.com/a/109806
  lines = iter(file)
  try:
    foundMarker = False
    while True:
      line = next(lines)
      # discard any empty lines trailing the previously found marker
      while foundMarker and removeTrailingEmptyLines and line == "\n":
        line = next(lines)
      foundMarker = False
      for marker, nmbLines in markers.items():
        if marker in line:
          foundMarker = True
          # discard nmbLines lines including the one with startMarker
          for _ in range(nmbLines - 1):
            next(lines)
          # if removeTrailingEmptyLines:
          #   while True:
          #     line = next(lines)
          #     if line != "\n":
          #       yield line
          #       break
          break
      if not foundMarker:
        yield line
  except StopIteration:
    return


def removeFromFiles(
  fileNamePattern,
  markers
):
  """
  removes all text blocks from file(s) matching `fileNamePattern`,
  which start with markers given by the keys of the `markers` dict;
  the dict values give the block length in number of lines (including
  the marker line)
  """
  fileNames = sorted(glob.glob(fileNamePattern))
  for fileName in fileNames:
    with open(fileName, 'r+') as file:
      print(f"Cleaning file '{fileName}'")
      filteredLines = list(filterLinesInFile(file, markers))
      file.seek(0)
      file.writelines(filteredLines)
      file.truncate()


BRANCH_ERROR_MARKERS = {
  "Error in <TBranch::WriteBasketImpl>: basket's WriteBuffer failed." : 1,
  "Error in <TBranch::TBranch::Fill>: Failed to write out basket." : 2,
  "Error in <TTree::Fill>: Failed filling" : 2,
  "This error is symptomatic of a Tree created as a memory-resident Tree" : 1,
  "Instead of doing:" : 1,
  "TTree *T = new TTree(...)" : 1,
  "TFile *f = new TFile(...)" : 1,
  "you should do:" : 1
}
CLING_ERROR_MARKERS = {
  "Error in <TCling::LoadPCM>: ROOT PCM" : 1,
  "Info in <TCling::LoadPCM>: In-memory ROOT PCM candidate" : 2
}


if __name__ == "__main__":
  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Cleans BruFit log files.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()

  markers = dict(BRANCH_ERROR_MARKERS, **CLING_ERROR_MARKERS)
  removeFromFiles(f"{args.outputDirName}/fitMissingMassSquared.log", markers)
  removeFromFiles(f"{args.outputDirName}/*/*/*.txt",                 markers)
