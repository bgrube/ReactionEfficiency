#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


def filterLinesInFile(
  file,
  startMarker,
  nmbLines = 1
):
  """
  generator that returns all lines from the given file handle, except
  those text blocks with nmbLines lines that start with startMarker
  """
  # see https://codereview.stackexchange.com/a/109806
  lines = iter(file)
  try:
    while True:
      line = next(lines)
      if startMarker in line:
        # discard nmbLines lines including the one with startMarker
        for _ in range(nmbLines - 1):
          next(lines)
      else:
        yield line
  except StopIteration:
    return


def removeFromFile(
  fileName,
  startMarker,
  nmbLines = 1
):
  """
  removes all text blocks with nmbLines lines that start with
  startMarker from given file
  """
  with open(fileName, 'r+') as file:
    filteredLines = list(filterLinesInFile(file, startMarker, nmbLines))
    file.seek(0)
    file.writelines(filteredLines)
    file.truncate()


if __name__ == "__main__":
  removeFromFile("./fitMissingMassSquared.log", "Error in <TBranch::WriteBasketImpl>: basket's WriteBuffer failed.", 13)
