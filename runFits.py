#!/usr/bin/env python3


from __future__ import annotations

from collections.abc import Generator
from contextlib import contextmanager
import functools
import os
import shutil
import subprocess
from typing import (
  Any,
  IO,
  TextIO,
)
import sys

import ROOT

from fitMissingMassSquared import fitMissingMassSquared


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def fileDescriptor(fileObjectOrDescriptor: IO[Any] | TextIO | int | None) -> int:
  """Returns the file descriptor for the file-like object of file descriptor `fileObjectOrDescriptor`"""
  # from https://stackoverflow.com/a/22434262
  # check if fileObjectOrDescriptor is a file-like object else assume it is a file descriptor
  fileDescriptor = getattr(fileObjectOrDescriptor, "fileno", lambda: fileObjectOrDescriptor)()
  if not isinstance(fileDescriptor, int):
    raise ValueError(f"Expected a file-like object with '.fileno()' method or a file descriptor, but got {fileObjectOrDescriptor}")
  return fileDescriptor


try:
  import ctypes
  from ctypes.util import find_library
except ImportError:
  libc = None
else:
  try:
    libc = ctypes.cdll.msvcrt  # Windows
  except OSError:
    libc = ctypes.cdll.LoadLibrary(find_library("c"))
def flush(stream: TextIO) -> None:
  """Flushes all libc buffers and the buffer of the given text stream"""
  # from https://stackoverflow.com/a/22434262
  try:
    libc.fflush(None)  # `fflush(NULL)` flushes all open output streams managed by libc
    stream.flush()
  except (AttributeError, ValueError, IOError):
    pass  # unsupported


@contextmanager
def redirect(
  destStream: TextIO | int | str = os.devnull,
  srcStream:  TextIO             = sys.stdout,
) -> Generator[TextIO | None, None, None]:
  """Redirects the output of the text stream `fromStream` to the text stream, file descriptor, or file name `toStream`"""
  # from https://stackoverflow.com/a/22434262
  srcStreamFd = fileDescriptor(srcStream)
  # save source-stream file object before it is overwritten
  #NOTE: `srcStreamCopy` is inheritable on Windows when duplicating a standard stream
  with os.fdopen(os.dup(srcStreamFd), "wb") as srcStreamCopy:
    flush(srcStream)  # flush library buffers that dup2() knows nothing about
    try:
      os.dup2(fileDescriptor(destStream), srcStreamFd)  # set srcStream to destStream
    except ValueError:  # destStream may be a file name
      with open(destStream, "wb") as destFile:
        os.dup2(destFile.fileno(), srcStreamFd)  # set srcStream to destStream
    try:
      yield srcStream  # run code with the redirected srcStream
    finally:
      flush(srcStream)
      # restore srcStream to its previous value
      #NOTE: dup2() makes stdout_fd inheritable unconditionally
      os.dup2(srcStreamCopy.fileno(), srcStreamFd)  # set srcStream to srcStreamCopy


@contextmanager
def mergeStderrIntoStdout() -> Generator[TextIO | None, None, None]:
  """Redirects stderr to stdout"""
  # from https://stackoverflow.com/a/22434262
  with redirect(destStream = sys.stdout, srcStream = sys.stderr):  # equivalent to $ exec 2>&1
    yield


if __name__ == "__main__":

  fitRootDir = "./fits"
  # fitRootDir = "./fits.pionComparison"
  # fitRootDir = "./fits.pionComparison.R6.28"
  dataPeriods = (
    # "2017_01-ver03",
    # "2018_01-ver02",
    # "2018_08-ver02",
    # "2019_11-ver01",
    "2017_01-ver03_goodToF",
    # "2018_01-ver02_goodToF",
    # "2018_08-ver02_goodToF",
  )
  useTotal = False

  dataSamples: list[dict[str, str]] = []
  for dataPeriod in dataPeriods:
    dataSamples += [
      # { # bggen MC
      #   **dict.fromkeys(["dataFileName", "bggenFileName"],
      #                   f"./data/MCbggen/{dataPeriod}/pippippimpimpmiss_flatTree.MCbggen_{dataPeriod}.root.brufit"),  # "dataFileName" and "bggenFileName" are set to identical values
      #   "dataPeriod" : dataPeriod,
      #   "dataLabel"  : f"bggen_{dataPeriod}",
      # },
      { # real data
        "dataFileName"  : f"./data/RD/{dataPeriod}/pippippimpimpmiss_flatTree.RD_{dataPeriod}.root.brufit",
        "bggenFileName" : f"./data/MCbggen/{dataPeriod}/pippippimpimpmiss_flatTree.MCbggen_{dataPeriod}.root.brufit",
        "dataPeriod"    : dataPeriod,
        "dataLabel"     : f"data_{dataPeriod}",
      },
    ]

  fits: list[list[dict[str, Any]]] = [[  # list of fits for each data sample
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFixed_noBkg",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "\"\"",  #TODO is it correct to define a string '""'? why not just empty string?
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFudge_noBkg",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "\"\"",
    # },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFixed",
      "kwargs" : {
        "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
        "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
      },
    },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "shift",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFudge",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgAllFudge",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFudge",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "Histogram",
    # },
  ] for dataSample in dataSamples]
  # add data period and input files (same for all fits of a given data sample)
  for index, dataSample in enumerate(dataSamples):
    for fit in fits[index]:
      fit.update({"dataPeriod" : dataSample["dataPeriod"]})
      fit["kwargs"]["dataFileName"]  = dataSample["dataFileName"]
      fit["kwargs"]["bggenFileName"] = dataSample["bggenFileName"]

  # pdfTypeBkg = "DoubleGaussian",
  # pdfTypeBkg = "DoubleGaussian_SameMean",
  # pdfTypeBkg = "SkewedGaussian_SkewNormal",
  # pdfTypeBkg = "SkewedGaussian_ExpMod",
  # pdfTypeBkg = "SkewedGaussian_Log",

  ROOT.gBenchmark.Start("Total processing time")
  for fitsForDataSample in fits:
    for fit in fitsForDataSample:
      fitDirectory = f"{fitRootDir}/{fit['dataPeriod']}/noShowers/{fit['fitDirectory']}"
      # recreate directories if already existing
      shutil.rmtree(fitDirectory, ignore_errors = True)
      os.makedirs(fitDirectory, exist_ok = True)
      print(f"Created directory '{fitDirectory}'")
      # run fits
      print(f"Starting fits ...")
      with mergeStderrIntoStdout():
        with redirect(destStream = f"{fitDirectory}/fitMissingMassSquared.log"):
          fitMissingMassSquared(
            outputDirName = fitDirectory,
            **fit["kwargs"],
          )
      # with open(f"{fitDirectory}/fitMissingMassSquared.log", "w") as logFile:  # write separate log file for each fit
      #   with contextlib.redirect_stdout(logFile), contextlib.redirect_stderr(logFile):
      #     fitMissingMassSquared(
      #       outputDirName = fitDirectory,
      #       **fit["kwargs"],
      #     )
      # cmdLineOptions = [(f"--{option} {fit[option]}" if option in fit else "") for option in ("pdfTypeSig", "fixParsSig", "pdfTypeBkg", "fixParsBkg")]
      # result = subprocess.run(
      #   f"./fitMissingMassSquared.py \"{fitDirectory}\" {fit['dataFileName']} {fit['bggenFileName']} {' '.join(cmdLineOptions)} &> \"{fitDirectory}/fitMissingMassSquared.log\"", shell = True)
      # if result.returncode != 0:
      #   raise RuntimeError(f"Fitting script failed with exit code '{result.returncode}'")
      # postprocess fit results
      subprocess.run(f"./cleanFitDir.sh \"{fitDirectory}\"", shell = True)
      print("Plotting fit results...")
      #TODO call python functions directly instead of making the detour via the command-line interface
      subprocess.run(f"./plotFitResults.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotFitResults.log\"", shell = True)
      print("Plotting efficiencies...")
      subprocess.run(f"./plotEfficiencies.py {'--useTotal' if useTotal else ''} \"{fitDirectory}\"  &> \"{fitDirectory}/plotEfficiencies.log\"", shell = True)
  ROOT.gBenchmark.Show("Total processing time")
