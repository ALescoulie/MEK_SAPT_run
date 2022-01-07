from os import path
import sys

from MDSAPT.mdsapt import Reader, Optimizer, TrajectorySAPT
from MDSAPT.mdsapt.reader import InputReader

if __name__ == '__main__':
    Settings = InputReader(sys.argv[1])

    Opt = Optimizer(Settings)

    SAPTRun = TrajectorySAPT(Settings, Optimizer).run(Settings.start, Settings.stop, Settings.step)

    SAPTRun.results.to_csv(sys.argv[2])
