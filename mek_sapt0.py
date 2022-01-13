from os import path
import sys

from mdsapt import InputReader, Optimizer, TrajectorySAPT

if __name__ == '__main__':
    Settings = InputReader('atp_glu_leu_met.yaml')

    Opt = Optimizer(Settings)

    SAPTRun = TrajectorySAPT(Settings, Optimizer).run(Settings.start, Settings.stop, Settings.step)

    SAPTRun.results.to_csv('results.csv')
