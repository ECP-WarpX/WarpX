"""
This script is a wrapper for plasma_density_oracle that helps automate use
of plasma_density_oracle. This script calculates the best parameters to use for
plasma_density_oracle and then uses those parameters to generate a prediction.
"""
import argparse
import bisect
import fileinput
import glob
import logging
import os

import numpy as np

from mewarpx.utils_store.plasma_density_oracle import PlasmaDensityOracle
import mewarpx.utils_store.util as mwxutil

logger = logging.getLogger(__name__)


class OracleControl:

    def __init__(self, diag_steps, dt, first_data, override_bins=None,
                 work_dir=None, run_script=None, save_dir=None, prefix=None,
                 species_names=None, coarsening_factor=None):
        """
        Arguments:
            diag_steps (int): number of time steps in a diagnostic interval
            dt (float): time of timestep in seconds
            first_data (float): time in seconds of the first data point used in
                the regression
            override_bins (list): list of 3 lists that contains override
                variables for determining prediction step parameters:
                - list 0: list of start time bins
                - list 1: list of possible interval values
                - list 2: list of prediction time values
                Defaults to
                    [
                        [1e-6, 3e-6, 15e-6],
                        [100e-9, 200e-9, 500e-9, 500e-9],
                        [1e-6, 2e-6, 4e-6, 8e-6]
                    ]
            work_dir (str): directory path where the runs are located, defaults
                to current directory
            save_dir (str): directory path to save generated prediction files,
                defaults to "work_dir/seed_density"
            prefix (str): prefix of run directories, defaults to "start"
            species_names (list): list of species names or identifying shared
                substrings of species names, eg. electrons or Xe_ion. The glob
                syntax using species name would be
                `{species_name}*_particle_density_{self.time_steps[i]:010d}.npy`.
                Defaults to ["electron", "ion"]
            coarsening_factor (int): the down-sampling integer factor along
                each axis for the grid. Defaults to 8
        """
        if work_dir is None:
            self.work_dir = os.path.abspath(os.curdir)
        else:
            self.work_dir = work_dir

        self.output_file = os.path.join(self.work_dir, "oracle_control.txt")

        # remove output file if it exists so that the autorun script will exit
        # if any unexpected errors are thrown in this script
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

        # set time values
        self.diag_interval = diag_steps
        self.dt = dt
        self.start_step = self._time_to_step(first_data)

        # set class variables
        if run_script is None:
            run_script = "run_simulation.py"
        if prefix is None:
            prefix = "start"
        if species_names is None:
            species_names = ["electron", "ion"]
        if coarsening_factor is None:
            coarsening_factor = 8
        self.run_script = os.path.join(self.work_dir, run_script)
        self.prefix = prefix
        self.species_names = species_names
        self.coarsening_factor = coarsening_factor
        self.save_dir = save_dir
        self.interval = None
        self.predict_step = None

        # value to be logged later
        self.absolute_predict_time = None
        self.new_run_step = None

        # value will be set if prediction is generated
        self.new_run = None

        # get directory of latest run
        dir_list = sorted(
            glob.glob(os.path.join(self.work_dir, f"{self.prefix}_*")),
            key=lambda fname: int(
                os.path.basename(fname).strip(f"{self.prefix}_")
                )
        )
        if len(dir_list) > 0:
            self.latest_run = os.path.abspath(dir_list[-1])
        # create new directory if no runs are found
        else:
            self.latest_run = os.path.join(self.work_dir, f"{self.prefix}_0")

        # check results.txt to see if latest run completed successfully
        self.latest_run_done = os.path.isfile(
            os.path.join(self.latest_run, "diags", "results.txt")
        )

        # set override bins or default bins
        if override_bins is None:
            override_bins = [
                [1e-6, 3e-6, 15e-6],
                [100e-9, 200e-9, 500e-9, 500e-9],
                [1e-6, 2e-6, 4e-6, 8e-6]
            ]

        self.time_bins = override_bins[0]
        self.interval_bins = override_bins[1]
        self.predict_bins = override_bins[2]

        if not (len(self.interval_bins) == len(self.time_bins) + 1 and
                len(self.predict_bins) == len(self.time_bins) + 1):
            raise ValueError("jump_vals and interval_vals must both be 1 "
                             "element longer than time_bins")
        if self.latest_run_done:
            self.get_prediction_parameters()
        else:
            self.new_run_step = int(
                os.path.basename(self.latest_run).strip(f"{self.prefix}_")
            )

    def get_prediction_parameters(self):
        """Uses the timing bins to determine prediction parameters for
        the latest run."""
        full_start_time = int(
            os.path.basename(self.latest_run).strip(f"{self.prefix}_")
        ) * self.dt

        bin_idx = bisect.bisect(self.time_bins, full_start_time)
        self.interval = self._time_to_step(self.interval_bins[bin_idx])
        jump_time = self.predict_bins[bin_idx]
        predict_time = (mwxutil.mwx_round(full_start_time + jump_time, 0.5e-6)
                        - full_start_time)
        self.predict_step = self._time_to_step(predict_time)

        self.absolute_predict_time = predict_time + full_start_time
        self.new_run_step = self._time_to_step(self.absolute_predict_time)
        run_dir = os.path.abspath(os.path.dirname(self.latest_run))

        # output new run directory
        self.new_run = os.path.join(
            run_dir, f"{self.prefix}_{str(self.new_run_step)}"
        )

    def predict_run(self):
        """Uses plasma density oracle to get plasma density prediction
        for the latest run. Should only be called if latest run completed
        successfully."""
        logger.info(f"Prediction control: Generating predictions for "
                    f"{self.absolute_predict_time * 1e6:.2f} us, using "
                    f"data from the run {self.prefix}_{self.latest_run}.")
        predictor = PlasmaDensityOracle(
            start_step=self.start_step, interval=self.interval,
            predict_step=self.predict_step, species_names=self.species_names,
            coarsening_factor=self.coarsening_factor, work_dir=self.latest_run,
            save_dir=self.save_dir
        )
        predictor.main()

        logger.info(f"Prediction control: Predictions generated for run "
                    f"{self.prefix}_{str(self.new_run_step)}.")

    def update_run_script(self):
        """If prediction was successful, update the run_simulation.py script's
        max time to allow for next prediction"""
        # find the bin corresponding to starting step of run
        bin_idx = bisect.bisect(self.time_bins, self.new_run_step * self.dt)
        new_run_total_time = (self.start_step * self.dt
                              + self.interval_bins[bin_idx] * 2)

        with fileinput.FileInput(self.run_script, inplace=True) as file:
            for line in file:
                spaces = len(line) - len(line.lstrip())
                line = line.rstrip()
                if "TOTAL_TIME = " in line:
                    contents = line.split(" = ")
                    contents[0] = f"{' ' * spaces}TOTAL_TIME"
                    contents[1] = str(new_run_total_time + 1 * self.dt)
                    line = " = ".join(contents)
                print(line)

    def _time_to_step(self, time):
        """Calculates the closest diag step corresponding to the given time.

        Arguments:
            time (float): simulation time to be converted into simulation diag
                step

        Returns:
            step (int): simulation diag step equivalent of time input
        """
        step = (int(np.round(time / (self.dt * self.diag_interval)))
                * self.diag_interval)
        return step

    def output_rundir(self):
        """Save run directory in output file so bash script knows where to run
        simulation."""
        if self.new_run is not None:
            output_dir = self.new_run
        else:
            output_dir = self.latest_run
        mwxutil.mkdir_p(os.path.dirname(self.output_file))
        with open(self.output_file, "w") as out_file:
            out_file.write(output_dir)

    def main(self):
        """Executes the complete logical steps of oracle control in order."""
        if self.latest_run_done:
            self.predict_run()
        self.update_run_script()
        self.output_rundir()


def entry():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--diag-steps", type=int, metavar="timesteps", required=True,
        dest="diag_steps",
        help="The number of timesteps in a diagnostic interval."
    )
    parser.add_argument(
        "--dt", type=float, metavar="sec", required=True, dest="dt",
        help="The time in seconds per timestep."
    )
    parser.add_argument(
        "--first-data", type=float, metavar="sec", required=True,
        dest="first_data",
        help="The time in seconds of the first data point used in the "
             "regression."
    )
    parser.add_argument(
        "-b", "--override-bins", type=float, metavar="list", nargs="+",
        dest="override_bins", action="append",
        help="Override bins for calculating prediction parameters. Must call 3 "
             "times to pass 'start time bins', 'interval values', "
             "'prediction time values', in that order. Defaults to "
             "["
                "[1e-6, 3e-6, 15e-6], "
                "[100e-9, 200e-9, 500e-9, 500e-9], "
                "[1e-6, 2e-6, 4e-6, 8e-6]"
             "]"
    )
    parser.add_argument(
        "--work-dir", type=str, metavar="path", dest="work_dir",
        help="Path to directory containing all run directories. Defaults to "
             "the current directory."
    )
    parser.add_argument(
        "--run-script", type=str, metavar="file", dest="run_script",
        help="Path of run script relative to the working directory. Defaults "
             "to 'run_simulation.py'"
    )
    parser.add_argument(
        "--save-dir", type=str, metavar="path", dest="save_dir",
        help="Path to where to save output files to. Defaults to "
        "work_dir/seed_density/"
    )
    parser.add_argument(
        "--prefix", type=str, metavar="prefix", dest="prefix",
        help="Prefix of run directories. Defaults to 'start'"
    )
    parser.add_argument(
        "--species", type=str, metavar="species_name", nargs="+",
        dest="species_names",
        help="Names of species or identifying shared substrings of species, "
             "eg. electrons or Xe_ion. Defaults to ['electron', 'ion']"
    )
    parser.add_argument(
        "--coarse-factor", type=int, metavar="n", dest="coarsening_factor",
        help="The down-sampling integer factor along each axis for the grid. "
             "Defaults to 8"
    )

    args = vars(parser.parse_args())

    controller = OracleControl(**args)
    controller.main()


if __name__ == "__main__":
    a = {'coarsening_factor': None,
         'diag_steps': 116129,
         'dt': 4.306e-13,
         'first_data': 1e-07,
         'override_bins': None,
         'prefix': None,
         'run_script': None,
         'save_dir': None,
         'species_names': None,
         'work_dir': None}
    controller = OracleControl(**a)
    controller.main()
