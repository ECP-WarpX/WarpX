"""
This script reads in plasma densities at 3 time points and uses exponential
regression, linear regression, or the mean, to predict the plasma density at a
future time point. This script is used to speed up the process of finding the
steady state of a simulation.
"""
import argparse
import glob
import logging
import os

import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.measure import block_reduce

import mewarpx.utils_store.util as mwxutil

logger = logging.getLogger(__name__)


class PlasmaDensityOracle:

    def __init__(self, start_step, interval, predict_step, species_names,
                 coarsening_factor=None, work_dir=None, save_dir=None):
        """
        Arguments:
            start_step (int): time step of first data point for exponential
                regression
            interval (int): number of time steps between data points
            predict_step (int): time step of prediction
            species_names (list): list of species names or identifying shared
                substrings of species names, eg. electrons or ion. The glob
                syntax using species name would be
                `{species_name}*_particle_density_{self.time_steps[i]:010d}.npy`
            coarsening_factor (int): the down-sampling integer factor along
                each axis for the grid, defaults to 8
            work_dir (str): path to parent directory of the diags directory,
                defaults to current directory
            save_dir (str): path to directory to save output files to, defaults
                to work_dir/../seed_density
        """
        self.time_steps = [
            start_step, start_step + interval, start_step + 2 * interval
        ]
        self.time_skip = predict_step

        if work_dir is None:
            work_dir = os.curdir
        if save_dir is None:
            save_dir = os.path.join(work_dir, "..", "seed_density")
        self.diags_dir = os.path.join(work_dir, "diags", "fields")
        self.output_dir = save_dir

        if coarsening_factor is None:
            coarsening_factor = 8
        self.coarsening_factor = coarsening_factor

        self.species_list = []
        for name in species_names:
            self.species_list.append(
                {"name": name}
            )

        logger.info(f"Generating prediction using data from steps "
                    f"{', '.join([str(x) for x in self.time_steps])}")
        logger.info(f"Prediction step is {self.time_skip}")

    def main(self):
        """Calls all functions in order to get raw data, generate prediction,
        and save prediction"""
        self.get_coarse_data()
        self.generate_regression()
        self.save_predictions()

    def get_coarse_data(self):
        """Gets data from diag files and coarsens them"""
        for species in self.species_list:
            for i in range(len(self.time_steps)):
                # pattern match files to combine densities from multiple files
                files = glob.glob(os.path.join(
                    self.diags_dir,
                    f"*{species['name']}*_particle_density_{self.time_steps[i]:010d}.npy"
                ))
                temp_data = np.load(files[0])
                for j in range(1, len(files)):
                    temp_data = temp_data + np.load(files[j])

                # slice off first value to remove ghost cell
                species[f"data_{i}"] = block_reduce(
                    temp_data[1:, 1:], self.coarsening_factor, np.mean
                )

    def generate_regression(self):
        """Analyzes regressions for each species and calculates the predicted
        density for each species at the time_skip timestep."""
        for species in self.species_list:
            num_rows = species["data_0"].shape[0]
            num_columns = species["data_0"].shape[1]
            density_max = np.amax(
                [species["data_0"], species["data_1"], species["data_2"]]
            )

            predict_array = np.zeros_like(species["data_0"])

            for i in range(num_rows):
                for j in range(num_columns):
                    # scale data to avoid overflow error
                    x_data = np.array(self.time_steps) / self.time_skip
                    y_data = np.array([species["data_0"][i, j],
                                       species["data_1"][i, j],
                                       species["data_2"][i, j]]) / density_max

                    # check if data is monotonically increasing or decreasing
                    if (np.all(y_data[1:] > y_data[:-1]) or
                       np.all(y_data[1:] < y_data[:-1])):
                        # use linear regression if plasma density is accelerating
                        # or if plasma density is linear
                        if (np.abs(y_data[2] - y_data[1])
                            >= np.abs(y_data[1] - y_data[0])):
                            m, b = self._solve_linear(x_data, y_data)
                            # time_skip is at x = 1, due to scaling down
                            prediction = m + b
                        # if density is decelerating
                        else:
                            constants = self._solve_inverse_exponential(x_data,
                                                                        y_data)
                            prediction = self._inverse_exponential(1,
                                                                   *constants)
                    # if data is not oscillating, use linear regression
                    elif (np.all(y_data[1:] >= y_data[:-1]) or
                          np.all(y_data[1:] <= y_data[:-1])):
                        m, b = self._solve_linear(x_data, y_data)
                        prediction = m + b
                    # use mean if density is oscillating
                    else:
                        prediction = np.mean(y_data)

                    predict_array[i, j] = prediction * density_max

            species["predict"] = gaussian_filter(predict_array, 1)

    def save_predictions(self):
        """Save the generated predictions in the output directory"""
        if not os.path.isdir(self.output_dir):
            mwxutil.mkdir_p(self.output_dir)

        for species in self.species_list:
            np.save(
                os.path.join(
                    self.output_dir,
                    f"{species['name']}_particle_density_prediction.npy"
                ),
                species["predict"]
            )

    @staticmethod
    def _solve_linear(x_data, y_data):
        """Solves for the constants in a linear regression, given 3 x and y
        values. Only the first and last x and y values will be used.

        Arguments:
            x_data (array like): three x values
            y_data (array like): three y values corresponding to the x values
        Returns:
            constants (list): list of constants that were solved for in the
                linear regression equation
        """
        m = (y_data[-1] - y_data[0]) / (x_data[-1] - x_data[0])
        b = y_data[0] - m * x_data[0]
        return [m, b]

    @staticmethod
    def _solve_inverse_exponential(x_data, y_data):
        """Solves for the constants in the inverse exponential equation using
        3 given x and y values. x_data and y_data must be monotonically
        increasing.

        Arguments:
            x_data (array like): three x values
            y_data (array like): three y values corresponding to the x values

        Returns:
            constants (list): list of constants that were solved for, in the
                inverse exponential equation
        """
        delta = (x_data[1] - x_data[0])

        b = -np.log(
            (y_data[2] - y_data[0]) / (y_data[1] - y_data[0]) - 1
            ) / delta
        a = ((y_data[1] - y_data[0]) / (1 - np.exp(-1 * delta * b))
             * np.exp(b * x_data[0]))
        c = y_data[0] - a * (1 - np.exp(-b * x_data[0]))

        return [a, b, c]

    @staticmethod
    def _inverse_exponential(x, a, b, c):
        """Negative inverse exponential equation used to for the regression."""
        return a * (1 - np.exp(-b * x)) + c


def entry():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--start-step", type=int, metavar="timestep", required=True,
        dest="start_step", help="The timestep of the first input data point."
    )
    parser.add_argument(
        "--interval", type=int, metavar="num_steps", required=True,
        help="The number of steps between input data points."
    )
    parser.add_argument(
        "--predict-step", type=int, metavar="timestep", required=True,
        dest="predict_step",
        help="The number of timestep intervals to predict the particle density "
             "at."
    )
    parser.add_argument(
        "--species", type=str, metavar="species_name", required=True,
        nargs="+", dest="species_names",
        help="Names of species or identifying shared substrings of species, "
             "eg. electrons or Xe_ion"
    )
    parser.add_argument(
        "--coarse-factor", type=int, metavar="n", dest="coarsening_factor",
        help="The down-sampling integer factor along each axis for the grid. "
             "Defaults to 8"
    )
    parser.add_argument(
        "--run-dir", type=str, metavar="path", dest="work_dir",
        help="Path to the parent directory of the diags directory. Defaults to "
             "the current directory."
    )
    parser.add_argument(
        "--save-dir", type=str, metavar="path", dest="save_dir",
        help="Path to where to save output files to. Defaults to "
             "work_dir/../seed_density/"
    )

    args = vars(parser.parse_args())

    predictor = PlasmaDensityOracle(**args)
    predictor.main()


if __name__ == "__main__":
    entry()
