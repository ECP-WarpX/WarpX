"""Test util.py functions."""
import os
import shutil

import numpy as np
import pytest

import mewarpx
from mewarpx.utils_store import (oracle_control, plasma_density_oracle,
                                 testing_util)
from mewarpx.utils_store import util as mwxutil


def test_utils_check_version(caplog):
    """Note the version passed in here is the "script" version."""
    with pytest.raises(
        ValueError, match="This version of mewarpx .* is older"
    ):
        mwxutil.check_version((np.inf, 0, 0), -1)

    with pytest.raises(
        ValueError, match="This physics version of mewarpx .* is older"
    ):
        mwxutil.check_version((0, 0, 9), np.inf)

    # Check warnings printed
    caplog.clear()
    mwxutil.check_version(
        (-1, 0, 0), mewarpx.__physics_version__)
    assert len(caplog.record_tuples) == 1
    assert "is a newer API version" in caplog.text

    caplog.clear()
    mwxutil.check_version(
        mewarpx.__version_info__, -1
    )
    assert len(caplog.record_tuples) == 1
    assert "is newer than the version" in caplog.text

    # Check an OK setup
    caplog.clear()
    mwxutil.check_version(
        (mewarpx.__version_info__[0],
         mewarpx.__version_info__[1] - 1,
         mewarpx.__version_info__[2] - 1,
         ),
        mewarpx.__physics_version__
    )
    assert len(caplog.record_tuples) == 0


def test_plasma_density_oracle():
    testing_util.initialize_testingdir("test_plasma_density_oracle")
    ref_dir = os.path.join(testing_util.test_dir, "utils", "plasma_density_oracle")

    plasma_oracle = plasma_density_oracle.PlasmaDensityOracle(
        start_step=100000,
        interval=100000,
        predict_step=700000,
        species_names=["ion"],
        coarsening_factor=2,
        work_dir=ref_dir,
        save_dir=os.curdir
    )

    plasma_oracle.get_coarse_data()
    plasma_oracle.generate_regression()

    ref_prediction = np.load(
        os.path.join(ref_dir, "ion_particle_density_prediction_700000.npy")
    )

    assert np.allclose(plasma_oracle.species_list[0]["predict"],
                       ref_prediction, rtol=1e-6)


def test_oracle_control():
    testing_util.initialize_testingdir("test_oracle_control")
    ref_dir = os.path.join(testing_util.test_dir, "utils", "plasma_density_oracle")
    shutil.copytree(ref_dir, "./start_0")
    runfile_path = "run_simulation.py"

    # generate fake runfile
    with open(os.path.join(runfile_path), "w") as runfile:
        runfile.write("TOTAL_TIME = WRONG_VALUE")

    # generate results.txt
    with open(os.path.join("start_0", "diags", "results.txt"), "w"):
        pass

    controller = oracle_control.OracleControl(
        diag_steps=50000,
        dt=5e-13,
        first_data=50e-9,
        override_bins=[[0], [50e-9, 50e-9], [350e-9, 350e-9]],
        run_script=runfile_path,
        species_names=["ar_ions"]
    )

    controller.main()

    with open(runfile_path, "r") as runfile:
        total_time = float(runfile.readline().rstrip().split(" = ")[1])
        assert np.allclose(total_time, 150e-9 + 5e-13, rtol=1e-8, atol=1e-15)

    with open("oracle_control.txt", "r") as control_file:
        run_dir_path = control_file.readline().rstrip()
        expected_path = os.path.abspath(os.path.join(os.curdir, "start_1000000"))
        assert run_dir_path == expected_path

    assert os.path.isfile(os.path.join(
        os.curdir, "seed_density", "ar_ions_particle_density_prediction.npy"
    ))
