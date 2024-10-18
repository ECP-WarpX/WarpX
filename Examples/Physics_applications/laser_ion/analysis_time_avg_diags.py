#!/usr/bin/env python3

import sys
import numpy as np
import openpmd_api as io


def load_field_from_iteration(series, iteration : int, field : str, coord : str = None) -> np.ndarray:
    """Load iteration of field data from file."""

    it = series.iterations[iteration]
    field_obj = it.meshes[f"{field}"]

    if field_obj.scalar:
        field_data = field_obj[io.Mesh_Record_Component.SCALAR].load_chunk()
    elif coord in [item[0] for item in list(field_obj.items())]:
        field_data = field_obj[coord].load_chunk()
    else:
        raise Exception(f"Specified coordinate: f{coord} is not available for field: f{field}.")
    series.flush()

    return field_data


def main():

    field = "E"
    coord = "z"
    avg_period_steps = 5
    avg_output_step = 100

    path_tpl_inst = "diags/diagInst/openpmd_%T.h5"

    dir_avg = sys.argv[1]
    path_tpl_avg = f"{dir_avg}/openpmd_%T.h5"

    si = io.Series(path_tpl_inst, io.Access.read_only)
    sa = io.Series(path_tpl_avg, io.Access.read_only)

    ii0 = si.iterations[0]
    fi0 = ii0.meshes[field][coord]
    shape = fi0.shape

    data_inst = np.zeros(shape)

    for i in np.arange(avg_output_step-avg_period_steps+1,avg_output_step+1):
        data_inst += load_field_from_iteration(si, i, field, coord)

    data_inst = data_inst / avg_period_steps

    data_avg = load_field_from_iteration(sa, avg_output_step, field, coord)

    # Compare the data
    if np.allclose(data_inst, data_avg):
        print("Test passed: actual data is close to expected data.")
    else:
        print("Test failed: actual data is not close to expected data.")
        sys.exit(1)


if __name__ == "__main__":
    main()