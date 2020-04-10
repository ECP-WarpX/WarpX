import os
import re
import scipy.constants as scc
import numpy as np


def _set_value(str_text, str_before, str_after, str_arg, val):
    str_line = str_before + str_arg + str_after
    str_text = re.sub('\n' + str_line + '.*',
                      '\n' + str_line + str(val),
                      str_text)
    return str_text


def write_sim_input(input_file, parameters):

    # Parameters exposed to optimization
    ramp_down_1 = parameters[0][0]
    ramp_down_2 = parameters[0][1]
    zlens_1 = parameters[0][2]
    adjust_factor = parameters[0][3]

    # Fixed parameters
    ramp_up_1 = 0.02
    plateau_1 = 0.297
    ramp_up_2 = ramp_up_1
    plateau_2 = plateau_1
    gap_12 = .0285

    # ramp_down_1 = 0.003
    # ramp_down_2 = 0.003
    # zlens_1 = 0.0175
    # adjust_factor = 1.

    end_stage_1 = ramp_up_1 + plateau_1 + ramp_down_1
    beg_stage_2 = end_stage_1 + gap_12
    end_stage_2 = beg_stage_2 + ramp_up_2 + plateau_2 + ramp_down_2

    # End simulation when beam has just escaped the last stage
    gamma_b = 30.
    zmax_stop_run = end_stage_2 - 55.e-6 * gamma_b**2 * 2.

    with open(input_file) as file_handler:
        output_text = file_handler.read()

    # Set end of stage 1
    output_text = _set_value(
        output_text, '', ' = ', 'electrons.zmax', str(end_stage_1))
    output_text = _set_value(
        output_text, '', ' = ', 'ions.zmax', str(end_stage_1))
    # Set length of final downramp of stage 1
    output_text = _set_value(
        output_text, '', ' = ', 'electrons.predefined_profile_params',
        '0.0 .02 .297 ' + str(ramp_down_1) + ' 40.e-6 1.7e23')
    output_text = _set_value(
        output_text, '', ' = ', 'ions.predefined_profile_params',
        '0.0 .02 .297 ' + str(ramp_down_1) + ' 40.e-6 1.7e23')
    # Set position of lens
    output_text = _set_value(
        output_text, '', ' = ', 'my_constants.zlen', str(end_stage_1 + zlens_1)
    )
    # Set beginning of stage 2
    output_text = _set_value(
        output_text, '', ' = ', 'electrons2.zmin', str(beg_stage_2))
    output_text = _set_value(
        output_text, '', ' = ', 'ions2.zmin', str(beg_stage_2))
    # Set end of stage 2
    output_text = _set_value(
        output_text, '', ' = ', 'electrons2.zmax', str(end_stage_2))
    output_text = _set_value(
        output_text, '', ' = ', 'ions2.zmax', str(end_stage_2))
    # Set length of final downramp of stage 2
    output_text = _set_value(
        output_text, '', ' = ', 'electrons2.predefined_profile_params',
        str(beg_stage_2) + ' .02 .297 ' + str(ramp_down_2) + ' 40.e-6 1.7e23')
    output_text = _set_value(
        output_text, '', ' = ', 'ions2.predefined_profile_params',
        str(beg_stage_2) + ' .02 .297 ' + str(ramp_down_2) + ' 40.e-6 1.7e23')
    # Set adjustment factor on lens strength
    output_text = _set_value(
        output_text, '', ' = ', 'my_constants.adjust_factor',
        str(adjust_factor))
    # Set when to stop the run
    output_text = _set_value(
        output_text, '', ' = ', 'warpx.zmax_plasma_to_compute_max_step',
        str(zmax_stop_run))

# Write new input file
    fout = open(input_file, 'w')
    fout.write(output_text)
    fout.close()
