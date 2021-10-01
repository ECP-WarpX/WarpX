"""File to hold utility functions associated with setting up and plotting time
varying voltage functions."""

import numpy as np
from functools import partial
from collections import namedtuple

from mewarpx.mwxrun import mwxrun
from pywarpx import callbacks


# namedtuple is a convenient way of making a class with given properties
PulseFunction = namedtuple(
    'PulseFunction',
    ['get_voltage', 'V_off', 'V_on', 'pulse_period', 'pulse_length',
     't_rise', 't_fall', 'wait_time']
)


def linear_pulse_function(V_off, V_on, pulse_period, pulse_length, t_rise=25e-9,
                          t_fall=25e-9, wait_time=5e-9, plot=True):
    r"""Function to construct an expression for a voltage pulse that can be
    parsed by the WarpX parser. This version uses linear segments to build a
    pulse with form
          ____                               ____
         /    \                             /    \
    ____/      \___________________________/      \_____

    Arguments:
        V_off (float): Electrode bias in the off-phase in Volt.
        V_on (float): Electrode voltage at the peak of the pulse in Volt.
        pulse_period (float): Full period from pulse start to the subsequent
            pulse start in seconds.
        pulse_length (float): Time duration of on-pulse phase in seconds. Does
            not include t_rise and t_fall.
        t_rise (float): Time for pulse to rise to its on value in
            seconds. Defaults to 25 ns.
        t_fall (float): Time for pulse to decrease to its off value in
            seconds. Defaults to 25 ns.
        wait_time (float): Time in seconds to offset the pulse function back
            to allow the simulation time to equilibrate before starting the
            pulse. Defaults to 5 ns.
        plot (bool): Optional parameter, if True the pulse waveform will be
            plotted.

    Returns:
        String expression of voltage function that can be parsed by the WarpX
        parser.
    """
    expr = (
        f"(modT=(t-{wait_time}); " # f"modT=(t-{wait_time})%{pulse_period}; " should insert modulo when supported
        f"p1=(modT>=0 and modT<{t_rise}); "
        f"p2=(modT>={t_rise} and modT<({t_rise}+{pulse_length})); "
        f"p3=(modT>=({t_rise}+{pulse_length}) and modT<({t_rise}+{pulse_length}+{t_fall})); "
        "p4=if(p1 or p2 or p3, 0, 1); "
        f"dV={V_on}-{V_off}; "
        f"p1*({V_off}+dV/{t_rise}*modT)+p2*{V_on}+"
        f"p3*({V_on}-dV/{t_fall}*(modT-{t_rise}-{pulse_length}))+p4*{V_off})"
    )

    if plot:
        def get_voltage(times):
            voltages = np.zeros_like(times)
            for ii, t in enumerate(times):
                voltages[ii] = mwxrun.eval_expression_t(expr, t)
            return voltages

        pulse_function = PulseFunction(
            get_voltage=get_voltage,
            V_off=V_off, V_on=V_on, pulse_period=pulse_period,
            pulse_length=pulse_length, t_rise=t_rise,
            t_fall=t_fall, wait_time=wait_time
        )
        f = partial(plot_pulse, pulse_function)
        f.__name__ = 'pulse_plot'
        callbacks.installafterinit(f)

    return expr


def plot_pulse(pulse_function, save_dir='diags/circuit'):
    """Utility function that plots the pulse shape applied to the device with
    annotations showing pulse parameters.

    Arguments:
        pulse_function (PulseFunction): Contains pulse details as well as a
            callable function `get_voltage` to get the pulse voltage at a given
            time.
        save_dir (str): Directory where the pulse schematic should be saved.
    """
    import os
    import matplotlib as ml
    import matplotlib.pyplot as plt

    from mewarpx.utils_store import util

    util.mkdir_p(save_dir)

    fig = plt.figure(figsize=(8, 7))
    font = {'size': 13}
    ml.rc('font', **font)
    fig.suptitle(
        'Schematic of pulse shape with V$_{on}$ = %.1f V, '
        'V$_{off}$ = %.1f V\nPulse period = %s, '
        'Pulse length = %.0f ns' % (
            pulse_function.V_on,
            pulse_function.V_off,
            (r"$\infty$" if
             np.isinf(pulse_function.pulse_period) else r"%.0f $\mu$s" %
             (pulse_function.pulse_period*1e6)),
            pulse_function.pulse_length*1e9
        )
    )

    plt.subplot(2, 1, 1)

    if hasattr(pulse_function, 'wait_time'):
        wait_time = pulse_function.wait_time
    else:
        wait_time = 0.0

    arrow_style = dict(arrowstyle="<->", connectionstyle="arc3")
    text_style = {
        'color': 'black', 'fontsize': 11, 'ha': 'center', 'va': 'center'
    }
    text_box_style = dict(boxstyle="round", fc="white", ec="black", pad=0.2)

    # Create top plot showing multiple pulse periods
    t_grid = np.linspace(
        max(-pulse_function.pulse_period*1.8, -10e-6),
        min(pulse_function.pulse_period*2.2, 50e-6), 50000
    ) + wait_time
    plt.plot(t_grid*1e6, pulse_function.get_voltage(t_grid))

    # Add arrow showing pulse period
    plt.annotate(
        "",
        xy=(wait_time*1e6, pulse_function.V_on*0.8),
        xytext=(
            min((pulse_function.pulse_period + wait_time), t_grid[-1])*1e6,
            pulse_function.V_on*0.8
        ),
        arrowprops=arrow_style
    )
    # Add text annotation to pulse period arrow
    plt.text(
        (wait_time + 0.5*min(pulse_function.pulse_period, t_grid[-1]))*1e6,
        pulse_function.V_on*0.8,
        'pulse period',
        text_style,
        bbox=text_box_style
    )

    # Add arrow showing pulse voltage
    plt.annotate(
        "",
        xy=((t_grid[0] - t_grid[0]*0.5)*1e6, pulse_function.V_off),
        xytext=((t_grid[0] - t_grid[0]/2.0)*1e6, pulse_function.V_on),
        arrowprops=arrow_style
    )
    # Add text annotation to pulse voltage arrow
    plt.text(
        (t_grid[0] - t_grid[0]/2.0)*1e6,
        (pulse_function.V_on + pulse_function.V_off)*0.5,
        'V$_{pulse}$',
        text_style,
        bbox=text_box_style
    )

    plt.ylabel('Applied voltage (V)')
    plt.xlabel(r'Time ($\mu s$)')
    plt.grid()

    plt.subplot(2, 1, 2)

    # Create bottom plot showing pulse details
    t_grid = np.linspace(
        -0.01e-6,
        pulse_function.t_rise + pulse_function.pulse_length +
        pulse_function.t_fall + 5e-9, 500
    )
    plt.plot(t_grid*1e9, pulse_function.get_voltage(t_grid + wait_time))

    # Add arrow showing rise time
    plt.annotate(
        "",
        xy=(0.0, pulse_function.V_on*0.5),
        xytext=(
            pulse_function.t_rise*1e9, pulse_function.V_on*0.5
        ),
        arrowprops=arrow_style
    )
    # Add text annotation to rise time arrow
    plt.text(
        (pulse_function.t_rise/2.)*1e9, pulse_function.V_on*0.5,
        't$_{rise}$',
        text_style,
        bbox=text_box_style
    )

    # Add arrow showing fall time
    plt.annotate(
        "",
        xy=(
            1e9*(pulse_function.t_rise+pulse_function.pulse_length),
            pulse_function.V_on * 0.5
        ),
        xytext=(
            (pulse_function.t_rise + pulse_function.pulse_length
             + pulse_function.t_fall)*1e9, pulse_function.V_on * 0.5
        ),
        arrowprops=arrow_style
    )
    # Add text annotation to fall time arrow
    plt.text(
        (pulse_function.t_rise + pulse_function.pulse_length
         + pulse_function.t_fall/2.)*1e9,
        pulse_function.V_on*0.5,
        't$_{fall}$',
        text_style,
        bbox=text_box_style
    )

    # Add arrow showing pulse length
    plt.annotate(
        "",
        xy=(
            pulse_function.t_rise*1e9, pulse_function.V_on*0.8
        ),
        xytext=(
            (pulse_function.t_rise+pulse_function.pulse_length)*1e9,
            pulse_function.V_on*0.8
        ),
        arrowprops=arrow_style
    )
    # Add text annotation to pulse length arrow
    plt.text(
        (pulse_function.t_rise
         + pulse_function.pulse_length/2.)*1e9,
        pulse_function.V_on*0.8,
        'pulse length',
        text_style,
        bbox=text_box_style
    )

    plt.ylabel('Applied voltage (V)')
    plt.xlabel('t - %.1e s (ns)' % wait_time)
    plt.grid()
    plt.subplots_adjust(top=0.91, hspace=0.22)
    plt.savefig(os.path.join(save_dir, 'pulse_schematic.png'), dpi=600)
    plt.close(fig)
