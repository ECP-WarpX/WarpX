"""Base diagnostic code used in many diagnostics."""
import logging
import time

import numpy as np
from pywarpx import callbacks

from mewarpx.mwxrun import mwxrun
import mewarpx.utils_store.util as mwxutil

# import psutil



# Get module-level logger
logger = logging.getLogger(__name__)


class WarpXDiagnostic(object):

    """Hold attributes common to all diagnostics, especially logic on when to
    run.
    """

    DIAG_DIR = "diags"

    def __init__(self, diag_steps, diag_step_offset=0,
                 extended_interval_level=None,
                 manual_timesteps=None,
                 ):
        """Initialize timing. Will often be called within a child object's
        init.

        Arguments:
            diag_steps (int): Run the diagnostic with this period.
            diag_step_offset (int): Run when ``simulation_step % diag_steps =
                diag_step_offset``. Default 0.
            extended_interval_level (int): Enable an exponentially-increasing
                interval as the run progresses. If None (default), calculate
                every time period. Otherwise, the interval increases on an
                exponential scale. 0 is most aggressive at skipping diagnostic
                timesteps; 1, 2, etc. are progressively less aggressive but
                still exponential. Specifically, 0 executes once during each
                power-of-2 set of diagnostic steps; 1 twice; 2 four times, etc.
            manual_timesteps (list of ints): A list of timesteps to execute on
                within each diag_steps. If this list contains steps for the
                entire run, set diag_steps to a very large number. Default
                None.
        """
        self.diag_steps = int(round(diag_steps))
        self.diag_step_offset = int(round(diag_step_offset))

        self.extended_interval_level = extended_interval_level
        if self.extended_interval_level is not None:
            self.extended_interval_level = int(round(
                self.extended_interval_level))

        self.manual_timesteps = manual_timesteps
        if self.manual_timesteps is not None:
            self.manual_timesteps = [int(round(x))
                                     for x in self.manual_timesteps]

        if (
            (self.manual_timesteps is not None)
            and (self.extended_interval_level is not None)
        ):
            raise ValueError("manual_timesteps and extended_interval_level "
                             "cannot be used together; one must be None.")

        # handle creating the folder to save diagnostics
        if mwxrun.initialized:
            self._create_diag_folder()
        else:
            callbacks.installafterinit(self._create_diag_folder)

    def _create_diag_folder(self):
        """Helper function to create folder in which to save diagnostics. Child
        classes that save diagnostic files should have a ``write_dir``
        attribute specifying the directory in which it will save diagnostic
        files."""
        if hasattr(self, "write_dir") and mwxrun.me == 0:
            mwxutil.mkdir_p(self.write_dir)

    def check_timestep(self):
        """Check if the diagnostic should run on this timestep.

        Returns:
            l_execute (bool): If True, run on this timestep. If False, don't.
        """
        it = mwxrun.get_it()
        if self.manual_timesteps is not None:
            return (it % self.diag_steps) in self.manual_timesteps

        if (it % self.diag_steps) == self.diag_step_offset:
            return (
                (self.extended_interval_level is None)
                or self._comp_calcinterval(it)
            )

        return False

    def _comp_calcinterval(self, it):
        """Decide whether to calculate current based on exponential scaling
        algorithm.

        The algorithm forces the number of diagnostic period to be divisible by
        some power of 2 below it. If calcinterval_level = 0, it must be a power
        of 2 (e.g. for diagnum between 64 and 127, it must be divisible by 64
        => 64 only). For each increase in calcinterval_level, relax the
        division by one power of 2. For example, if calcinterval_level = 3 and
        diagnum is between 128 and 255, diagnum must be divisible by 128/2^3 =
        16.

        Also sets self.stride_multiplier.

        Arguments:
            it (int): The current integer step number
        """
        diagnum = it // self.diag_steps
        diagnum_log2 = np.floor(np.log2(diagnum))
        diagnum_divisor = int(max(
            2**(diagnum_log2 - self.extended_interval_level), 1))
        self.stride_multiplier = diagnum_divisor

        if diagnum % diagnum_divisor == 0:
            return True

        else:
            next_diagnum = ((diagnum // diagnum_divisor) + 1)*diagnum_divisor
            next_itnum = next_diagnum * self.diag_steps + self.diag_step_offset

            logger.info(
                f"Waiting until diagnostic period {next_diagnum} (step {next_itnum}) to run "
                f"diagnostic due to extended_interval_level = {self.extended_interval_level}"
            )

            return False


class TextDiag(WarpXDiagnostic):

    """Output diagnostics every certain number of steps.

    Contains:
        text_diag (function): Function to do the write-out, already inserted in
        warp.installafterstep(). Use only if needed for a future reference in
        script.
    """

    def __init__(self, diag_steps, preset_string='default',
                 custom_string=None, install=True, **kwargs):
        """Generate and install function to write out step #.

        Arguments:
            diag_steps (int): Number of steps between each output
            simulation (mespecies.Simulation): Main simulation object
            preset_string (str): Defaults to choose between:

              - ``default`` - just the step number and total particle num
              - ``perfdebug`` - like particledebug, plus interval wall time,
                step rate, and particle-step rate
              - ``memdebug`` - print out verbose memory usage information

            custom_string (str): Overrides preset_string if not None. The full
                string to output, with:

              - ``{step}`` formatting symbol for where the step number should
                go
              - ``{wall_time}`` run time of the last diag_steps steps
              - ``{step_rate}`` diag_steps / wall_time
              - ``{particle_step_rate}`` nplive * diag_steps / wall_time
              - ``{nplive}`` for number of live particles (global).
              - ``{npperspecies}`` for number of particles per species
                (global).
              - ``{iproc}`` for the current processor number
              - ``{system_memory}`` for verbose information on system memory
                usage.
              - ``{memory_usage}`` for memory usage of the current process
                only.

            install (bool): If False, don't actually install this into WarpX.
                Use if you want to call manually for debugging.
            kwargs: See :class:`mewarpx.mewarpx.diags_store.diag_base.WarpXDiagnostic`
                for more timing options.
        """
        # In warp we used a specific walltime counter it had
        # (warp.top.steptime). Not sure what issues we'll hit just using time
        # here.
        self.prev_time = time.time()
        self.start_time = self.prev_time
        self.prev_step = mwxrun.get_it()
        self.defaults_dict = {
            'default': "Step #{step:6d}; {nplive:8d} particles",
            'perfdebug': ("Step #{step:6d}; {nplive:8d} particles "
                          "{npperspecies} "
                          "{wall_time:6.1f} s wall time; "
                          "{step_rate:4.2f} steps/s; "
                          "{particle_step_rate:4.2f} particle*steps/s in the last {diag_steps} steps; "
                          "{particle_step_rate_total:4.2f} particle*steps/s overall"),
            'memdebug': ("{system_memory}\n"
                         "Proc {iproc} usage:\n{memory_usage}"),
        }
        if custom_string is not None:
            self.diag_string = custom_string
        else:
            if preset_string not in self.defaults_dict:
                logger.warning(("Preset {} not found for set_step_diag, "
                                "using default").format(preset_string))
                preset_string = 'default'
            self.diag_string = self.defaults_dict[preset_string]

        super(TextDiag, self).__init__(diag_steps=diag_steps, **kwargs)

        if install:
            callbacks.installafterstep(self.text_diag)

        self.particle_steps_total = 0
        self.previous_particle_steps_total = 0

    def text_diag(self):
        """Write requested information to output."""

        if self.check_timestep():
            live_parts, parts_per_species_str = self._get_part_nums()
            # Loop everything else so run doesn't crash on diag error
            try:
                wall_time = time.time() - self.prev_time
                steps = mwxrun.get_it() - self.prev_step
                # This isn't perfectly accurate, but it's a good approximation
                self.particle_steps_total += live_parts * steps
                if wall_time > 0:
                    particle_step_rate = (
                        self.particle_steps_total
                        - self.previous_particle_steps_total
                    ) / wall_time
                    step_rate = steps / wall_time
                else:
                    step_rate = 0
                    particle_step_rate = 0

                total_elapsed_time = time.time() - self.start_time
                particle_step_rate_total = self.particle_steps_total / total_elapsed_time

                self.status_dict = {
                    'step': mwxrun.get_it(),
                    'nplive': live_parts,
                    'npperspecies': parts_per_species_str,
                    'wall_time': wall_time,
                    'step_rate': step_rate,
                    'particle_step_rate': particle_step_rate,
                    "diag_steps": self.diag_steps,
                    "particle_step_rate_total" : particle_step_rate_total,
                    # TODO: Reimplement when we have new parallel comm hook.
                    # 'iproc': mwxutil.iproc,
                    'iproc': None,
                }

                # Iff memory usage is requested, compute it.
                if (("system_memory" in self.diag_string)
                        or ("memory_usage" in self.diag_string)):
                    self.update_memory()

                # Support child objects by having arbitrary updates
                self._update_status_dict()

                # Print the string with everything that ended up in the
                # dictionary.
                logger.info(self.diag_string.format(**self.status_dict))
                self.previous_particle_steps_total = self.particle_steps_total

            except Exception as err:
                logger.error(
                    f"Failed to output diag_string {self.diag_string} "
                    f"with error {err}"
                )

            self.prev_time = time.time()
            self.prev_step = mwxrun.get_it()

    def _update_status_dict(self):
        """Child objects can override this to make additional modifications to
        status dict.

        Do this by assigning key/value pairs into self.status_dict.
        """

    def _get_part_nums(self):
        """Handle fetching of particle numbers."""
        live_parts = mwxrun.get_npart()
        npart_dict = mwxrun.get_npart_species_dict()

        parts_per_species_str = '[{}]'.format(
            ', '.join(
                "{}: {}".format(key, val) for key, val in npart_dict.items()
            )
        )

        return live_parts, parts_per_species_str

    def update_memory(self):
        """Update memory usage information with psutil."""
        raise NotImplementedError("Need new iproc implementation to use.")
        # # See psutil/scripts/meminfo.py and
        # # https://stackoverflow.com/questions/276052/how-to-get-current-cpu-and-ram-usage-in-python
        # # for some inspiration.
        # if mwxutil.iproc == 0:
        #     sysmem = "SYSTEM MEMORY USAGE\n-------------------\n"
        #     sysmem += self._prettyprint_mem(psutil.virtual_memory())
        #     sysmem += "\nSWAP USAGE\n----------\n"
        #     sysmem += self._prettyprint_mem(psutil.swap_memory())
        #     self.status_dict['system_memory'] = sysmem
        # else:
        #     self.status_dict['system_memory'] = ""

        # # Gets current process w/ no argument
        # proc = psutil.Process()
        # procmem = self._prettyprint_mem(proc.memory_info())
        # self.status_dict['memory_usage'] = procmem

    @staticmethod
    def _prettyprint_mem(namedtuple):
        """Helper function for memory usage printing. From
        psutil/scripts/meminfo.py.

        Returns string.
        """
        raise NotImplementedError("Need new iproc implementation to use.")
        # memstr = ""
        # for name in namedtuple._fields:
        #     value = getattr(namedtuple, name)
        #     if name != 'percent':
        #         value = psutil._common.bytes2human(value)
        #     memstr += '%-10s : %7s\n' % (name.capitalize(), value)
        # return memstr

    def print_performance_summary(self):
        total_time = time.time() - self.start_time
        total_timesteps = mwxrun.get_it()
        steps_per_second = total_timesteps / total_time

        steps_per_second_per_proc = steps_per_second / mwxrun.n_procs

        particle_steps_per_second = self.particle_steps_total / total_time
        particle_steps_per_second_per_proc = (
            particle_steps_per_second / mwxrun.n_procs
        )
        logger.info("### Run Summary ###")

        logger.info(f"steps / second : {steps_per_second:.4f}")
        logger.info(f"steps / second / proc : {steps_per_second_per_proc:.4f}")
        logger.info(f"particle * steps / second : {particle_steps_per_second:.4f}")
        logger.info(f"particle * steps / second / proc : {particle_steps_per_second_per_proc:.4f}")
