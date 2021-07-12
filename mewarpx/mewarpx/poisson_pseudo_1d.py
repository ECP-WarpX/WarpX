import numpy as np
from scipy.sparse import csc_matrix, linalg as sla

from mewarpx.mwxrun import mwxrun
from pywarpx import callbacks
from pywarpx.picmi import Cartesian2DGrid, ElectrostaticSolver, constants

class PoissonSolverPseudo1D(ElectrostaticSolver):

    def __init__(self, grid, **kwargs):
        """Direct solver for the Poisson equation using superLU. This solver is
        useful for pseudo 1D cases i.e. diode simulations with small x extent.

        Arguments:
            grid (picmi.Cartesian2DGrid): Instance of the grid on which the
            solver will be installed.
        """
        # Sanity check that this solver is appropriate to use
        if not isinstance(grid, Cartesian2DGrid):
            raise RuntimeError('Direct solver can only be used on a 2D grid.')

        super(PoissonSolverPseudo1D, self).__init__(
            grid=grid, method=kwargs.pop('method', 'Multigrid'),
            required_precision=1, **kwargs
        )

    def initialize_inputs(self):
        """Grab geometrical quantities from the grid. The boudary potentials
        are also obtained from the grid using 'warpx_potential_zmin' for the
        left_voltage and 'warpx_potential_zmax' for the right_voltage.
        These can be given as floats or callable functions that return the
        potential. If callable that function should accept simulation time
        as an argument."""
        # grab the boundary potentials from the grid object
        left_voltage = self.grid.potential_zmin
        if left_voltage is None:
            left_voltage = 0.0
        if isinstance(left_voltage, str):
            raise AttributeError(
                "The direct solver requires a float or a callable as argument "
                "for 'warpx_potential_zmin'."
            )
        if callable(left_voltage):
            self.left_voltage = left_voltage
        else:
            self.left_voltage = lambda x: left_voltage

        right_voltage = self.grid.potential_zmax
        if right_voltage is None:
            right_voltage = 0.0
        if isinstance(right_voltage, str):
            raise AttributeError(
                "The direct solver requires a float or a callable as argument "
                "for 'warpx_potential_zmax'."
            )
        if callable(right_voltage):
            self.right_voltage = right_voltage
        else:
            self.right_voltage = lambda x: right_voltage

        # set WarpX boundary potentials to None since we will handle it
        # ourselves in this solver
        self.grid.potential_xmin = None
        self.grid.potential_xmax = None
        self.grid.potential_ymin = None
        self.grid.potential_ymax = None
        self.grid.potential_zmin = None
        self.grid.potential_zmax = None

        super(PoissonSolverPseudo1D, self).initialize_inputs()

        self.nx = self.grid.nx
        self.nz = self.grid.ny
        self.dx = (self.grid.xmax - self.grid.xmin) / self.nx
        self.dz = (self.grid.ymax - self.grid.ymin) / self.nz

        if not np.isclose(self.dx, self.dz):
            raise RuntimeError('Direct solver requires dx = dz.')

        self.nxguardrho = 2
        self.nzguardrho = 2
        self.nxguardphi = 1
        self.nzguardphi = 1

        self.phi = np.zeros(
            (self.nx + 1 + 2*self.nxguardphi,
            self.nz + 1 + 2*self.nzguardphi)
        )

        self.decompose_matrix()

        print('Using direct solver.')
        callbacks.installfieldsolver(self._run_solve)

    def decompose_matrix(self):
        """Function to build the superLU object used to solve the linear
        system."""
        self.nzsolve = self.nz + 1
        self.nxsolve = self.nx + 3

        # Set up the computation matrix in order to solve A*phi = rho
        A = np.zeros(
            (self.nzsolve*self.nxsolve, self.nzsolve*self.nxsolve)
        )
        kk = 0
        for ii in range(self.nxsolve):
            for jj in range(self.nzsolve):
                temp = np.zeros((self.nxsolve, self.nzsolve))

                if jj == 0 or jj == self.nzsolve - 1:
                    temp[ii, jj] = 1.
                elif jj == 1:
                    temp[ii, jj] = -2.0
                    temp[ii, jj-1] = 1.0
                    temp[ii, jj+1] = 1.0
                elif jj == self.nzsolve - 2:
                    temp[ii, jj] = -2.0
                    temp[ii, jj+1] = 1.0
                    temp[ii, jj-1] = 1.0
                elif ii == 0:
                    temp[ii, jj] = 1.0
                    temp[-3, jj] = -1.0
                elif ii == self.nxsolve - 1:
                    temp[ii, jj] = 1.0
                    temp[2, jj] = -1.0
                else:
                    temp[ii, jj] = -4.0
                    temp[ii+1, jj] = 1.0
                    temp[ii-1, jj] = 1.0
                    temp[ii, jj-1] = 1.0
                    temp[ii, jj+1] = 1.0

                A[kk] = temp.flatten()
                kk += 1

        A = csc_matrix(A, dtype=np.float32)
        self.lu = sla.splu(A)

    def _run_solve(self):
        """Function run on every step to perform the required steps to solve
        Poisson's equation."""
        # get rho from WarpX
        self.rho_data = mwxrun.get_gathered_rho_grid()
        # run superLU solver only on the root processor
        if mwxrun.me == 0:
            self.rho_data = self.rho_data[0][:,:,0]
            self.solve()
        # write phi to WarpX
        mwxrun.set_phi_grid(self.phi)

    def solve(self):
        """The solution step. Includes getting the boundary potentials and
        calculating phi from rho."""

        left_voltage = self.left_voltage(mwxrun.get_t())
        right_voltage = self.right_voltage(mwxrun.get_t())

        rho = -self.rho_data[
            self.nxguardrho:-self.nxguardrho, self.nzguardrho:-self.nzguardrho
        ] / constants.ep0

        # Construct b vector
        nx, nz = np.shape(rho)
        source = np.zeros((nx+2, nz), dtype=np.float32)
        source[1:-1,:] = rho * self.dx**2

        source[:,0] = left_voltage
        source[:,-1] = right_voltage

        # Construct b vector
        b = source.flatten()

        flat_phi = self.lu.solve(b)
        self.phi[:, self.nzguardphi:-self.nzguardphi] = (
            flat_phi.reshape(np.shape(source))
        )

        self.phi[:,:self.nzguardphi] = left_voltage
        self.phi[:,-self.nzguardphi:] = right_voltage

        # the electrostatic solver in WarpX keeps the ghost cell values as 0
        self.phi[:self.nxguardphi,:] = 0
        self.phi[-self.nxguardphi:,:] = 0
