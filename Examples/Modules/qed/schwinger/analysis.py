import os
import glob
import yt
import numpy as np

# define some parameters

c = 299792458.
m_e = 9.10938356e-31;
e = 1.6021766208e-19;
hbar = 1.054571800e-34;
E_S = m_e**2*c**3/e/hbar

numcells = 512 # total number of cells in the test case
dV = (1.e-6)**3 # total volume
dt = 2.407291502e-16
filename = "diags/plotfiles/plt00001"

Ex_test1 = 1.e16
Ey_test1 = 0.
Ez_test1 = 0.
Bx_test1 = 16792888.570516706
By_test1 = 5256650.141557486
Bz_test1 = 18363530.799561853

Ex_test2 = 1.e18
Ey_test2 = 0.
Ez_test2 = 0.
Bx_test2 = 1679288857.0516706
By_test2 = 525665014.1557486
Bz_test2 = 1836353079.9561853

Ex_test3 = 0.
Ey_test3 = 1.0321239524474501e+17
Ez_test3 = 0.
Bx_test3 = 0. 
By_test3 = 0.
Bz_test3 = 0.

Ex_test4 = 0.
Ey_test4 = 0.
Ez_test4 = 2.5e+20
Bx_test4 = 0. 
By_test4 = 833910154604.3563
Bz_test4 = 0.

def calculate_rate(Ex,Ey,Ez,Bx,By,Bz):

    E_squared = Ex**2 + Ey**2 + Ez**2
    H_squared = c**2*(Bx**2 + By**2 + Bz**2)

    F = (E_squared - H_squared)/2.
    G = c*(Ex*Bx + Ey*By + Ez*Bz)

    epsilon = np.sqrt(np.sqrt(F**2+G**2)+F)/E_S
    eta = np.sqrt(np.sqrt(F**2+G**2)-F)/E_S

    if(epsilon != 0. and eta != 0.):
        return  e**2*E_S**2/4./np.pi**2/c/hbar**2*epsilon*eta/np.tanh(np.pi*eta/epsilon)*np.exp(-np.pi/epsilon)
    elif (epsilon == 0.):
        return 0.
    else:
        return  e**2*E_S**2/4./np.pi**2/c/hbar**2*epsilon**2/np.pi*np.exp(-np.pi/epsilon)


def do_analysis(Ex,Ey,Ez,Bx,By,Bz):

    expected_total_physical_pairs_created = dV*dt*calculate_rate(Ex,Ey,Ez,Bx,By,Bz)

    if expected_total_physical_pairs_created < 0.01:
        assert(not os.path.isdir(filename+"/ele_schwinger/Level_0/"))    
        ## Assert whether pairs are created or not. Is there a cleaner way to do that?
    
    else:
        data_set = yt.load(filename)

        all_data = data_set.all_data()

        ele_data = all_data["ele_schwinger",'particle_weight']
        pos_data = all_data["pos_schwinger",'particle_weight']

        std_total_physical_pairs_created = np.sqrt(expected_total_physical_pairs_created)
        assert(np.array_equal(ele_data,pos_data))
        assert(np.abs(np.sum(ele_data)-expected_total_physical_pairs_created)<5*std_total_physical_pairs_created)


def launch_analysis(executable):
    # First test with "weak" EM field. No pair should be created.
    os.system("./" + executable + " inputs_3d_schwinger 'warpx.E_external_grid = %.15f %.15f %.15f' 'warpx.B_external_grid = %.15f %.15f %.15f'" \
               % (Ex_test1, Ey_test1, Ez_test1, Bx_test1, By_test1, Bz_test1) )
    do_analysis(Ex_test1, Ey_test1, Ez_test1, Bx_test1, By_test1, Bz_test1)

    # Second test with stronger EM field. Many pairs are created and a Gaussian 
    # distribution is used to get the weights of the particles.
    os.system("./" + executable + " inputs_3d_schwinger 'warpx.E_external_grid = %.15f %.15f %.15f' 'warpx.B_external_grid = %.15f %.15f %.15f'" \
               % (Ex_test2, Ey_test2, Ez_test2, Bx_test2, By_test2, Bz_test2) )
    do_analysis(Ex_test2, Ey_test2, Ez_test2, Bx_test2, By_test2, Bz_test2)

    # Third test with intermediate electric field such that average created pair per cell is 1. A 
    # Poisson distribution is used to obtain the weights of the particles.
    os.system("./" + executable + " inputs_3d_schwinger 'warpx.E_external_grid = %.15f %.15f %.15f' 'warpx.B_external_grid = %.15f %.15f %.15f'" \
               % (Ex_test3, Ey_test3, Ez_test3, Bx_test3, By_test3, Bz_test3) )
    do_analysis(Ex_test3, Ey_test3, Ez_test3, Bx_test3, By_test3, Bz_test3)

    # Fourth test with extremely strong EM field but with E and B perpendicular and nearly equal so 
    # that the pair production rate is fairly low. A poisson distribution is used in this case.
    os.system("./" + executable + " inputs_3d_schwinger 'warpx.E_external_grid = %.15f %.15f %.15f' 'warpx.B_external_grid = %.15f %.15f %.15f'" \
               % (Ex_test4, Ey_test4, Ez_test4, Bx_test4, By_test4, Bz_test4) )
    do_analysis(Ex_test4, Ey_test4, Ez_test4, Bx_test4, By_test4, Bz_test4)


def main() :
    executables = glob.glob("main3d*")
    if len(executables) == 1 :
        launch_analysis(executables[0])
    else :
        assert(False)
    print('Passed')

if __name__ == "__main__":
    main()
