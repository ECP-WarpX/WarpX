# Copyright 2018-2019 Andrew Myers, Luca Fedeli, Maxence Thevenet
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script modifies `WarpX-test.ini` (which is used for nightly builds)
# and creates the file `travis-test.ini` (which is used for continous
# integration on TravisCI (https://travis-ci.org/)
# The subtests that are selected are controlled by WARPX_TEST_DIM
# The architecture (CPU/GPU) is selected by WARPX_TEST_ARCH
import re
import os
# Get relevant environment variables
dim = os.environ.get('WARPX_CI_DIM', None)
qed = os.environ.get('WARPX_CI_QED', 'FALSE')
arch = os.environ.get('WARPX_TEST_ARCH', 'CPU')
single_precision = os.environ.get('WARPX_CI_SINGLE_PRECISION', 'FALSE')
electrostatic = os.environ.get('WARPX_CI_ELECTROSTATIC', 'FALSE')
python_main = os.environ.get('WARPX_CI_PYTHON_MAIN', 'FALSE')
psatd = os.environ.get('WARPX_CI_PSATD', 'FALSE')
# Find the directory in which the tests should be run
current_dir = os.getcwd()
test_dir = re.sub('warpx/Regression', '', current_dir )

with open('WarpX-tests.ini') as f:
    text = f.read()

# Replace default folder name
text = re.sub('/home/regtester/AMReX_RegTesting', test_dir, text)
# Remove the web directory
text = re.sub('[\w\-\/]*/web', '', text)

# Add doComparison = 0 for each test
text = re.sub( '\[(?P<name>.*)\]\nbuildDir = ',
               '[\g<name>]\ndoComparison = 0\nbuildDir = ', text )

# Change compile options when running on GPU
if arch == 'GPU':
    text = re.sub( 'addToCompileString =',
                    'addToCompileString = USE_GPU=TRUE USE_OMP=FALSE USE_ACC=TRUE ', text)
    text = re.sub( 'COMP\s*=.*', 'COMP = pgi', text )
print('Compiling for %s' %arch)

# Add runtime option: crash for unused variables
text = re.sub('runtime_params =',
              'runtime_params = amrex.abort_on_unused_inputs=1 ',
              text)

# Use only 2 cores for compiling
text = re.sub( 'numMakeJobs = \d+', 'numMakeJobs = 2', text )
# Use only 1 MPI and 1 thread proc for tests
text = re.sub( 'numprocs = \d+', 'numprocs = 2', text)
text = re.sub( 'numthreads = \d+', 'numthreads = 1', text)
# Prevent emails from being sent
text = re.sub( 'sendEmailWhenFail = 1', 'sendEmailWhenFail = 0', text )

# Remove Python test (does not compile)
text = re.sub( '\[Python_Langmuir\]\n(.+\n)*', '', text)

# Remove Langmuir_x/y/z test (too long; not that useful)
text = re.sub( '\[Langmuir_[xyz]\]\n(.+\n)*', '', text)

# Select the tests to be run
# --------------------------

# - Extract test blocks (they are identified by the fact that they contain "inputFile")
select_test_regex = r'(\[(.+\n)*inputFile(.+\n)*)'
test_blocks =  [ match[0] for match in re.findall(select_test_regex, text) ]
# - Remove the test blocks from `text` (only the selected ones will be added back)
text = re.sub( select_test_regex, '', text )

# Keep tests that have the right dimension
if dim is not None:
    print('Selecting tests with dim = %s' %dim)
    # Cartesian tests
    if dim in ['2', '3']:
        test_blocks = [ block for block in test_blocks \
             if ('dim = %s'%dim in block) and not ('USE_RZ' in block) ]
    elif dim == 'RZ':
        test_blocks = [ block for block in test_blocks if 'USE_RZ' in block ]
    else:
        raise ValueError('Unkown dimension: %s' %dim)

def select_tests(blocks, do_test, test_env_variable, match_string_list):
    # Remove or keep tests according to do_test variable
    if do_test not in ['TRUE', 'FALSE']:
        raise ValueError(test_env_variable + ' must be TRUE or FALSE')
    if (do_test == "FALSE"):
        for match_string in match_string_list:
            print('Selecting tests without ' + match_string)
            blocks = [ block for block in blocks if not match_string in block ]
    else:
        for match_string in match_string_list:
            print('Selecting tests with ' + match_string)
            blocks = [ block for block in blocks if match_string in block ]
    return blocks

test_blocks = select_tests(test_blocks, qed, 'WARPX_CI_QED',
                           ['USE_QED=TRUE'])
test_blocks = select_tests(test_blocks, single_precision, 'WARPX_CI_SINGLE_PRECISION',
                           ['PRECISION=FLOAT', 'USE_SINGLE_PRECISION_PARTICLES=TRUE'])
test_blocks = select_tests(test_blocks, electrostatic, 'WARPX_CI_ELECTROSTATIC',
                           ['DO_ELECTROSTATIC=TRUE'])
test_blocks = select_tests(test_blocks, python_main, 'WARPX_CI_PYTHON_MAIN',
                           ['PYTHON_MAIN=TRUE'])
test_blocks = select_tests(test_blocks, psatd, 'WARPX_CI_PSATD',
                           ['USE_PSATD=TRUE'])

# - Add the selected test blocks to the text
text = text + '\n' + '\n'.join(test_blocks)

with open('travis-tests.ini', 'w') as f:
    f.write(text)
