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
dim = os.environ.get('WARPX_TEST_DIM', None)
qed = os.environ.get('HAS_QED', 'FALSE')
arch = os.environ.get('WARPX_TEST_ARCH', 'CPU')
single_precision = os.environ.get('SINGLE_PRECISION', 'FALSE')
electrostatic = os.environ.get('ELECTROSTATIC', 'FALSE')
python_main = os.environ.get('PYTHON_MAIN', 'FALSE')
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

# Remove or keep QED tests according to 'qed' variable
if qed not in ['TRUE', 'FALSE']:
    raise ValueError('HAS_QED must be TRUE or FALSE')
print('Selecting tests with HAS_QED = %s' %qed)
if (qed == "FALSE"):
    test_blocks = [ block for block in test_blocks if not 'QED=TRUE' in block ]
else:
    test_blocks = [ block for block in test_blocks if 'QED=TRUE' in block ]

# Remove or keep SINGLE_PRECISION tests according to 'single_precision' variable
if single_precision not in ['TRUE', 'FALSE']:
    raise ValueError('SINGLE_PRECISION must be TRUE or FALSE')
print('Selecting tests with SINGLE_PRECISION = %s' %qed)
if (qed == "FALSE"):
    test_blocks = [ block for block in test_blocks if not (
            'PRECISION=FLOAT' in block and
            'USE_SINGLE_PRECISION_PARTICLES=TRUE' in block) ]
else:
    test_blocks = [ block for block in test_blocks if (
            'PRECISION=FLOAT' in block and
            'USE_SINGLE_PRECISION_PARTICLES=TRUE' in block) ]

# Remove or keep ELECTROSTATIC tests according to 'electrostatic' variable
if electrostatic not in ['TRUE', 'FALSE']:
    raise ValueError('ELECTROSTATIC must be TRUE or FALSE')
print('Selecting tests with ELECTROSTATIC = %s' %qed)
if (electrostatic == "FALSE"):
    test_blocks = [ block for block in test_blocks if not 'DO_ELECTROSTATIC=TRUE' in block ]
else:
    test_blocks = [ block for block in test_blocks if 'DO_ELECTROSTATIC=TRUE' in block ]

# Remove or keep PYTHON_MAIN tests according to 'python_main' variable
if python_main not in ['TRUE', 'FALSE']:
    raise ValueError('PYTHON_MAIN must be TRUE or FALSE')
print('Selecting tests with PYTHON_MAIN = %s' %qed)
if (python_main == "FALSE"):
    test_blocks = [ block for block in test_blocks if not 'PYTHON_MAIN=TRUE' in block ]
else:
    test_blocks = [ block for block in test_blocks if 'PYTHON_MAIN=TRUE' in block ]

# - Add the selected test blocks to the text
text = text + '\n' + '\n'.join(test_blocks)

with open('travis-tests.ini', 'w') as f:
    f.write(text)
