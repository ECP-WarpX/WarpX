Unit Tests
==========

Structure
---------
To test functionality within mewarpx, a unit test is written for the new module
using the following format

:code:`test_<module_name>.py`

Within each unit test are functions in the underlying module that check whether
each aspect of the module is working properly.  Each test function follows the
format

:code:`def test_<function_name>():`

and must :code:`assert True` if the test should pass or :code:`assert False` if the test
should fail at the last statement of the function.  The assertions should compare
some result of the function (current, charge density, electric potential etc.) to
a known result located in :file:`tests/test_files/Result_Stats`.  Typically, a
single unit test should not run for too long, 30 seconds to 1 minute or less.

All test files for mewarpx are located in the ``tests`` directory.  For each pull request
into :file:`ModernElectron/WarpX`, CircleCi will execute all tests with pytest to determine
if all tests pass.

Examples
--------
`Test Assemblies <https://github.com/ModernElectron/WarpX/blob/memaster/mewarpx/tests/test_assemblies.py>`_

`Test Emission <https://github.com/ModernElectron/WarpX/blob/memaster/mewarpx/tests/test_emission.py>`_

`Test Flux Diagnostic <https://github.com/ModernElectron/WarpX/blob/memaster/mewarpx/tests/test_flux_diagnostic.py>`_

`Test mwxrun <https://github.com/ModernElectron/WarpX/blob/memaster/mewarpx/tests/test_mwxrun.py>`_