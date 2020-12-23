.. _developers-workflow:

Workflow
========

Create a new Github release
---------------------------

WarpX has one release per month.
In order to create a release, you need to:

 1. Create a new branch from ``development`` and update the version number in all source files.
    There is a script for that, so you can do:

    .. code-block:: sh

       cd Tools/DevUtils/
       ./update_release.sh # This replaces the old version number with the new one.

    Then open a PR, as usual. NOTE: do not merge this PR before step 2 is completed.

 2. Click the ``Draft a new release`` button at https://github.com/ECP-WarpX/WarpX/releases and follow instructions.
    Please specify the compatible versions of dependencies (see previous releases), and provide info on the content of the release.
    In order to get a list of PRs merged since last release, you may run

    .. code-block:: sh

       git log --since=<date> | grep -A 3 "Author: " | grep -B 1 "\-\-" | sed '/--/d' | sed -e 's/^    /- /'

    where ``<date>`` is the date of the last release, say ``2020-05-01`` if the last release was on May 1, 2020.

 3. Optional: create a ``release-<version>`` branch, write a changelog, and backport bug-fixes for a few days.
