.. _developers-contributing:

Contribute to WarpX
===================

We welcome new contributors!
Here is how to participate to the WarpX development.

Git workflow
------------

The WarpX project uses `git <https://git-scm.com>`_ for version control.
If you are new to git, you can follow one of these tutorials:

- `Learn git with bitbucket <https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud>`_
- `git - the simple guide <http://rogerdudler.github.io/git-guide/>`_

Configure your GitHub Account & Development Machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, let's setup your Git environment and GitHub account.

1. Go to https://github.com/settings/profile and add your real name and affiliation
2. Go to https://github.com/settings/emails and add & verify the professional e-mails you want to be associated with.
3. Configure ``git`` on the machine you develop on to *use the same spelling of your name and email*:

   - ``git config --global user.name "FIRSTNAME LASTNAME"``
   - ``git config --global user.email EMAIL@EXAMPLE.com``
4. Go to https://github.com/settings/keys and add the SSH public key of the machine you develop on. (Check out the GitHub guide to `generating SSH keys <https://docs.github.com/articles/generating-an-ssh-key/>`__ or `troubleshoot common SSH problems <https://docs.github.com/ssh-issues/>`__. )

Make your own fork
^^^^^^^^^^^^^^^^^^

First, fork the WarpX `"mainline" repo on GitHub <https://github.com/ECP-WarpX/WarpX>`__ by pressing the *Fork* button on the top right of the page.
A fork is a copy of WarpX on GitHub, which is under your full control.

Then, we create local copies, for development:

.. code-block:: sh

   # Clone the mainline WarpX source code to your local computer.
   # You cannot write to this repository, but you can read from it.
   git clone git@github.com:ECP-WarpX/WarpX.git
   cd WarpX

   # rename what we just cloned: call it "mainline"
   git remote rename origin mainline

   # Add your own fork. You can get this address on your fork's Github page.
   # Here is where you will publish new developments, so that they can be
   # reviewed and integrated into "mainline" later on.
   # "myGithubUsername" needs to be replaced with your user name on GitHub.
   git remote add myGithubUsername git@github.com:myGithubUsername/WarpX.git

Now you are free to play with your fork (for additional information, you can visit the
`Github fork help page <https://help.github.com/en/articles/fork-a-repo>`__).

.. note::

   We only need to do the above steps for the first time.

Let's Develop
^^^^^^^^^^^^^

You are all set!
Now, the basic WarpX development workflow is:

1. Implement your changes and push them on a new branch ``branch_name`` on your fork.
2. Create a Pull Request from branch ``branch_name`` on your fork to branch ``development`` on the main WarpX repo.

Create a branch ``branch_name`` (the branch name should reflect the piece of code you want to add, like ``fix-spectral-solver``) with

.. code-block:: sh

   # start from an up-to-date development branch
   git checkout development
   git pull mainline development

   # create a fresh branch
   git checkout -b branch_name

and do the coding you want.

It is probably a good time to look at the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/>`_ and at the Doxygen reference pages:

* WarpX Doxygen: https://warpx.readthedocs.io/en/latest/_static/doxyhtml
* AMReX Doxygen: https://amrex-codes.github.io/amrex/doxygen
* PICSAR Doxygen: (todo)

Once you are done developing, add the files you created and/or modified to the ``git`` *staging area* with

.. code-block:: sh

   git add <file_I_created> <and_file_I_modified>


Build your changes
^^^^^^^^^^^^^^^^^^

If you changed C++ files, then now is a good time to test those changes by compiling WarpX locally.
Follow the `developer instructions in our manual <https://warpx.readthedocs.io/en/latest/install/cmake.html>`__ to set up a local development environment, then compile and `run <https://warpx.readthedocs.io/en/latest/usage/how_to_run.html>`__ WarpX.


Commit & push your changes
^^^^^^^^^^^^^^^^^^^^^^^^^^

Periodically commit your changes with

.. code-block:: sh

   git commit

The commit message (between quotation marks) is super important in order to follow the developments during code-review and identify bugs.
A typical format is:

.. code-block:: text

   This is a short, 40-character title

   After a newline, you can write arbitray paragraphs. You
   usually limit the lines to 70 characters, but if you don't, then
   nothing bad will happen.

   The most important part is really that you find a descriptive title
   and add an empty newline after it.

For the moment, commits are on your local repo only.
You can push them to your fork with

.. code-block:: sh

   git push -u myGithubUsername branch_name

If you want to synchronize your branch with the ``development`` branch (this is useful when the ``development`` branch is being modified while you are working on ``branch_name``), you can use

.. code-block:: sh

   git pull mainline development

and fix any conflict that may occur.

Submit a Pull Request
^^^^^^^^^^^^^^^^^^^^^

A Pull Request (PR) is the way to efficiently visualize the changes you made and to propose your new feature/improvement/fix to the WarpX project.
Right after you push changes, a banner should appear on the Github page of your fork, with your ``<branch_name>``.

- Click on the ``compare & pull request`` button to prepare your PR.
- It is time to communicate your changes: write a title and a description for your PR.
  People who review your PR are happy to know

  * what feature/fix you propose, and why
  * how you made it (added new/edited files, created a new class than inherits from...)
  * how you tested it and what was the output you got
  * and anything else relevant to your PR (attach images and scripts, link papers, *etc.*)
- Press ``Create pull request``.
  Now you can navigate through your PR, which highlights the changes you made.

Please DO NOT write large pull requests, as they are very difficult and time-consuming to review.
As much as possible, split them into small, targeted PRs.
For example, if find typos in the documentation open a pull request that only fixes typos.
If you want to fix a bug, make a small pull request that only fixes a bug.

If you want to implement a feature and are not too sure how to split it, just open an issue about your plans and ping other WarpX developers on it to chime in.
Generally, write helper functionality first, test it and then write implementation code.
Submit tests, documentation changes and implementation of a feature together for pull request review.

Even before your work is ready to merge, it can be convenient to create a PR (so you can use Github tools to visualize your changes).
In this case, please put the ``[WIP]`` tag (for Work-In-Progress) at the beginning of the PR title.
You can also use the GitHub project tab in your fork to organize the work into separate tasks/PRs and share it with the WarpX community to get feedback.

Include a test to your PR
"""""""""""""""""""""""""

A new feature is great, a **working** new feature is even better!
Please test your code and add your test to the automated test suite.
It's the way to protect your work from adventurous developers.
Instructions are given in the :ref:`testing section <developers-testing>` of our `developer's documentation <https://warpx.readthedocs.io/en/latest/developers/testing.html>`_.

Include documentation about your PR
"""""""""""""""""""""""""""""""""""

Now, let users know about your new feature by describing its usage in the `WarpX documentation <https://warpx.readthedocs.io>`_.
Our documentation uses `Sphinx <http://www.sphinx-doc.org/en/master/usage/quickstart.html>`_, and it is located in ``Docs/source/``.
For instance, if you introduce a new runtime parameter in the input file, you can add it to :ref:`Docs/source/running_cpp/parameters.rst <running-cpp-parameters>`.
If Sphinx is installed on your computer, you should be able to generate the html documentation with

.. code-block:: sh

   make html

in ``Docs/``. Then open ``Docs/build/html/index.html`` with your favorite web browser and look
for your changes.

Once your code is ready with documentation and automated test, congratulations!
You can create the PR (or remove the ``[WIP]`` tag if you already created it).
Reviewers will interact with you if they have comments/questions.


.. _developers-contributing-style-conventions:

Style and conventions
---------------------

- For indentation, WarpX uses four spaces (no tabs)

- Some text editors automatically modify the files you open. We recommend to turn on to remove trailing spaces and replace Tabs with 4 spaces.

- The number of characters per line should be <100

- Exception: in documentation files (``.rst``/``.md``) use one sentence per line independent of its number of characters, which will allow easier edits.

- Space before and after assignment operator (``=``)

- To define a function , for e.g., ``myfunction()`` use a space between the name of the function and the paranthesis - ``myfunction ()``.
  To call the function, the space is not required, i.e., just use ``myfunction()``.

- The reason this is beneficial is that when we do a ``git grep`` to search for ``myfunction ()``, we can clearly see the locations where ``myfunction ()`` is defined and where ``myfunction()`` is called.

- Also, using ``git grep "myfunction ()"`` searches for files only in the git repo, which is more efficient compared to the ``grep "myfunction ()"`` command that searches through all the files in a directory, including plotfiles for example.

- It is recommended that style changes are not included in the PR where new code is added.
  This is to avoid any errors that may be introduced in a PR just to do style change.

- WarpX uses ``CamelCase`` convention for file names and class names, rather than ``snake_case``.

- The names of all member variables should be prefixed with ``m_``.
  This is particularly useful to avoid capturing member variables by value in a lambda function, which causes the whole object to be copied to GPU when running on a GPU-accelerated architecture.
  This convention should be used for all new piece of code, and it should be applied progressively to old code.

- ``#include`` directives in C++ have a distinct order to avoid bugs, see :ref:`the WarpX repo structure <developers-repo-structure>` for details

- For all new code, we should avoid relying on ``using namespace amrex;`` and all amrex types should be prefixed with `amrex::`.
  Inside limited scopes, AMReX type literals can be included with ``using namespace amrex::literals;``.
  Ideally, old code should be modified accordingly.
