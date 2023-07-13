How to contribute
=================

.. important::

    The markup used for the documentation is `reStructuredText`_.

.. _reStructuredText: https://docutils.sourceforge.io/rst.html

We maintain the EEW Simulation Package for Python as an open-source project for common interest. To allow your contributions if provided via a `pull-request`_, you should follow the following rules:

1. Make sure the package can be installed from source (requires to remove the pre-installed version as bellow or to increment version number):
   
    .. code:: sh

        python -m pip uninstall eewsimpy 
        python -m pip install ./  

2. Make sure all dependencies are listed in `requirements.txt`.
3. Make sure version number is upgraded in `setup.py` in order to trigger `pip` update and consistent with `docs/conf.py` (variable `release`).
4. Make sure all features are documented and test documentation with (requires `sphinx` installed) and check the result (open `_build/html/index.html`):
   
    .. code:: sh

        sphinx-apidoc -f -o docs eewsimpy 
        cd docs 
        make clean html

5. Submit your contribution via a pull request from a branch of a fork, see the example below.
6. Please make sure your commits and pull requests are properly described following common best practices such as https://www.atlassian.com/blog/git/written-unwritten-guide-pull-requests.
7. Know that the title of your pull request will be used as the default merge message for all squashed commits (https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/configuring-pull-request-merges/configuring-commit-squashing-for-pull-requests).

Please also be patient with the maintenance team, and remember that your help is very welcomed!

.. _pull-request: https://docs.github.com/en/get-started/quickstart/contributing-to-projects


How to PR from fork branch
--------------------------

1. Create a fork of the main project online at https://github.com/ATTAC-EEW-in-Central-America/EEWsimpy
2. Clone the main project, and go in the clone:

    .. code:: sh

        git clone https://github.com/ATTAC-EEW-in-Central-America/EEWsimpy
        cd EEWsimpy

3. Setup pushing to your fork (e.g., with mine, have your own instead):

    .. code:: sh
    
        git remote set-url --add --push origin https://github.com/<YOUR USER NAME>/EEWsimpy

4. Push your contribution to a branch of your fork:

    .. code:: sh

        git checkout -b <CONTRIB NAME>
        git add <files with contrib>
        git commit 
        git push origin  <CONTRIB NAME>

.. note::

    Note that git provides a direct link to create your pull request, e.g.:

        .. code:: sh       
                         
            ...
            remote: Create a pull request for '<CONTRIB NAME>' on GitHub by visiting:
            remote:      https://github.com/ATTAC-EEW-in-Central-America/eewsimpy

Bonus
-----

Project owners might push minor contribution from their master branch to both their forks and the main project repo: 

.. code:: sh    

    git push all master
