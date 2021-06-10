.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests with default options"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	# @echo "install - install the package to the active Python's site-packages"
	@echo "devel -- install as a user development package"
	@echo "devel-all -- install as a user development package with all optional components"
	@echo "devel-system -- install as a development package while inside a virtual environment"
	@echo "devel-all-system -- install as a development package with all optional components while inside a virtual environment"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -rf tests/temp_files/*

lint:
	flake8 --max-line-length=79 mewarpx tests

test:
	pytest

coverage:
	pytest --cov=mewarpx --runall tests/
	coverage report -m
	coverage html
	@echo "Now open htmlcov/index.html in your browser"

docs:
	rm -f docs/mewarpx.*
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ mewarpx
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	# $(MAKE) -C docs latexpdf
	@echo "Now open docs/_build/html/index.html in your browser"

# install: clean
# 	python setup.py install

devel: clean
	pip install --user -e "."

devel-all: clean
	pip install --user -e ".[complete]"

devel-system: clean
	pip install -e "."

devel-all-system: clean
	pip install -e ".[complete]"
