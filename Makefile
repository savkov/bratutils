NAME=bratutils
VERSION=$(shell git rev-parse HEAD)

clean:
	rm -rf .cache
	find . -name "*pyc" -delete
	find . -name ".coverage" -delete

build: clean
	pip install -r test_requirements.txt

test:
	pytest -v --cov-config .coveragerc --cov .
	coverage xml

lint:
	flake8 src/bratutils

release: build
	pip install twine
	python setup.py sdist bdist_wheel
	twine upload -u ${PYPI_USER} -p ${PYPI_PASS} dist/*