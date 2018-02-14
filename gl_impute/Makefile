init:
	pip install -r requirements.txt
	pip install git+https://github.com/bmilius-nmdp/pyglstring.git

test:
	nosetests -v tests

install:
	python setup.py install

upgrade:
	pip install . --upgrade

clean:
	pip uninstall gl-impute
	python setup.py clean --all

examples:
	python -m doctest -v examples.rst