CFLAGS := -O3 -fPIC -g -std=gnu99 -pedantic -Wall -Wfatal-errors
LIBS := -Llib -lhdf5 -lcube -lpumas -lm

.PHONY: build clean

build: build-c build-cython

build-c:
	@$(CC) PUMAS_physics_dump.c -o PUMAS_physics_dump $(CFLAGS) $(LIBS)
	./PUMAS_physics_dump
	@$(CC) geometry_double_cube.c -o geometry_double_cube $(CFLAGS) $(LIBS)

build-cython:
	pip install -r requirements.txt
	python setup.py build_ext --inplace -i

clean: clean-c clean-build clean-pyc clean-cython

clean-c:
	rm -rf geometry_double_cube

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

clean-cython:
	find . -name 'cgeometry_double_cube.c' -exec rm -f {} +
	find . -name '*.so' -exec rm -f {} +
