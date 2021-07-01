CFLAGS := -O3 -fPIC -g -std=gnu99 -pedantic -Wall -Wfatal-errors
INCLUDES := -Iinclude
LIBS := -Llib -lhdf5 -lcube -lpumas -lm

.PHONY: build clean

build: build-cython

lib: lib/libgeometry_multi_cube.so
	@rm -f *.o

lib/libgeometry_double_cube.so: src/geometry_multi_cube.c include/geometry_multi_cube.h
	@mkdir -p lib
	@$(CC) -o $@ $(CFLAGS) -shared $(INCLUDES) $(LIBS) $<

build-cython:
	pip install -r requirements_install.txt
	python setup.py build_ext --inplace -i
	pip install -e .

clean: clean-c clean-build clean-pyc clean-cython

clean-c:
	@rm -rf lib bin

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

clean-cython:
	find . -name 'cgeometry_multi_cube.c' -exec rm -f {} +
	find . -name '*.so' -exec rm -f {} +
