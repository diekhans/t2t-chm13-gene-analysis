PYTHON = python3
pyprogs = $(shell file -F $$'\t' bin/* tests/bin tests/*/bin/* tests/*/*/bin/* | awk '/Python script/{print $$1}')

all::

lint:
	${PYTHON} -m flake8 ${pyprogs} bin/geneBoundsLib.py

