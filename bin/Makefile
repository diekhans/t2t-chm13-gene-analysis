FLAKE8 = python3 -m flake8


pyprogs = $(shell file -F $$'\t' * | awk '/Python script/{print $$1}')

lint:
	${FLAKE8} ${pyprogs} *.py
