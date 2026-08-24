.PHONY: all extensions clean test

all: extensions

extensions:
	python tools/build_extensions.py

clean:
	python tools/build_extensions.py --clean
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	rm -rf build dist *.egg-info

test:
	python -m unittest discover -s tests -v
