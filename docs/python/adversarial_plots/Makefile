PYTHON?=python3
VENV?=venv
REQUIREMENTS?=requirements.txt

.PHONY: install reinstall

install:
	$(PYTHON) -m pip install virtualenv
	$(PYTHON) -m virtualenv $(VENV)
	. $(VENV)/bin/activate; $(PYTHON) -m pip install -r $(REQUIREMENTS)

reinstall:
	rm -rf $(VENV)
	$(MAKE) install

notebook:
	. $(VENV)/bin/activate; /home/vesely/.local/bin/jupyter notebook --notebook-dir=notebooks
