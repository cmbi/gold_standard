PROJECT = gold-standard
PACKAGE = gold_standard_src/gold_standard
PYTHON_FILES = $(PACKAGE)/*.py

# Build environment
VIRTUALENV = .venv
VENVACTIVATE = $(VIRTUALENV)/bin/activate

# Output
TARGET = target/
SLOCCOUNT = $(TARGET)sloccount.sc
PYLINT = $(TARGET)pylint.txt
PYFLAKES = $(TARGET)pyflakes.txt
PEP8 = $(TARGET)pep8.txt
NOSE = $(TARGET)nosetests.xml
COVERAGE = $(TARGET)coverage.xml
EPYDOC = $(TARGET)epydoc
VULTURE=$(TARGET)vulture.txt
CLONEDIGGER=$(TARGET)clonedigger.xml


all: clean $(SLOCCOUNT) $(PYLINT) $(PYFLAKES) $(PEP8) $(NOSE) $(COVERAGE) $(EPYDOC) $(VULTURE) $(CLONEDIGGER)

clean:
	-rm -r $(TARGET) $(PACKAGE)/*.pyc

$(TARGET): $(BASEDIR)
	-@mkdir -p $@

$(VIRTUALENV): $(VENVACTIVATE)
$(VENVACTIVATE): requirements.txt
	test -d $(VIRTUALENV) || virtualenv --system-site-packages $(VIRTUALENV)
	. $@; \
	pip install -Ur $<
	touch $@

$(SLOCCOUNT): $(PACKAGE) $(TARGET)
	sloccount --duplicates --details $</../ > $@

$(PYLINT): $(PACKAGE) $(VIRTUALENV) $(TARGET)
	-. $(VENVACTIVATE) && pylint --rcfile=.pylintrc $(PYTHON_FILES) > $@

$(PYFLAKES): $(PACKAGE) $(TARGET)
	-. $(VENVACTIVATE) && pyflakes $< > $@

$(PEP8): $(PACKAGE) $(VIRTUALENV) $(TARGET)
	-. $(VENVACTIVATE) && pep8 --max-line-length=150 $< > $@

$(NOSE): $(PACKAGE) $(VIRTUALENV) $(TARGET)
	. $(VENVACTIVATE) && nosetests --traverse-namespace  gold_standard_src/tests --verbosity=3 --with-coverage --cover-xml --with-xunit --xunit-file=$@

$(COVERAGE): $(NOSE) $(VIRTUALENV) $(TARGET)
	-. $(VENVACTIVATE) && coverage xml -o $@

$(EPYDOC): $(PACKAGE) $(VIRTUALENV) $(TARGET)
	-. $(VENVACTIVATE) && epydoc -v --graph all --output $@ --name $(PROJECT) --url https://www.bio-prodict.nl/ $(PYTHON_FILES)

$(CLONEDIGGER): $(PACKAGE) $(VIRTUALENV) $(TARGET)
	-. $(VENVACTIVATE) && clonedigger --fast --ignore-dir=$(VIRTUALENV) --cpd-output --output=$@ $<

$(VULTURE): $(PACKAGE) $(VIRTUALENV) $(TARGET)
	-. $(VENVACTIVATE) && vulture $< > $@
