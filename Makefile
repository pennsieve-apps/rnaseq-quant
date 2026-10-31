.PHONY: help test fixture run clean

SERVICE_NAME ?= rnaseq-quant
PYTHON ?= python3

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test     - pytest on the synthetic fixture (no Docker); the end-to-end tests need salmon and fastp on PATH"
	@echo "                (or SALMON_BIN / FASTP_BIN) and are skipped without"
	@echo "make fixture  - write the synthetic input to data/fixture/"
	@echo "make run      - build the image and run it in Docker over data/fixture/, writing to data/output/"
	@echo "make clean    - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	rm -rf data/fixture
	$(PYTHON) -m tests.make_fixture data/fixture

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output
