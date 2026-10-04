# Factor project commands use PyPy unless PYTHON is explicitly overridden.
PYTHON ?= pypy3
N ?=
SEED ?= 7
REPETITIONS ?= 9
WARMUP_SECONDS ?= 3
MILESTONE ?= phase_1
RUN_ID := $(shell date +%Y%m%dT%H%M%S)
RUNTIME_NAME := $(notdir $(PYTHON))
VALIDATION_OUTPUT ?= audit/validation_$(RUN_ID)_$(RUNTIME_NAME).json
BENCHMARK_OUTPUT ?= benchmarks/runtime_$(RUN_ID)_$(RUNTIME_NAME).json

# Resolve v2 imports even when make is invoked with -C from another directory.
MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
REPO_ROOT := $(abspath $(MAKEFILE_DIR)/..)
export PYTHONPATH := $(REPO_ROOT)$(if $(PYTHONPATH),:$(PYTHONPATH))
export PYTHONDONTWRITEBYTECODE := 1

.PHONY: help runtime run run-bounded test validate benchmark benchmark-phase-two lint
.PHONY: benchmark-phase-three-reference benchmark-phase-three-collector
.PHONY: benchmark-phase-three-pipeline

help:
	@echo "Default runtime: $(PYTHON)"
	@echo "make -C v2 run N=626100403 SEED=7"
	@echo "make -C v2 test | validate | benchmark | runtime | lint"
	@echo "Supported runtime: PyPy implementing Python 3.11."

runtime:
	$(PYTHON) --version

run:
	$(PYTHON) factor.py $(N) --seed $(SEED)

run-bounded:
	$(PYTHON) factor.py $(N) --seed $(SEED) --bounded

test:
	$(PYTHON) -m unittest discover -s tests -v

validate:
	$(PYTHON) -m v2.benchmarks.validate_phase_one \
	  --milestone "$(MILESTONE)" \
	  --output "$(VALIDATION_OUTPUT)"

benchmark:
	$(PYTHON) -u -m v2.benchmarks.phase_one \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-two:
	$(PYTHON) -u -m v2.benchmarks.phase_two \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-reference:
	$(PYTHON) -u -m v2.benchmarks.phase_three_reference \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-collector:
	$(PYTHON) -u -m v2.benchmarks.phase_three_collector \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-pipeline:
	$(PYTHON) -u -m v2.benchmarks.phase_three_pipeline \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

# Tools are installed in the ignored, project-local PyPy development venv.
lint:
	.venv/bin/ruff check --no-cache .
	.venv/bin/ruff format --check --no-cache .
	.venv/bin/python -m pycodestyle --exclude=audit,.venv .
