# Factor project commands use PyPy unless PYTHON is explicitly overridden.
PYTHON ?= pypy3
N ?=
SEED ?= 7
REPETITIONS ?= 9
WARMUP_SECONDS ?= 3
MILESTONE ?= phase_1
RUN_ID := $(shell date +%Y%m%dT%H%M%S)
RUNTIME_NAME := $(notdir $(PYTHON))
VALIDATION_OUTPUT ?= audit/results/validation_$(RUN_ID)_$(RUNTIME_NAME).json
BENCHMARK_OUTPUT ?= benchmarks/results/runtime_$(RUN_ID)_$(RUNTIME_NAME).json

# Resolve v2 imports even when make is invoked with -C from another directory.
MAKEFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
REPO_ROOT := $(abspath $(MAKEFILE_DIR)/..)
export PYTHONPATH := $(REPO_ROOT)$(if $(PYTHONPATH),:$(PYTHONPATH))
export PYTHONDONTWRITEBYTECODE := 1

.PHONY: help runtime run run-bounded test validate benchmark benchmark-phase-two lint
.PHONY: benchmark-phase-three-reference benchmark-phase-three-collector
.PHONY: benchmark-phase-three-pipeline benchmark-phase-three-families
.PHONY: benchmark-phase-three-siqs benchmark-phase-three-large
.PHONY: benchmark-phase-three-sss
.PHONY: benchmark-phase-three-parallel
.PHONY: benchmark-backends benchmark-backend-sizes
BACKEND_SUITES ?= arithmetic ecm qs
BACKEND_SPLIT ?= screen
LARGE_PHASE ?= all

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
	@mkdir -p "$(dir $(VALIDATION_OUTPUT))"
	$(PYTHON) -m v2.benchmarks.validate_phase_one \
	  --milestone "$(MILESTONE)" \
	  --output "$(VALIDATION_OUTPUT)"

benchmark:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_one \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-two:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_two \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-reference:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_reference \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-collector:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_collector \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-pipeline:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_pipeline \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-families:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_families \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-siqs:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_siqs \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-large:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_large \
	  --phase $(LARGE_PHASE) --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-sss:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_sss \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-phase-three-parallel:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.phase_three_parallel \
	  --phase all --repetitions $(REPETITIONS) \
	  --warmup-seconds $(WARMUP_SECONDS) --output "$(BENCHMARK_OUTPUT)"

benchmark-backends:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.p43_backend \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

benchmark-backend-sizes:
	@mkdir -p "$(dir $(BENCHMARK_OUTPUT))"
	$(PYTHON) -u -m v2.benchmarks.p43_sizes \
	  --suites $(BACKEND_SUITES) --split $(BACKEND_SPLIT) \
	  --repetitions $(REPETITIONS) --warmup-seconds $(WARMUP_SECONDS) \
	  --output "$(BENCHMARK_OUTPUT)"

# Tools are installed in the ignored, project-local PyPy development venv.
lint:
	.venv/bin/ruff check --no-cache .
	.venv/bin/ruff format --check --no-cache .
	.venv/bin/python -m pycodestyle --exclude=audit,.venv .
