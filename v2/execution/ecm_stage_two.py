"""Compatibility import; implementation lives in v2.ecm.stage_two."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("v2.ecm.stage_two")
