"""Compatibility import; implementation lives in v2.ecm.presets."""

import importlib
import sys

sys.modules[__name__] = importlib.import_module("v2.ecm.presets")
