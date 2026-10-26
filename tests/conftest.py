"""Acceso diferido exclusivamente a la API pública, sin leer implementación."""

import importlib

import pytest


@pytest.fixture(scope="session")
def q():
    """La ausencia del motor falla explícitamente; nunca se marca como skip."""
    return importlib.import_module("quantum1d")
