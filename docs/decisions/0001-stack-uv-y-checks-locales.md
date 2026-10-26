# 0001. Stack declarado para uv, checks locales con PYTHONPATH=src
- Date: 2026-10-26
- Status: accepted
- Decision: pyproject.toml declara requires-python ">=3.12", dependencias numpy, scipy y matplotlib, pytest en [dependency-groups] dev, backend hatchling con src layout y [tool.pytest.ini_options] pythonpath = ["src"]. Todos los checks de aceptación corren con el Python local como `PYTHONPATH=src python ...`, sin instalar nada (ni uv ni pip install -e). La FFT es numpy.fft.
- Reason: en la máquina no hay uv y está prohibido instalar o descargar; el Python 3.14 local ya trae numpy 2.4.4, scipy 1.18.0, matplotlib y pytest. El código no usa nada posterior a 3.12, así que funciona igual con uv.
