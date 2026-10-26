# Reglas del proyecto GBQuantum (Fase 0)

## Reglas permanentes

1. Unidades naturales: hbar = m = 1. Todo es `float64` / `complex128`.
2. Nunca declares algo "correcto" sin un test que pase.
3. Cada tolerancia numérica lleva un comentario que la justifica. Si un check falla, no se relaja la tolerancia: se reporta la causa.
4. Sin features extra: lo que sobra va a "Mejoras no implementadas" en `STATUS.md`.
5. Al terminar, `STATUS.md` incluye una sección "Lo que NO verifiqué".
6. No gastes más del 10% del esfuerzo en documentación; docstrings y comentarios en español, identificadores en inglés.

## Reparto de dueños

- ROL A: `src/quantum1d/`, `experiments/`, `outputs/`, `pyproject.toml`, `STATUS.md`.
- ROL B: `tests/` y `docs/verification.md` (verifica a ciegas contra la API).
- Arquitecto: `docs/physics-spec.md`, `docs/plan.md`, `docs/decisions/`.

## Cómo correr

Desde la raíz del repositorio, sin instalar nada:

```
PYTHONPATH=src python experiments/01_tunneling.py
PYTHONPATH=src python experiments/02_euler_vs_splitstep.py
PYTHONPATH=src python -m pytest
```

(PowerShell: `$env:PYTHONPATH='src'; python ...`.)

## Referencias

- Física y API vinculante: `docs/physics-spec.md` (sección 12).
- Decisiones: `docs/decisions/`.
