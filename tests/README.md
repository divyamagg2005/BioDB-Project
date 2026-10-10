# Tests

Run from the repository root:

```bash
python -m pytest
```

- Unit tests (`test_common.py`, `test_phaseN_*.py` logic tests) use small in-memory fixtures and need no network or data files.
- Tests marked `data` check a phase's generated outputs under `data/processed/`. They are skipped when those outputs do not exist yet.
