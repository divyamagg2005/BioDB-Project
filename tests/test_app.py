from pathlib import Path

import pytest

from biodb import app_data, config

pytestmark = pytest.mark.skipif(
    not (config.METRICS_DIR / "phase7" / "baseline_oof_predictions.csv").exists(), reason="Phase 7 outputs missing"
)
APP = str(Path(__file__).parent.parent / "streamlit_app.py")


def test_parse_variant():
    assert app_data.parse_variant("TP53", "Y220C") == "TP53:p.Y220C"
    assert app_data.parse_variant("TP53", "p.y220c") == "TP53:p.Y220C"
    assert app_data.parse_variant("TP53", "R175H") is None  # not in the cohort


def test_predictions_cover_both_models_and_are_probabilities():
    preds = app_data.predictions("TP53:p.Y220C")
    assert [p.experiment for p in preds] == ["baseline", "baseline_ddg"]
    assert all(0 <= p.out_of_fold <= 1 and 0 <= p.final_model <= 1 for p in preds)


def _run(monkeypatch, **inputs):
    from streamlit.testing.v1 import AppTest

    monkeypatch.setattr(app_data, "fetch_structure", lambda url: "data_empty\n")  # no network in tests
    app = AppTest.from_file(APP, default_timeout=60)
    app.run()
    for key, value in inputs.items():
        getattr(app, key)[0].set_value(value)
    app.run()
    return app


def test_app_renders_pathogenic_and_benign_examples(monkeypatch):
    app = _run(monkeypatch)
    assert not app.exception
    app.sidebar.selectbox[1].set_value("TP53:p.Y220C").run()
    assert "TP53 p.Y220C" in app.subheader[0].value
    assert any("Pathogenic" == m.value for m in app.metric)
    app.sidebar.selectbox[1].set_value("TP53:p.P8L").run()
    assert not app.exception


def test_app_handles_unknown_variant(monkeypatch):
    app = _run(monkeypatch, text_input="R175H")
    assert not app.exception
    assert any("not in the 5,000-variant cohort" in e.value for e in app.error)


def test_typed_change_clears_when_gene_changes(monkeypatch):
    app = _run(monkeypatch, text_input="Y220C")
    assert "TP53 p.Y220C" in app.subheader[0].value
    app.sidebar.selectbox[0].set_value("SOX9").run()
    assert not app.error, "a TP53 change must not be looked up in SOX9"
    assert app.subheader[0].value.startswith("SOX9 p.")
    assert app.text_input[0].value == ""
