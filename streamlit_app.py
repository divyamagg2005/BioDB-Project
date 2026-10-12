"""BioDB variant explorer: pathogenicity predictions with AlphaFold structure and FoldX ΔΔG.

Run locally with ``streamlit run streamlit_app.py``.

This module is the frontend only. All data comes through ``biodb.app_data``; nothing
here trains a model or touches the Phase 1-8 pipeline, ``data/``, ``models/`` or
``results/``. The look is a dark, glass-card dashboard with animated ECharts.
"""

from __future__ import annotations

import html
import json
import sys
from pathlib import Path

import numpy as np
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "src"))

from biodb import app_data  # noqa: E402

# streamlit-echarts powers the animated gauges/bars/distribution. It is a new
# dependency (see requirements.txt); if it is missing the page still renders and
# each chart degrades to a small notice instead of crashing.
try:
    from streamlit_echarts import JsCode, st_echarts
except Exception:  # noqa: BLE001
    st_echarts = None
    JsCode = None

st.set_page_config(page_title="BioDB Variant Explorer", layout="wide")

# --- Biological palette: "AlphaFold blues" ---------------------------------
BG = "#0A0F1E"           # deep navy
CARD = "#131A2B"         # navy card
INSET = "#1A2338"        # inset cells (stat grid, expanders)
BENIGN = "#4ADE80"       # benign / stabilising
PATHOGENIC = "#FB7185"   # pathogenic / destabilising
THIS = "#38BDF8"         # the selected variant / mutated residue (cyan)
TEXT = "#E8EDF7"
MUTED = "#8893A8"
TRACK = "#2A3350"        # inactive gauge/axis track


def inject_css() -> None:
    """Dark glassmorphism theme, hidden Streamlit chrome, and a fade-in for cards."""
    st.markdown(
        f"""
        <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Space+Grotesk:wght@500;600;700&display=swap');

        html, body, .stApp, [class*="css"] {{ font-family: 'Inter', system-ui, sans-serif; }}
        /* Hide the menu/deploy/status clutter and footer, but KEEP the toolbar and
           header — the sidebar-expand button lives inside the toolbar, and hiding the
           toolbar (display:none) would remove even its fixed-positioned descendants. */
        #MainMenu, footer, [data-testid="stDecoration"], [data-testid="stToolbarActions"],
        [data-testid="stStatusWidget"], [data-testid="stAppDeployButton"],
        [data-testid="stMainMenu"] {{ display: none !important; }}
        header[data-testid="stHeader"] {{ background: transparent; }}
        /* Button that re-opens a COLLAPSED sidebar. Streamlit 1.6x renders it in the
           header as data-testid="stExpandSidebarButton" but with faint (fadedText60)
           colour, so it is invisible on our dark theme. Make it a clear, fixed button
           on BOTH desktop and mobile so a collapsed sidebar can always be reopened. */
        /* The expand control IS the button element (the test-id sits on the button
           itself, with no child <button>), so size and paint it directly. */
        [data-testid="stExpandSidebarButton"] {{
            display: inline-flex !important; align-items: center; justify-content: center;
            visibility: visible !important; opacity: 1 !important;
            position: fixed !important; z-index: 1000000 !important;
            top: calc(.6rem + env(safe-area-inset-top)); left: .7rem;
            width: 42px; height: 42px; border-radius: 12px;
            background: {CARD} !important; border: 1px solid rgba(255,255,255,.18) !important;
            box-shadow: 0 6px 18px rgba(0,0,0,.5); color: {TEXT} !important; }}
        [data-testid="stSidebarCollapseButton"] button {{
            background: {CARD} !important; border: 1px solid rgba(255,255,255,.18) !important;
            color: {TEXT} !important; border-radius: 12px; }}
        [data-testid="stExpandSidebarButton"] *,
        [data-testid="stSidebarCollapseButton"] * {{ color: {TEXT} !important; fill: {TEXT} !important; }}
        .stApp {{
            background:
                radial-gradient(1100px 520px at 12% -8%, #16203A 0%, rgba(22,32,58,0) 55%),
                radial-gradient(900px 500px at 100% 0%, #10243E 0%, rgba(16,36,62,0) 50%),
                {BG};
        }}
        .block-container {{ padding-top: 2.2rem; padding-bottom: 3rem; }}

        /* Hero */
        .biodb-eyebrow {{ color: {MUTED}; font-size: .8rem; letter-spacing: .22em;
            text-transform: uppercase; margin: 0 0 .35rem; }}
        .biodb-hero-title {{ font-family: 'Space Grotesk', sans-serif; font-weight: 700;
            font-size: clamp(2rem, 6vw, 3.1rem); line-height: 1.02; letter-spacing: -1.5px;
            color: #fff; margin: 0; }}
        .biodb-hero-title .accent {{ color: {BENIGN}; }}

        /* The hero already shows "GENE p.Change"; the real st.subheader is kept for
           programmatic tests but hidden so it isn't shown twice. */
        [data-testid="stHeading"] {{ display: none !important; }}

        /* Links */
        .stApp a {{ color: {BENIGN} !important; text-decoration: none; }}
        .stApp a:hover {{ text-decoration: underline; }}

        /* Metric cards (KPI chips + prediction cards) */
        [data-testid="stMetric"] {{
            background: {CARD}; border: 1px solid rgba(255,255,255,.06);
            border-radius: 16px; padding: 14px 18px;
            box-shadow: 0 10px 28px rgba(0,0,0,.38); }}
        [data-testid="stMetricLabel"] p {{ color: {MUTED} !important;
            font-size: .78rem !important; letter-spacing: .04em; text-transform: uppercase; }}
        [data-testid="stMetricValue"] {{ font-family: 'Space Grotesk', sans-serif;
            font-weight: 700; color: {TEXT}; }}

        /* Glass cards from st.container(border=True) */
        [data-testid="stVerticalBlockBorderWrapper"] {{
            background: {CARD}; border: 1px solid rgba(255,255,255,.06) !important;
            border-radius: 22px; padding: 12px 18px 20px; margin-bottom: 16px;
            box-shadow: 0 16px 40px rgba(0,0,0,.45);
            animation: fadeInUp .55s cubic-bezier(.2,.7,.2,1) both;
            transition: border-color .25s ease, transform .25s ease; }}
        [data-testid="stVerticalBlockBorderWrapper"]:hover {{
            border-color: rgba(56,189,248,.40) !important; transform: translateY(-2px); }}
        @keyframes fadeInUp {{ from {{ opacity: 0; transform: translateY(16px); }}
            to {{ opacity: 1; transform: translateY(0); }} }}

        .biodb-card-title {{ font-family: 'Space Grotesk', sans-serif; font-weight: 600;
            font-size: 1.05rem; color: {TEXT}; margin: .2rem 0 .1rem; }}

        /* KPI chips (value font size is set inline so long text never clips) */
        .biodb-kpi {{ background: {CARD}; border: 1px solid rgba(255,255,255,.06);
            border-radius: 16px; padding: 14px 18px; height: 100%; margin-bottom: 22px;
            box-shadow: 0 10px 28px rgba(0,0,0,.38); }}
        .biodb-kpi .k {{ color: {MUTED}; font-size: .78rem; letter-spacing: .04em;
            text-transform: uppercase; }}
        .biodb-kpi .v {{ font-family: 'Space Grotesk', sans-serif; font-weight: 700;
            color: {TEXT}; line-height: 1.12; margin-top: 4px;
            white-space: normal; overflow-wrap: anywhere; }}

        /* Structural-context / VEP stat grid */
        .biodb-stat-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px;
            margin-top: 6px; margin-bottom: 4px; }}
        .biodb-stat {{ background: {INSET}; border: 1px solid rgba(255,255,255,.05);
            border-radius: 12px; padding: 10px 12px; }}
        .biodb-stat .k {{ color: {MUTED}; font-size: .72rem; text-transform: uppercase;
            letter-spacing: .04em; }}
        .biodb-stat .v {{ color: {TEXT}; font-family: 'Space Grotesk', sans-serif;
            font-weight: 600; font-size: 1.05rem; margin-top: 2px; }}

        /* Sidebar */
        section[data-testid="stSidebar"] {{ background: #0D1424;
            border-right: 1px solid rgba(255,255,255,.06); }}
        section[data-testid="stSidebar"] .stSelectbox label,
        section[data-testid="stSidebar"] .stTextInput label {{ color: {MUTED}; }}

        /* Expanders */
        [data-testid="stExpander"] {{ border: 1px solid rgba(255,255,255,.06);
            border-radius: 14px; background: {INSET}; }}

        /* Scrollbars */
        ::-webkit-scrollbar {{ width: 10px; height: 10px; }}
        ::-webkit-scrollbar-thumb {{ background: {TRACK}; border-radius: 8px; }}

        /* "Results loaded" toast: a glass pill pinned top-center that slides in and
           fades itself out. Shown only when the sidebar selection changes.
           Two identical keyframes (A/B) are alternated on every render: Streamlit
           reuses the toast's DOM node across reruns, and a CSS animation only
           (re)starts when the animation-name changes — so flipping A<->B forces the
           entrance to replay for every gene/variant/mutation change. */
        @keyframes biodbToastA {{
            0%   {{ opacity: 0; transform: translate(-50%, -14px); }}
            10%  {{ opacity: 1; transform: translate(-50%, 0); }}
            82%  {{ opacity: 1; transform: translate(-50%, 0); }}
            100% {{ opacity: 0; transform: translate(-50%, -14px); visibility: hidden; }}
        }}
        @keyframes biodbToastB {{
            0%   {{ opacity: 0; transform: translate(-50%, -14px); }}
            10%  {{ opacity: 1; transform: translate(-50%, 0); }}
            82%  {{ opacity: 1; transform: translate(-50%, 0); }}
            100% {{ opacity: 0; transform: translate(-50%, -14px); visibility: hidden; }}
        }}
        .biodb-toast {{
            position: fixed; top: calc(.9rem + env(safe-area-inset-top)); left: 50%;
            z-index: 2000000; display: inline-flex; align-items: center; gap: .55rem;
            background: {CARD}; border: 1px solid rgba(74,222,128,.42);
            color: {TEXT}; padding: 10px 18px; border-radius: 12px;
            font-size: .9rem; font-weight: 500; white-space: nowrap;
            box-shadow: 0 12px 34px rgba(0,0,0,.55);
            animation-duration: 3.4s; animation-timing-function: cubic-bezier(.2,.7,.2,1);
            animation-fill-mode: forwards;
        }}
        .biodb-toast.tA {{ animation-name: biodbToastA; }}
        .biodb-toast.tB {{ animation-name: biodbToastB; }}
        .biodb-toast .dot {{ width: 8px; height: 8px; border-radius: 50%;
            background: {BENIGN}; box-shadow: 0 0 10px {BENIGN}; }}
        .biodb-toast .muted {{ color: {MUTED}; font-weight: 400; }}

        /* Full-screen loading veil shown while a newly selected variant's data and
           AlphaFold structure load, before the results render. */
        @keyframes biodbSpin {{ to {{ transform: rotate(360deg); }} }}
        @keyframes biodbVeilIn {{ from {{ opacity: 0; }} to {{ opacity: 1; }} }}
        .biodb-loading {{
            position: fixed; inset: 0; z-index: 1500000;
            display: flex; flex-direction: column; align-items: center; justify-content: center;
            gap: 14px; background: rgba(10,15,30,.72); backdrop-filter: blur(3px);
            animation: biodbVeilIn .15s ease both; }}
        .biodb-loading .ring {{
            width: 46px; height: 46px; border-radius: 50%;
            border: 4px solid {TRACK}; border-top-color: {THIS};
            animation: biodbSpin .8s linear infinite; }}
        .biodb-loading .msg {{ color: {MUTED}; font-size: .85rem; letter-spacing: .04em; }}

        /* Phones: stack columns, breathing room clear of the notch/URL bar, and a
           always-visible sidebar hamburger. */
        @media (max-width: 640px) {{
            .block-container {{
                padding: calc(3.8rem + env(safe-area-inset-top)) 1rem
                         calc(3rem + env(safe-area-inset-bottom)); }}
            [data-testid="stHorizontalBlock"] {{ flex-direction: column; gap: .9rem; }}
            [data-testid="stColumn"] {{ width: 100% !important; flex: 1 1 100% !important; }}
            [data-testid="stVerticalBlockBorderWrapper"] {{ margin-bottom: 12px; }}
            .biodb-kpi {{ margin-bottom: 10px; }}
            .biodb-stat-grid {{ grid-template-columns: 1fr; }}
            .biodb-hero-title {{ font-size: 2.1rem; }}
            /* Slightly larger tap target for the expand-sidebar button on phones. */
            [data-testid="stExpandSidebarButton"] {{ width: 46px; height: 46px; }}
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )


def card_title(text: str) -> None:
    st.markdown(f'<div class="biodb-card-title">{text}</div>', unsafe_allow_html=True)


def kpi(label: str, value: str) -> None:
    """A KPI chip whose value font shrinks for long text (e.g. 'Uncertain significance') so it never clips."""
    value = str(value)
    budget = 11  # characters that comfortably fit at the base size
    size = 30 if len(value) <= budget else max(14, int(30 * budget / len(value)))
    st.markdown(
        f'<div class="biodb-kpi"><div class="k">{html.escape(label)}</div>'
        f'<div class="v" style="font-size:{size}px">{html.escape(value)}</div></div>',
        unsafe_allow_html=True,
    )


def notify(message: str) -> None:
    """Top-center glass toast that slides in and fades itself out via CSS.

    Streamlit reuses this markdown element's DOM node across reruns, so the entrance
    animation would only ever play once. Alternating the tA/tB class (and thus the
    CSS animation-name) on each call forces the browser to restart the animation, so
    the toast re-appears for every sidebar change.
    """
    n = st.session_state.get("_toast_n", 0) + 1
    st.session_state["_toast_n"] = n
    cls = "tA" if n % 2 else "tB"
    st.markdown(
        f'<div class="biodb-toast {cls}"><span class="dot"></span>{message}</div>',
        unsafe_allow_html=True,
    )


def stat_grid(values: dict[str, str]) -> None:
    """Small two-column stat grid used for structural context and VEP scores."""
    cells = "".join(
        f'<div class="biodb-stat"><div class="k">{k}</div><div class="v">{v}</div></div>'
        for k, v in values.items()
    )
    st.markdown(f'<div class="biodb-stat-grid">{cells}</div>', unsafe_allow_html=True)


def structure_viewer(cif: str, position: int, height: int = 420) -> None:
    """3Dmol.js cartoon coloured by AlphaFold pLDDT, mutated residue highlighted, dark background."""
    # Escape "</" so the structure text can never close the <script> element.
    cif_json = json.dumps(cif).replace("</", "<\\/")
    # Zero the iframe body margin and hide overflow; the default 8px margin made the frame scroll.
    html = f"""
    <style>html, body {{ margin: 0; padding: 0; overflow: hidden; background: {BG}; }}</style>
    <div id="viewer" style="width:100%;height:{height}px;position:relative;"></div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "{BG}"}});
      viewer.addModel({cif_json}, "cif");
      const plddt = (atom) => atom.b > 90 ? "#0053D6" : atom.b > 70 ? "#65CBF3" : atom.b > 50 ? "#FFDB13" : "#FF7D45";
      viewer.setStyle({{}}, {{cartoon: {{colorfunc: plddt}}}});
      viewer.addStyle({{resi: {position}}}, {{stick: {{colorscheme: "cyanCarbon", radius: 0.3}}}});
      viewer.addStyle({{resi: {position}, atom: "CA"}}, {{sphere: {{color: "{THIS}", radius: 1.2}}}});
      viewer.zoomTo({{within: {{distance: 20, sel: {{resi: {position}}}}}}});
      viewer.render();
    </script>"""
    if hasattr(st, "iframe"):  # Streamlit >= 1.5x; components.html is deprecated there
        st.iframe(html, height=height + 10)
    else:
        import streamlit.components.v1 as components

        components.html(html, height=height + 10)


# ---------------------------------------------------------------------------
# Animated ECharts helpers. Each returns an options dict; render_chart draws it.
# ---------------------------------------------------------------------------
def _stagger(step_ms: int) -> object:
    """Per-element entrance delay (JsCode) so bars grow in one after another."""
    if JsCode is None:
        return 0
    return JsCode(f"function (idx) {{ return idx * {step_ms}; }}").js_code


def render_chart(options: dict, key: str, height: str = "300px") -> None:
    if st_echarts is None:
        st.info("Install `streamlit-echarts` to see this chart.")
        return
    # Transparent canvas so the chart blends with the dark card in any environment.
    options.setdefault("backgroundColor", "transparent")
    try:
        st_echarts(options=options, height=height, key=key)
    except Exception as error:  # noqa: BLE001 - a chart hiccup must not break the page
        st.warning(f"Chart unavailable ({error}).")


def gauge_option(prob: float, label: str) -> dict:
    """Animated gauge with a needle: pathogenic probability as a 0-100% dial."""
    pct = round(prob * 100, 1)
    colour = PATHOGENIC if prob >= 0.5 else BENIGN
    return {
        "animationDuration": 1600,
        "animationEasing": "cubicOut",
        "series": [{
            "type": "gauge",
            "startAngle": 210,
            "endAngle": -30,
            "min": 0,
            "max": 100,
            "radius": "92%",
            "progress": {"show": True, "width": 13, "roundCap": True, "itemStyle": {"color": colour}},
            "axisLine": {"roundCap": True, "lineStyle": {"width": 13, "color": [[1, TRACK]]}},
            "axisTick": {"show": False},
            "splitLine": {"show": False},
            "axisLabel": {"show": False},
            "pointer": {"length": "52%", "width": 5, "itemStyle": {"color": colour}},
            "anchor": {"show": True, "showAbove": True, "size": 11,
                       "itemStyle": {"color": CARD, "borderColor": colour, "borderWidth": 3}},
            "title": {"show": False},
            # The number sits in the dial's bottom opening so the needle/hub never overlaps it.
            "detail": {"valueAnimation": True, "formatter": "{value}%", "color": TEXT,
                       "fontSize": 26, "fontWeight": "bolder", "offsetCenter": [0, "40%"]},
            "data": [{"value": pct, "name": label}],
        }],
    }


def gene_bars_option(gene_rows, sel_position: int) -> dict:
    """ΔΔG of every variant in this gene by position; selected one highlighted magenta."""
    rows = gene_rows.sort_values("protein_position")
    categories, data = [], []
    for pos, lab, ddg, wt, mut in zip(
        rows["protein_position"], rows["label"], rows["ddg"], rows["wt_aa"], rows["mut_aa"]
    ):
        selected = int(pos) == sel_position
        colour = THIS if selected else (PATHOGENIC if lab == 1 else BENIGN)
        categories.append(f"{wt}{int(pos)}{mut}")
        data.append({
            "value": round(float(ddg), 2),
            "itemStyle": {"color": colour, "borderRadius": [5, 5, 0, 0],
                          "shadowBlur": 10 if selected else 0, "shadowColor": THIS},
        })
    many = len(data) > 22
    return {
        "animationDuration": 1200,
        "animationDelay": _stagger(35),
        "grid": {"left": 44, "right": 16, "top": 20, "bottom": 70 if many else 54},
        "tooltip": {"trigger": "axis", "backgroundColor": CARD,
                    "borderColor": "rgba(255,255,255,.1)", "textStyle": {"color": TEXT}},
        "xAxis": {
            "type": "category", "data": categories,
            "axisLabel": {"color": MUTED, "fontSize": 10, "rotate": 60, "interval": 0 if not many else "auto"},
            "axisLine": {"lineStyle": {"color": TRACK}},
            "axisTick": {"show": False},
        },
        "yAxis": {
            "type": "value", "name": "ΔΔG",
            "nameTextStyle": {"color": MUTED}, "axisLabel": {"color": MUTED},
            "splitLine": {"lineStyle": {"color": "rgba(255,255,255,.05)"}},
        },
        "dataZoom": ([{"type": "inside"}, {"type": "slider", "height": 14, "bottom": 10,
                      "borderColor": TRACK, "fillerColor": "rgba(56,189,248,.15)",
                      "handleStyle": {"color": THIS}, "textStyle": {"color": MUTED}}]
                     if many else []),
        "series": [{"type": "bar", "data": data, "barMaxWidth": 26,
                    "markLine": {"silent": True, "symbol": "none",
                                 "lineStyle": {"color": "rgba(255,255,255,.25)", "type": "dashed"},
                                 "data": [{"yAxis": 0}]}}],
    }


def distribution_option(cohort, this_ddg: float) -> dict:
    """Smooth ΔΔG density for benign vs pathogenic, with a marker at this variant."""
    lo, hi = float(cohort["ddg"].quantile(0.01)), float(cohort["ddg"].quantile(0.99))
    grid = np.linspace(lo, hi, 140)
    series = []
    from scipy.stats import gaussian_kde

    for value, colour, name in ((0, BENIGN, "Benign"), (1, PATHOGENIC, "Pathogenic")):
        vals = cohort.loc[cohort["label"] == value, "ddg"].dropna().to_numpy()
        density = gaussian_kde(vals)(grid)
        series.append({
            "name": name, "type": "line", "smooth": True, "showSymbol": False,
            "lineStyle": {"color": colour, "width": 2},
            "areaStyle": {"color": colour, "opacity": 0.18},
            "data": [[round(float(x), 3), round(float(y), 4)] for x, y in zip(grid, density)],
        })
    series[0]["markLine"] = {
        "silent": True, "symbol": "none",
        "label": {"formatter": f"This variant ({this_ddg:+.2f})", "color": THIS, "position": "insideEndTop"},
        "lineStyle": {"color": THIS, "width": 2.4},
        "data": [{"xAxis": round(float(this_ddg), 3)}],
    }
    return {
        "animationDuration": 1400,
        "legend": {"data": ["Benign", "Pathogenic"], "textStyle": {"color": MUTED}, "top": 2},
        "grid": {"left": 44, "right": 18, "top": 36, "bottom": 40},
        "tooltip": {"trigger": "axis", "backgroundColor": CARD,
                    "borderColor": "rgba(255,255,255,.1)", "textStyle": {"color": TEXT}},
        "xAxis": {"type": "value", "name": "ΔΔG (kcal/mol)", "nameLocation": "middle", "nameGap": 26,
                  "nameTextStyle": {"color": MUTED}, "axisLabel": {"color": MUTED},
                  "axisLine": {"lineStyle": {"color": TRACK}},
                  "splitLine": {"lineStyle": {"color": "rgba(255,255,255,.04)"}}},
        "yAxis": {"type": "value", "name": "density", "nameTextStyle": {"color": MUTED},
                  "axisLabel": {"color": MUTED},
                  "splitLine": {"lineStyle": {"color": "rgba(255,255,255,.05)"}}},
        "series": series,
    }


def predictor_bars_option(row, cohort) -> dict:
    """This variant vs cohort means for five predictors, each min-max scaled (higher = more damaging)."""
    features = [("AlphaMissense", "alphamissense", False), ("REVEL", "revel", False),
                ("PolyPhen-2", "polyphen_score", False), ("CADD", "cadd_phred", False),
                ("SIFT (1−score)", "sift_score", True)]
    names, var_v, ben_v, path_v = [], [], [], []
    for name, col, invert in features:
        s = cohort[col].dropna()
        lo, hi = s.min(), s.max()

        def scale(v: float) -> float | None:
            if v != v:  # NaN
                return None
            x = 0.5 if hi <= lo else (v - lo) / (hi - lo)
            x = 1 - x if invert else x
            return round(float(x), 3)

        names.append(name)
        var_v.append(scale(row[col]))
        ben_v.append(scale(cohort.loc[cohort["label"] == 0, col].mean()))
        path_v.append(scale(cohort.loc[cohort["label"] == 1, col].mean()))

    def bar(name, values, colour):
        return {"name": name, "type": "bar", "data": values,
                "barMaxWidth": 11, "itemStyle": {"color": colour, "borderRadius": [0, 4, 4, 0]},
                "animationDelay": _stagger(60)}

    return {
        "animationDuration": 1200,
        "legend": {"data": ["This variant", "Pathogenic mean", "Benign mean"],
                   "textStyle": {"color": MUTED}, "top": 2},
        "grid": {"left": 110, "right": 20, "top": 40, "bottom": 20},
        "tooltip": {"trigger": "axis", "backgroundColor": CARD,
                    "borderColor": "rgba(255,255,255,.1)", "textStyle": {"color": TEXT}},
        "xAxis": {"type": "value", "min": 0, "max": 1, "axisLabel": {"color": MUTED},
                  "splitLine": {"lineStyle": {"color": "rgba(255,255,255,.05)"}}},
        "yAxis": {"type": "category", "data": names, "inverse": True,
                  "axisLabel": {"color": TEXT}, "axisLine": {"lineStyle": {"color": TRACK}},
                  "axisTick": {"show": False}},
        "series": [bar("This variant", var_v, THIS),
                   bar("Pathogenic mean", path_v, PATHOGENIC),
                   bar("Benign mean", ben_v, BENIGN)],
    }


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------
inject_css()

with st.sidebar:
    gene = st.selectbox("Gene", app_data.gene_list(), index=app_data.gene_list().index("TP53"))
    variants = app_data.variants_for(gene)
    example = variants[0].split(":p.")[1]
    # Keyed by gene so the box clears when the gene changes; a change typed for
    # one gene would otherwise be looked up in the next gene.
    typed = st.text_input(f"Protein change (e.g. {example})", "", key=f"change_{gene}")
    if typed:
        variant_id = app_data.parse_variant(gene, typed)
        if variant_id is None:
            st.error(f"{gene} p.{typed.strip().removeprefix('p.').upper()} is not in the 5,000-variant cohort. "
                     "Only cohort variants have precomputed VEP and FoldX features. "
                     "Clear the box to pick a variant from the list.")
            st.stop()
    else:
        variant_id = st.selectbox("Variant", variants, format_func=lambda v: v.split(":")[1])
    st.markdown(
        "**Colour key:** pLDDT > 90 dark blue, 70–90 light blue, 50–70 yellow, < 50 orange. "
        "The mutated residue is shown in cyan."
    )

# Full-screen loading veil: a fixed overlay (high z-index) that covers the page
# while this run fetches data and the AlphaFold structure and renders the charts.
# It stays up as the script runs top-to-bottom and is cleared at the very end
# (_veil.empty()), so the finished results are revealed all at once.
_veil = st.empty()
_veil.markdown(
    '<div class="biodb-loading"><div class="ring"></div>'
    '<div class="msg">Loading variant…</div></div>',
    unsafe_allow_html=True,
)

details = app_data.variant_details(variant_id)
row = details["row"]
position = int(row["protein_position"])
preds = app_data.predictions(variant_id)
cohort = app_data.tables()["dataset"]
gene_rows = cohort[cohort["gene"] == gene]

# --- Hero ------------------------------------------------------------------
change = f"{row['wt_aa']}{position}{row['mut_aa']}"

# Any sidebar change (gene, variant, or a typed mutation) lands on a different
# variant_id and reruns the script. When the new results are ready, flash a toast —
# skipping the very first load, which is not a "change".
_prev_variant = st.session_state.get("_last_variant")
if _prev_variant is not None and _prev_variant != variant_id:
    notify(f'Results loaded <span class="muted">· {gene} p.{change}</span>')
st.session_state["_last_variant"] = variant_id

st.markdown(
    f'<p class="biodb-eyebrow">BioDB Variant Explorer</p>'
    f'<h1 class="biodb-hero-title">{gene} <span class="accent">p.{change}</span></h1>',
    unsafe_allow_html=True,
)
# Kept as a real st.subheader so AppTest can assert the variant string renders.
st.subheader(f"{gene} p.{row['wt_aa']}{position}{row['mut_aa']}")
st.write(
    f"UniProt [{details['uniprot']}](https://www.uniprot.org/uniprotkb/{details['uniprot']}) · "
    f"AlphaFold [{details['alphafold_entry']}](https://alphafold.ebi.ac.uk/entry/{details['uniprot']}) · "
    f"ClinVar [VCV{details['variation_id']:09d}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{details['variation_id']}/)"
)
st.caption(
    "Class project: does protein stability (FoldX ΔΔG on AlphaFold structures) improve pathogenicity "
    "prediction for ClinVar missense variants? Research demonstration only, not for clinical use."
)

# --- KPI chips -------------------------------------------------------------
headline = preds[-1]  # Baseline + FoldX ΔΔG, the experiment this project is about
k1, k2, k3, k4 = st.columns(4)
with k1:
    kpi("Pathogenic prob (OOF)", f"{headline.out_of_fold:.0%}")
with k2:
    kpi("ΔΔG (kcal/mol)", f"{float(row['ddg']):+.2f}")
with k3:
    kpi("pLDDT", f"{row['plddt']:.0f}")
with k4:
    kpi("ClinVar", details["clinvar_significance"])

# --- Row 1: predictions (gauges) | structure ------------------------------
left, right = st.columns([1, 1])
with left:
    with st.container(border=True):
        card_title("Pathogenicity prediction")
        st.caption("Out-of-fold: each probability is from an XGBoost model that never saw this gene.")
        gcols = st.columns(2)
        for col, prediction in zip(gcols, preds):
            with col:
                # st.metric value is the predicted class — AppTest asserts a "Pathogenic" metric exists.
                st.metric(prediction.label, prediction.predicted_class,
                          help="Threshold 0.5 on the out-of-fold probability")
                render_chart(gauge_option(prediction.out_of_fold, "pathogenic"),
                             key=f"gauge_{prediction.experiment}", height="200px")
                st.caption("pathogenic probability")
        with st.expander("Final model (trained on all 5,000 variants, including this one)"):
            for prediction in preds:
                st.write(f"{prediction.label}: pathogenic probability {prediction.final_model:.3f}")
            st.caption("This variant was in the training data, so these numbers are optimistic.")
        with st.expander("ClinVar classification (the training label)"):
            st.write(f"{details['clinvar_significance']} ({details['review_status']})")
            st.caption(details["hgvs"])

with right:
    with st.container(border=True):
        card_title("AlphaFold predicted structure")
        try:
            structure_viewer(app_data.fetch_structure(details["cif_url"]), position)
        except Exception as error:  # noqa: BLE001 - network problems should not break the page
            st.warning(f"Could not load the AlphaFold structure ({error}).")

# --- Row 2: FoldX stability + context | VEP features ----------------------
low, vep = st.columns([1, 1])
with low:
    with st.container(border=True):
        card_title("FoldX: estimated stability change")
        st.metric("ΔΔG (kcal/mol)", f"{float(row['ddg']):+.2f}",
                  help="FoldX 5.1 BuildModel on the repaired AlphaFold model, mean of 3 runs. Positive = destabilising.")
        if row["plddt"] < 70:
            st.warning(f"pLDDT {row['plddt']:.0f} < 70: the structure is low-confidence here, so ΔΔG is unreliable.")
        card_title("Structural context")
        stat_grid({
            "pLDDT": f"{row['plddt']:.1f}",
            "Rel. solvent acc.": f"{row['rsa']:.2f}",
            "Neighbours (10 Å)": str(int(row["neighbor_count"])),
            "Secondary structure": "helix" if row["ss_helix"] else "strand" if row["ss_strand"] else "coil",
        })

with vep:
    with st.container(border=True):
        card_title("Conventional features (VEP)")

        def fmt(value, digits: int = 3) -> str:
            return "missing" if value != value else f"{value:.{digits}g}"

        stat_grid({
            "AlphaMissense": fmt(row["alphamissense"]),
            "REVEL": fmt(row["revel"]),
            "CADD (phred)": fmt(row["cadd_phred"]),
            "SIFT": fmt(row["sift_score"]),
            "PolyPhen-2": fmt(row["polyphen_score"]),
            "GERP": fmt(row["gerp_conservation"]),
            "gnomAD AF": "absent" if row["gnomad_absent"] else fmt(row["gnomad_af_max"]),
            "BLOSUM62": str(int(row["blosum62"])),
        })

# --- Row 3: stability across the gene (full width) ------------------------
with st.container(border=True):
    card_title(f"Stability change across {gene}")
    st.caption("FoldX ΔΔG for every cohort variant in this gene. Magenta is the selected variant.")
    render_chart(gene_bars_option(gene_rows, position), key="gene_bars", height="320px")

# --- Row 4: distribution | predictor comparison ---------------------------
dist, pred = st.columns([1, 1])
with dist:
    with st.container(border=True):
        card_title("ΔΔG: benign vs pathogenic")
        st.caption("Density across all 5,000 variants; the line marks this variant.")
        render_chart(distribution_option(cohort, float(row["ddg"])), key="ddg_dist", height="300px")
with pred:
    with st.container(border=True):
        card_title("Predictor scores vs cohort")
        st.caption("Min-max scaled so higher = more damaging.")
        render_chart(predictor_bars_option(row, cohort), key="predictor_bars", height="300px")

# Everything is rendered: drop the loading veil to reveal the results.
_veil.empty()
