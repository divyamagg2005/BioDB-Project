"""BioDB variant explorer: pathogenicity predictions with AlphaFold structure and FoldX ΔΔG.

Run locally with ``streamlit run streamlit_app.py``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "src"))

from biodb import app_data  # noqa: E402

st.set_page_config(page_title="BioDB Variant Explorer", layout="wide")


def structure_viewer(cif: str, position: int, height: int = 480) -> None:
    """3Dmol.js cartoon coloured by AlphaFold pLDDT, with the mutated residue highlighted."""
    # Escape "</" so the structure text can never close the <script> element.
    cif_json = json.dumps(cif).replace("</", "<\\/")
    # Zero the iframe body margin and hide overflow; the default 8px margin made the frame scroll.
    html = f"""
    <style>html, body {{ margin: 0; padding: 0; overflow: hidden; }}</style>
    <div id="viewer" style="width:100%;height:{height}px;position:relative;"></div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
      viewer.addModel({cif_json}, "cif");
      const plddt = (atom) => atom.b > 90 ? "#0053D6" : atom.b > 70 ? "#65CBF3" : atom.b > 50 ? "#FFDB13" : "#FF7D45";
      viewer.setStyle({{}}, {{cartoon: {{colorfunc: plddt}}}});
      viewer.addStyle({{resi: {position}}}, {{stick: {{colorscheme: "magentaCarbon", radius: 0.3}}}});
      viewer.addStyle({{resi: {position}, atom: "CA"}}, {{sphere: {{color: "magenta", radius: 1.2}}}});
      viewer.zoomTo({{within: {{distance: 20, sel: {{resi: {position}}}}}}});
      viewer.render();
    </script>"""
    if hasattr(st, "iframe"):  # Streamlit >= 1.5x; components.html is deprecated there
        st.iframe(html, height=height + 10)
    else:
        import streamlit.components.v1 as components

        components.html(html, height=height + 10)


def show_table(values: dict) -> None:
    """Two-column table with every value as text, so Arrow serialisation never sees mixed types."""
    st.table(pd.DataFrame({"Value": [str(v) for v in values.values()]}, index=list(values)))


st.title("BioDB Variant Explorer")
st.caption(
    "Class project: does protein stability (FoldX ΔΔG on AlphaFold structures) improve pathogenicity "
    "prediction for ClinVar missense variants? Research demonstration only, not for clinical use."
)

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
        "The mutated residue is shown in magenta."
    )

details = app_data.variant_details(variant_id)
row = details["row"]
position = int(row["protein_position"])

st.subheader(f"{gene} p.{row['wt_aa']}{position}{row['mut_aa']}")
st.write(
    f"UniProt [{details['uniprot']}](https://www.uniprot.org/uniprotkb/{details['uniprot']}) · "
    f"AlphaFold [{details['alphafold_entry']}](https://alphafold.ebi.ac.uk/entry/{details['uniprot']}) · "
    f"ClinVar [VCV{details['variation_id']:09d}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{details['variation_id']}/)"
)

st.markdown("### Machine learning: pathogenicity prediction")
columns = st.columns(2)
for column, prediction in zip(columns, app_data.predictions(variant_id)):
    with column:
        st.metric(prediction.label, prediction.predicted_class, help="Threshold 0.5 on the out-of-fold probability")
        st.write(f"Pathogenic probability **{prediction.out_of_fold:.3f}** · benign probability **{1 - prediction.out_of_fold:.3f}**")
st.caption(
    "Probabilities are out-of-fold: they come from XGBoost models trained in gene-grouped cross-validation "
    f"without any {gene} variants, so they show how the model performs on an unseen gene."
)
with st.expander("Final model (trained on all 5,000 variants, including this one)"):
    for prediction in app_data.predictions(variant_id):
        st.write(f"{prediction.label}: pathogenic probability {prediction.final_model:.3f}")
    st.caption("This variant was in the training data, so these numbers are optimistic. Use the out-of-fold values above.")
with st.expander("ClinVar classification (the training label)"):
    st.write(f"{details['clinvar_significance']} ({details['review_status']})")
    st.caption(details["hgvs"])

left, right = st.columns([3, 2])
with left:
    st.markdown("### AlphaFold: predicted structure")
    try:
        structure_viewer(app_data.fetch_structure(details["cif_url"]), position)
    except Exception as error:  # noqa: BLE001 - network problems should not break the page
        st.warning(f"Could not load the AlphaFold structure ({error}).")
with right:
    st.markdown("### FoldX: estimated stability change")
    ddg = float(row["ddg"])
    st.metric("ΔΔG (kcal/mol)", f"{ddg:+.2f}",
              help="FoldX 5.1 BuildModel on the repaired AlphaFold model, mean of 3 runs. Positive = destabilising.")
    if row["plddt"] < 70:
        st.warning(f"pLDDT {row['plddt']:.0f} < 70: the structure is low-confidence here, so ΔΔG is unreliable.")
    st.markdown("### Structural context")
    show_table({
        "pLDDT": f"{row['plddt']:.1f}",
        "Relative solvent accessibility": f"{row['rsa']:.2f}",
        "Neighbour count (10 Å)": int(row["neighbor_count"]),
        "Secondary structure": "helix" if row["ss_helix"] else "strand" if row["ss_strand"] else "coil",
    })

st.markdown("### Conventional features (VEP)")


def fmt(value, digits: int = 3) -> str:
    return "missing" if value != value else f"{value:.{digits}g}"


show_table({
    "AlphaMissense": fmt(row["alphamissense"]),
    "REVEL": fmt(row["revel"]),
    "CADD (phred)": fmt(row["cadd_phred"]),
    "SIFT": fmt(row["sift_score"]),
    "PolyPhen-2": fmt(row["polyphen_score"]),
    "GERP": fmt(row["gerp_conservation"]),
    "gnomAD AF": "absent from gnomAD" if row["gnomad_absent"] else fmt(row["gnomad_af_max"]),
    "BLOSUM62": int(row["blosum62"]),
})

# ---------------------------------------------------------------------------
# Graphs: where this induced mutation sits relative to the rest of the cohort.
# Read-only; everything comes from the data app_data already loaded.
# ---------------------------------------------------------------------------
BENIGN, PATHOGENIC, THIS = "#9E9E9E", "#C62828", "#D81B9A"
sns.set_theme(style="whitegrid", context="notebook")

cohort = app_data.tables()["dataset"]
gene_rows = cohort[cohort["gene"] == gene]


def ddg_landscape() -> plt.Figure:
    """FoldX ΔΔG for every variant in this gene by position; this mutation starred."""
    fig, ax = plt.subplots(figsize=(8, 3.2))
    for mask, colour, label in ((gene_rows["label"] == 0, BENIGN, "Benign"),
                                (gene_rows["label"] == 1, PATHOGENIC, "Pathogenic")):
        pts = gene_rows[mask]
        ax.scatter(pts["protein_position"], pts["ddg"], c=colour, s=28, alpha=0.7,
                   edgecolor="white", linewidth=0.4, label=f"{label} ({len(pts)})")
    ax.scatter([position], [row["ddg"]], c=THIS, s=320, marker="*", zorder=5,
               edgecolor="black", linewidth=0.6, label="This variant")
    ax.axhline(0, color="black", lw=0.6, ls="--")
    ax.set_xlabel("Protein position")
    ax.set_ylabel("ΔΔG (kcal/mol)")
    ax.set_title(f"Stability change across {gene}")
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    return fig


def ddg_distribution() -> plt.Figure:
    """ΔΔG distribution for benign vs pathogenic across the whole cohort, with this variant marked."""
    fig, ax = plt.subplots(figsize=(8, 3.2))
    clip = (cohort["ddg"].quantile(0.01), cohort["ddg"].quantile(0.99))
    for value, colour, label in ((0, BENIGN, "Benign"), (1, PATHOGENIC, "Pathogenic")):
        sns.kdeplot(cohort.loc[cohort["label"] == value, "ddg"], ax=ax, fill=True,
                    alpha=0.4, color=colour, clip=clip, label=label)
    ax.axvline(float(row["ddg"]), color=THIS, lw=2.2, label=f"This variant ({row['ddg']:+.2f})")
    ax.set_xlabel("ΔΔG (kcal/mol)")
    ax.set_ylabel("Density")
    ax.set_title("ΔΔG: benign vs pathogenic (5,000 variants)")
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    return fig


def predictor_comparison() -> plt.Figure:
    """This variant's predictor scores vs cohort means, each min-max scaled so higher = more damaging."""
    features = [("AlphaMissense", "alphamissense", False), ("REVEL", "revel", False),
                ("PolyPhen-2", "polyphen_score", False), ("CADD", "cadd_phred", False),
                ("SIFT (1−score)", "sift_score", True)]
    names, var_v, ben_v, path_v = [], [], [], []
    for name, col, invert in features:
        series = cohort[col].dropna()
        lo, hi = series.min(), series.max()

        def scale(v: float) -> float:
            if v != v:  # NaN
                return np.nan
            s = 0.5 if hi <= lo else (v - lo) / (hi - lo)
            return 1 - s if invert else s

        names.append(name)
        var_v.append(scale(row[col]))
        ben_v.append(scale(cohort.loc[cohort["label"] == 0, col].mean()))
        path_v.append(scale(cohort.loc[cohort["label"] == 1, col].mean()))

    y = np.arange(len(names))
    height = 0.26
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.barh(y + height, var_v, height, color=THIS, label="This variant")
    ax.barh(y, path_v, height, color=PATHOGENIC, alpha=0.8, label="Pathogenic mean")
    ax.barh(y - height, ben_v, height, color=BENIGN, alpha=0.8, label="Benign mean")
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Scaled score (higher = more damaging)")
    ax.set_title("Predictor scores vs cohort")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    return fig


st.markdown("---")
st.markdown("### Graphs: this mutation in context")
st.caption(
    "How this variant compares with the rest of the gene and the 5,000-variant cohort. "
    "Magenta marks the variant selected above."
)

figure = ddg_landscape()
st.pyplot(figure)
plt.close(figure)

graph_left, graph_right = st.columns(2)
with graph_left:
    figure = ddg_distribution()
    st.pyplot(figure)
    plt.close(figure)
with graph_right:
    figure = predictor_comparison()
    st.pyplot(figure)
    plt.close(figure)
