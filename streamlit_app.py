"""BioDB variant explorer: pathogenicity predictions with AlphaFold structure and FoldX ΔΔG.

Run locally with ``streamlit run streamlit_app.py``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

sys.path.insert(0, str(Path(__file__).parent / "src"))

from biodb import app_data  # noqa: E402

st.set_page_config(page_title="BioDB Variant Explorer", layout="wide")


def structure_viewer(cif: str, position: int, height: int = 480) -> None:
    """3Dmol.js cartoon coloured by AlphaFold pLDDT, with the mutated residue highlighted."""
    html = f"""
    <div id="viewer" style="width:100%;height:{height}px;position:relative;"></div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
      viewer.addModel({json.dumps(cif)}, "cif");
      const plddt = (atom) => atom.b > 90 ? "#0053D6" : atom.b > 70 ? "#65CBF3" : atom.b > 50 ? "#FFDB13" : "#FF7D45";
      viewer.setStyle({{}}, {{cartoon: {{colorfunc: plddt}}}});
      viewer.addStyle({{resi: {position}}}, {{stick: {{colorscheme: "magentaCarbon", radius: 0.3}}}});
      viewer.addStyle({{resi: {position}, atom: "CA"}}, {{sphere: {{color: "magenta", radius: 1.2}}}});
      viewer.zoomTo({{within: {{distance: 20, sel: {{resi: {position}}}}}}});
      viewer.render();
    </script>"""
    components.html(html, height=height + 10)


st.title("BioDB Variant Explorer")
st.caption(
    "Class project: does protein stability (FoldX ΔΔG on AlphaFold structures) improve pathogenicity "
    "prediction for ClinVar missense variants? Research demonstration only, not for clinical use."
)

with st.sidebar:
    gene = st.selectbox("Gene", app_data.gene_list(), index=app_data.gene_list().index("TP53"))
    variants = app_data.variants_for(gene)
    typed = st.text_input("Protein change (e.g. Y220C)", "")
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
    st.table({
        "pLDDT": f"{row['plddt']:.1f}",
        "Relative solvent accessibility": f"{row['rsa']:.2f}",
        "Neighbour count (10 Å)": int(row["neighbor_count"]),
        "Secondary structure": "helix" if row["ss_helix"] else "strand" if row["ss_strand"] else "coil",
    })

st.markdown("### Conventional features (VEP)")


def fmt(value, digits: int = 3) -> str:
    return "missing" if value != value else f"{value:.{digits}g}"


st.table({
    "AlphaMissense": fmt(row["alphamissense"]),
    "REVEL": fmt(row["revel"]),
    "CADD (phred)": fmt(row["cadd_phred"]),
    "SIFT": fmt(row["sift_score"]),
    "PolyPhen-2": fmt(row["polyphen_score"]),
    "GERP": fmt(row["gerp_conservation"]),
    "gnomAD AF": "absent from gnomAD" if row["gnomad_absent"] else fmt(row["gnomad_af_max"]),
    "BLOSUM62": int(row["blosum62"]),
})
