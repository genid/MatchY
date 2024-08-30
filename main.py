import streamlit as st
from models import Pedigree, MarkerSet, Simulation
from simulation import run_simulation

st.set_page_config(
        page_title="fraternitY",
        page_icon="🧬",
        layout="wide"
    )

marker_set = MarkerSet()
marker_set.read_marker_set_from_file(r"examples/RM/mutation_rates.csv")

pedigree = Pedigree()
pedigree.read_pedigree_from_file("examples/pedigree_large.tgf")

pedigree.read_known_haplotype_from_file("George", "examples/RM/George.csv", marker_set)
# pedigree.read_known_haplotype_from_file("William", "examples/RM/William.csv", marker_set)
# pedigree.read_known_haplotype_from_file("known2", "manuscript/haplotypes/known_plus1.csv", marker_set)
# pedigree.read_known_haplotype_from_file("known3", "manuscript/haplotypes/known_0.csv", marker_set)

# suspect = st.selectbox("Select a suspect", pedigree.get_known_individuals_names())
suspect = "George"

pedigree.reroot_pedigree(suspect)
selected_node_id = pedigree.visualize_pedigree()

simulation = run_simulation(pedigree, marker_set, suspect, number_of_iterations=100000)
