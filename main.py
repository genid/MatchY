import streamlit as st
from models import Pedigree, MarkerSet, Simulation
from simulation import run_simulation

st.set_page_config(
        page_title="fraternitY",
        page_icon="🧬",
        layout="wide"
    )

marker_set = MarkerSet()
marker_set.read_marker_set_from_file("examples/mutation_rates.csv")

pedigree = Pedigree()
pedigree.read_pedigree_from_file("examples/pedigree_example2.tgf")

# pedigree.read_known_haplotype_from_file("Archie", "examples/Archie.csv", marker_set)
# pedigree.read_known_haplotype_from_file("Edward", "examples/Edward.csv", marker_set)
pedigree.read_known_haplotype_from_file("George", "examples/George.csv", marker_set)
pedigree.read_known_haplotype_from_file("William", "examples/William.csv", marker_set)
# pedigree.read_known_haplotype_from_file("test2", "examples/Louis.csv", marker_set)

suspect = st.selectbox("Select a suspect", pedigree.get_known_individuals_names())

pedigree.reroot_pedigree(suspect)
selected_node_id = pedigree.visualize_pedigree()

run_simulation(pedigree, marker_set, suspect, number_of_iterations=100000)

