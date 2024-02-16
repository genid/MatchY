import streamlit as st
from models import Pedigree, MarkerSet, Simulation
from simulation import run_simulation

st.set_page_config(
        page_title="fraternitY",
        page_icon="🧬",
        layout="wide"
    )

marker_set = MarkerSet()
marker_set.read_marker_set_from_file("mutation_rates.csv")

pedigree = Pedigree()
pedigree.read_pedigree_from_file("pedigree_large.tgf")

pedigree.read_known_haplotype_from_file("Archie", "Archie.csv", marker_set)
pedigree.read_known_haplotype_from_file("Edward", "Edward.csv", marker_set)
pedigree.read_known_haplotype_from_file("George", "George.csv", marker_set)
pedigree.read_known_haplotype_from_file("Louis", "Louis.csv", marker_set)

suspect = st.selectbox("Select a suspect", pedigree.get_known_individuals_names())

pedigree.reroot_pedigree(suspect)
selected_node_id = pedigree.visualize_pedigree()

run_simulation(pedigree, marker_set, suspect, number_of_iterations=1000)

