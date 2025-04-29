import uuid
from pathlib import Path

import streamlit as st
from streamlit_agraph import Config, Edge, Node, agraph
from configparser import ConfigParser
from pedigree_lr.models import Individual, Pedigree, SimulationResult
import matplotlib.pyplot as plt
import glob

NODE_COLOR_KNOWN_HAPLOTYPE = "#b2d3c2"
NODE_COLOR_UNKNOWN_HAPLOTYPE = "#eeeeee"
NODE_COLOR_SUSPECT = "#ff0000"
NODE_COLOR_EXCLUDED = "#888888"
EDGE_COLOR = "#aaaaaa"


def _get_node_color(individual: Individual,
                    global_config: ConfigParser) -> str:
    if individual.exclude:
        return global_config["graph"]["NODE_COLOR_EXCLUDED"]
    elif individual.haplotype_class == "known":
        return global_config["graph"]["NODE_COLOR_KNOWN_HAPLOTYPE"]
    elif individual.haplotype_class == "unknown":
        return global_config["graph"]["NODE_COLOR_UNKNOWN_HAPLOTYPE"]
    elif individual.haplotype_class == "suspect":
        return global_config["graph"]["NODE_COLOR_SUSPECT"]
    raise ValueError(f"Unknown haplotype class {individual.haplotype_class}")


def st_print_pedigree(pedigree: Pedigree) -> None:
    for individual in pedigree.individuals:
        for allele in individual.haplotype.alleles.values():
            allele_str = f"{allele.value}.{allele.intermediate_value}" if allele.intermediate_value is not None else str(
                allele.value)
            parent_str = f"{allele.parent_value}.{allele.parent_intermediate_value}" if allele.parent_intermediate_value is not None else str(
                allele.parent_value)

            st.write(
                f"{individual.name}, {individual.haplotype_class}, {allele.marker.name}, "
                f"{allele_str}, {parent_str}, {allele.mutation_value}, "
                f"{allele.mutation_probability}\n"
            )


def st_visualize_pedigree(pedigree: Pedigree,
                          global_config: ConfigParser) -> int:
    graph = Config(
        width=global_config["pedigree_window"]["width"],
        height=global_config["pedigree_window"]["height"],
        directed=True,
        hierarchical=True,
        direction="UD",
        sortMethod="directed",
        physics=False,
        nodeSpacing=150,
        key=str(uuid.uuid4()),
    )

    nodes = [
        Node(id=individual.id, label=individual.name, color=_get_node_color(individual, global_config))
        for individual in pedigree.individuals
    ]

    edges = [
        Edge(
            source=relationship.parent_id,
            target=relationship.child_id,
            color=global_config["graph"]["EDGE_COLOR"],
        )
        for relationship in pedigree.relationships
    ]

    selected_node_id = agraph(nodes=nodes, edges=edges, config=graph)

    return selected_node_id


def plot_probabilities(
        results_path: Path,
        l_list: list[int],
) -> None:
    """
    Plot the probabilities for each l in l_list.
    """
    files = glob.glob(f'{results_path}/average_pedigree_probabilities_*.txt')
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
            probabilities = [float(line) for line in lines[:]]
            plt.plot(probabilities, label="average_pedigree_probabilities")
        if "True" in file:
            plt.title(f'Average pedigree probabilities')
            plt.savefig(f'{results_path}/average_pedigree_probabilities_outside_pedigree.png')
        else:
            plt.title(f'Average outside pedigree probabilities')
            plt.savefig(f'{results_path}/average_pedigree_probabilities.png')
        plt.clf()

    for l in l_list:
        if l == 0:
            continue
        files = glob.glob(f"{results_path}/{l}_pedigree_probabilities_model_*False*.txt")

        for file in files:
            with open(file, "r") as f:
                lines = f.readlines()
                probabilities = [float(line) for line in lines[:]]
                plt.plot(probabilities, label=file)
        plt.title(f"Probabilities for l={l}")
        plt.savefig(f"{results_path}/l_{l}_probabilities.png")
        plt.clf()

    files = glob.glob(f"{results_path}/1_pedigree_probabilities_model_*True*.txt")

    for file in files:
        with open(file, "r") as f:
            lines = f.readlines()
            probabilities = [float(line) for line in lines[:]]
            plt.plot(probabilities, label=file)
    plt.title(f"Outside match probability")
    plt.savefig(f"{results_path}/outside_match_probability.png")
    plt.clf()