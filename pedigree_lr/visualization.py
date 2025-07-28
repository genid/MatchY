import uuid
from decimal import Decimal
from pathlib import Path
from typing import Union
import streamlit as st
from streamlit_agraph import Config, Edge, Node, agraph
from configparser import ConfigParser
from pedigree_lr.models import Individual, Pedigree, SimulationResult
import matplotlib.pyplot as plt
import glob
import networkx as nx


def _get_node_color(individual: Individual, global_config: ConfigParser) -> str:
    def clean(color: str) -> str:
        return color.strip('"').strip("'")

    if individual.exclude and individual.haplotype_class == "known":
        return clean(global_config["graph"]["NODE_COLOR_EXCLUDED_KNOWN"])
    elif individual.exclude and individual.haplotype_class == "unknown":
        return clean(global_config["graph"]["NODE_COLOR_EXCLUDED_UNKNOWN"])
    elif individual.haplotype_class == "known" and not individual.exclude:
        return clean(global_config["graph"]["NODE_COLOR_KNOWN_HAPLOTYPE"])
    elif individual.haplotype_class == "unknown" and not individual.exclude:
        return clean(global_config["graph"]["NODE_COLOR_UNKNOWN_HAPLOTYPE"])
    elif individual.haplotype_class == "suspect":
        return clean(global_config["graph"]["NODE_COLOR_SUSPECT"])
    raise ValueError(f"Unknown haplotype class {individual.haplotype_class}")



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


def make_plot(files: list[str],
              average: Decimal,
              results_path: Path,
              label: str,
              title: str,
              out_file_name: str) -> None:
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
            probabilities = [float(line) for line in lines[100:]]
            x_values = [i / 10 for i in range(len(probabilities))]
            plt.plot(x_values, probabilities, label=label)
    plt.axhline(y=float(average), color='r', linestyle='--', label='Average')
    plt.xlabel('iterations (x1000)')
    plt.ylabel('probability')
    plt.title(title)
    plt.savefig(f'{results_path}/{out_file_name}.png')
    plt.clf()


def plot_probabilities(
        simulation_result: SimulationResult,
        results_path: Path,
) -> None:
    """
    Plot the probabilities for each l in l_list.
    """
    average_inside_pedigree_probability_files = glob.glob(
        f'{results_path}/average_pedigree_probabilities_m_*_outside_False_*.txt')
    make_plot(
        files=average_inside_pedigree_probability_files,
        average=simulation_result.average_pedigree_probability,
        results_path=results_path,
        label="average_pedigree_probabilities",
        title='Average inside pedigree probabilities',
        out_file_name='average_inside_pedigree_probabilities'
    )

    average_outside_pedigree_probability_files = glob.glob(
        f'{results_path}/average_pedigree_probabilities_m_*_outside_True_*.txt')
    make_plot(
        files=average_outside_pedigree_probability_files,
        average=simulation_result.extended_average_pedigree_probability,
        results_path=results_path,
        label="average_outside_pedigree_probabilities",
        title='Average outside pedigree probabilities',
        out_file_name='average_outside_pedigree_probabilities'
    )

    inside_match_probability_files = glob.glob(
        f'{results_path}/match_probabilities_model_*_outside_False_*.txt')
    make_plot(
        files=inside_match_probability_files,
        average=simulation_result.inside_match_probability[1],
        results_path=results_path,
        label="match_probabilities",
        title='Match probabilities inside pedigree',
        out_file_name='inside_match_probabilities'
    )

    outside_match_probability_files = glob.glob(
        f'{results_path}/match_probabilities_model_*_outside_True_*.txt')
    make_plot(
        files=outside_match_probability_files,
        average=simulation_result.outside_match_probability,
        results_path=results_path,
        label="match_probabilities",
        title='Match probabilities outside pedigree',
        out_file_name='outside_match_probabilities'
    )


def save_pedigree_to_png(pedigree: Pedigree,
                         global_config: ConfigParser,
                         results_path: Path,
                         pedigree_name: str) -> None:
    G = nx.DiGraph()

    # Add nodes with attributes (color, label)
    for individual in pedigree.individuals:
        G.add_node(
            individual.id,
            label=individual.name,
            color=_get_node_color(individual, global_config)
        )

    # Add edges
    for rel in pedigree.relationships:
        G.add_edge(rel.parent_id, rel.child_id)

    # Extract colors and labels
    node_colors = [G.nodes[n]["color"] for n in G.nodes()]
    node_labels = {n: G.nodes[n]["label"] for n in G.nodes()}

    # Layout (spring layout is okay, but hierarchy needs workarounds)
    pos = nx.nx_agraph.graphviz_layout(G, prog='dot')  # Needs pygraphviz installed
    # pos = nx.spring_layout(G, seed=42)  # Use spring layout for better visualization

    # Draw the graph
    plt.figure(figsize=(10, 6))
    nx.draw(G, pos, with_labels=False, arrows=True, node_color=node_colors, node_size=1500)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10, font_color="black")

    out_path = Path(f'{results_path}/{pedigree_name}.png')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, format='png')
    plt.close()