import uuid

import streamlit as st
from streamlit_agraph import Config, Edge, Node, agraph

from pedigree_lr.models import Individual, Pedigree, SimulationResult

# Constants for color codes
NODE_COLOR_KNOWN_HAPLOTYPE = "#b2d3c2"
NODE_COLOR_UNKNOWN_HAPLOTYPE = "#eeeeee"
NODE_COLOR_SUSPECT = "#ff0000"

EDGE_COLOR = "#aaaaaa"


def _get_node_color(individual: Individual) -> str:
    if individual.haplotype_class == "known":
        return NODE_COLOR_KNOWN_HAPLOTYPE
    elif individual.haplotype_class == "unknown":
        return NODE_COLOR_UNKNOWN_HAPLOTYPE
    elif individual.haplotype_class == "suspect":
        return NODE_COLOR_SUSPECT
    raise ValueError(f"Unknown haplotype class {individual.haplotype_class}")


def st_print_pedigree(pedigree: Pedigree) -> None:
    for individual in pedigree.individuals:
        for allele in individual.haplotype.alleles.values():
            allele_str = f"{allele.value}.{allele.intermediate_value}" if allele.intermediate_value is not None else str(allele.value)
            parent_str = f"{allele.parent_value}.{allele.parent_intermediate_value}" if allele.parent_intermediate_value is not None else str(allele.parent_value)

            st.write(
                f"{individual.name}, {individual.haplotype_class}, {allele.marker.name}, "
                f"{allele_str}, {parent_str}, {allele.mutation_value}, "
                f"{allele.mutation_probability}\n"
            )


def st_visualize_pedigree(pedigree: Pedigree) -> int:
    config = Config(
        width=700,
        height=700,
        directed=True,
        hierarchical=True,
        direction="UD",
        sortMethod="directed",
        physics=False,
        nodeSpacing=150,
        key=str(uuid.uuid4()),
    )

    nodes = [
        Node(id=individual.id, label=individual.name, color=_get_node_color(individual))
        for individual in pedigree.individuals
    ]

    edges = [
        Edge(
            source=relationship.parent_id,
            target=relationship.child_id,
            color=EDGE_COLOR,
        )
        for relationship in pedigree.relationships
    ]

    selected_node_id = agraph(nodes=nodes, edges=edges, config=config)

    return selected_node_id
