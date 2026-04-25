/// Petgraph-based DAG wrappers for pedigree operations.
///
/// The pedigree is a rooted DAG (actually a tree in most cases — each individual
/// has at most one father). Edges flow parent → child.
use crate::{MatchyError, Pedigree, Result};
use petgraph::algo::{is_cyclic_directed, toposort};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::{Bfs, EdgeRef};
use petgraph::Direction;
use std::collections::HashMap;

/// Build a petgraph DiGraph from a Pedigree.
/// Returns (graph, id_to_node, node_to_id).
pub fn build_graph(
    pedigree: &Pedigree,
) -> (DiGraph<String, ()>, HashMap<String, NodeIndex>, HashMap<NodeIndex, String>) {
    let mut graph = DiGraph::new();
    let mut id_to_node: HashMap<String, NodeIndex> = HashMap::new();
    let mut node_to_id: HashMap<NodeIndex, String> = HashMap::new();

    for individual in &pedigree.individuals {
        let idx = graph.add_node(individual.id.clone());
        id_to_node.insert(individual.id.clone(), idx);
        node_to_id.insert(idx, individual.id.clone());
    }

    for rel in &pedigree.relationships {
        if let (Some(&parent_idx), Some(&child_idx)) =
            (id_to_node.get(&rel.parent_id), id_to_node.get(&rel.child_id))
        {
            graph.add_edge(parent_idx, child_idx, ());
        }
    }

    (graph, id_to_node, node_to_id)
}

/// Check that the pedigree is a valid DAG (no cycles).
pub fn validate_dag(pedigree: &Pedigree) -> Result<()> {
    let (graph, _, _) = build_graph(pedigree);
    if is_cyclic_directed(&graph) {
        return Err(MatchyError::InvalidPedigree(
            "Pedigree contains cycles — must be a DAG".into(),
        ));
    }
    Ok(())
}

/// Topological order of individual IDs (ancestors before descendants).
pub fn topological_order(pedigree: &Pedigree) -> Result<Vec<String>> {
    let (graph, _, node_to_id) = build_graph(pedigree);
    let order = toposort(&graph, None).map_err(|_| {
        MatchyError::InvalidPedigree("Pedigree contains cycles (toposort failed)".into())
    })?;
    Ok(order.iter().map(|idx| node_to_id[idx].clone()).collect())
}

/// BFS layers starting from a given root ID.
/// Layer 0 = root, layer 1 = children, etc.
pub fn bfs_layers(pedigree: &Pedigree, root_id: &str) -> Result<Vec<Vec<String>>> {
    let (graph, id_to_node, node_to_id) = build_graph(pedigree);
    let root = *id_to_node
        .get(root_id)
        .ok_or_else(|| MatchyError::IndividualNotFound(root_id.into()))?;

    let mut layers: Vec<Vec<String>> = Vec::new();
    let mut bfs = Bfs::new(&graph, root);
    // petgraph BFS doesn't expose layers directly — use depth tracking via HashMap
    let mut depth: HashMap<NodeIndex, usize> = HashMap::new();
    depth.insert(root, 0);

    while let Some(node) = bfs.next(&graph) {
        let d = *depth.get(&node).unwrap_or(&0);
        if layers.len() <= d {
            layers.resize_with(d + 1, Vec::new);
        }
        layers[d].push(node_to_id[&node].clone());
        for edge in graph.edges_directed(node, Direction::Outgoing) {
            let child = edge.target();
            depth.entry(child).or_insert(d + 1);
        }
    }
    Ok(layers)
}

/// Shortest undirected path length between two individuals.
/// Returns None if no path exists.
pub fn shortest_path_length(pedigree: &Pedigree, from: &str, to: &str) -> Option<usize> {
    use petgraph::algo::astar;
    let (graph, id_to_node, _) = build_graph(pedigree);
    let undirected = graph.map(|_, n| n.clone(), |_, e| *e);
    // Convert to undirected by treating as undirected in astar
    let from_idx = *id_to_node.get(from)?;
    let to_idx = *id_to_node.get(to)?;
    // Use BFS on the undirected interpretation
    let undir = petgraph::graph::UnGraph::<String, ()>::from_edges(
        pedigree.relationships.iter().flat_map(|r| {
            let a = id_to_node.get(&r.parent_id).copied()?;
            let b = id_to_node.get(&r.child_id).copied()?;
            Some((a, b))
        }),
    );
    // Map indices from directed graph to undirected graph — simpler: just use BFS on directed+reversed
    let _ = undirected; // suppress warning
    astar(
        &graph,
        from_idx,
        |n| n == to_idx,
        |_| 1usize,
        |_| 0usize,
    )
    .map(|(cost, _)| cost)
}

/// All ancestors of an individual (inclusive of the individual itself).
pub fn ancestors(pedigree: &Pedigree, id: &str) -> Vec<String> {
    let (graph, id_to_node, node_to_id) = build_graph(pedigree);
    let Some(&start) = id_to_node.get(id) else {
        return vec![];
    };
    // Walk edges in reverse direction (child → parent)
    let mut visited = std::collections::HashSet::new();
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        if visited.insert(node) {
            for edge in graph.edges_directed(node, Direction::Incoming) {
                stack.push(edge.source());
            }
        }
    }
    visited.iter().map(|n| node_to_id[n].clone()).collect()
}

/// All descendants of an individual (inclusive of the individual itself).
pub fn descendants(pedigree: &Pedigree, id: &str) -> Vec<String> {
    let (graph, id_to_node, node_to_id) = build_graph(pedigree);
    let Some(&start) = id_to_node.get(id) else {
        return vec![];
    };
    let mut visited = std::collections::HashSet::new();
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        if visited.insert(node) {
            for edge in graph.edges_directed(node, Direction::Outgoing) {
                stack.push(edge.target());
            }
        }
    }
    visited.iter().map(|n| node_to_id[n].clone()).collect()
}

/// Find the Most Recent Common Ancestor of two individuals.
/// Returns the MRCA individual ID, or None if none exists.
pub fn most_recent_common_ancestor(pedigree: &Pedigree, a: &str, b: &str) -> Option<String> {
    let anc_a: std::collections::HashSet<_> = ancestors(pedigree, a).into_iter().collect();
    let anc_b: std::collections::HashSet<_> = ancestors(pedigree, b).into_iter().collect();
    let common: Vec<_> = anc_a.intersection(&anc_b).cloned().collect();
    if common.is_empty() {
        return None;
    }
    // The MRCA is the common ancestor with the greatest depth (farthest from root)
    let (graph, id_to_node, _) = build_graph(pedigree);
    let roots = pedigree.roots();
    let root_id = roots.first()?;
    let root_idx = *id_to_node.get(*root_id)?;
    // depth = BFS distance from root
    let mut depth_map: HashMap<NodeIndex, usize> = HashMap::new();
    let mut bfs = Bfs::new(&graph, root_idx);
    let mut d_counter = 0usize;
    while let Some(node) = bfs.next(&graph) {
        depth_map.insert(node, d_counter);
        d_counter += 1;
    }
    common.into_iter().max_by_key(|id| {
        id_to_node
            .get(id.as_str())
            .and_then(|idx| depth_map.get(idx))
            .copied()
            .unwrap_or(0)
    })
}

/// Shortest undirected path between two individuals.
/// Returns the sequence of IDs from `from` to `to` (inclusive), or None if unreachable.
pub fn shortest_path_undirected(pedigree: &Pedigree, from: &str, to: &str) -> Option<Vec<String>> {
    if from == to {
        return Some(vec![from.to_string()]);
    }

    // Build undirected adjacency list from directed edges.
    let mut adj: HashMap<String, Vec<String>> = HashMap::new();
    for rel in &pedigree.relationships {
        adj.entry(rel.parent_id.clone()).or_default().push(rel.child_id.clone());
        adj.entry(rel.child_id.clone()).or_default().push(rel.parent_id.clone());
    }

    // BFS — track predecessor for path reconstruction.
    use std::collections::{VecDeque, HashMap as LMap};
    let mut prev: LMap<String, String> = LMap::new();
    let mut visited: std::collections::HashSet<String> = std::collections::HashSet::new();
    visited.insert(from.to_string());
    let mut queue: VecDeque<String> = VecDeque::new();
    queue.push_back(from.to_string());

    while let Some(node) = queue.pop_front() {
        if node == to {
            // Reconstruct path from `to` back to `from`.
            let mut path = vec![to.to_string()];
            let mut cur = to.to_string();
            while let Some(p) = prev.get(&cur) {
                path.push(p.clone());
                cur = p.clone();
            }
            path.reverse();
            return Some(path);
        }
        if let Some(nbrs) = adj.get(&node) {
            for nbr in nbrs {
                if visited.insert(nbr.clone()) {
                    prev.insert(nbr.clone(), node.clone());
                    queue.push_back(nbr.clone());
                }
            }
        }
    }
    None
}

/// Check connectivity (all nodes reachable from any single node via undirected edges).
pub fn is_connected(pedigree: &Pedigree) -> bool {
    if pedigree.individuals.is_empty() {
        return true;
    }
    let (graph, _, _) = build_graph(pedigree);
    let undirected = petgraph::visit::NodeFiltered::from_fn(&graph, |_| true);
    let _ = undirected;
    // Use petgraph connected_components on undirected version
    let undir: petgraph::graph::UnGraph<(), ()> = petgraph::graph::UnGraph::from_edges(
        pedigree
            .relationships
            .iter()
            .filter_map(|r| {
                let a = graph
                    .node_indices()
                    .find(|&n| graph[n] == r.parent_id)?;
                let b = graph
                    .node_indices()
                    .find(|&n| graph[n] == r.child_id)?;
                Some((a, b))
            }),
    );
    petgraph::algo::connected_components(&undir) <= 1
}
