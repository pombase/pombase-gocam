//! Graph algorithms

use std::collections::HashSet;

use petgraph::{Direction, Graph};
use petgraph::{graph::NodeIndex, visit::EdgeRef};
//use petgraph::algo::is_isomorphic_subgraph;

use anyhow::{Result, anyhow};

/// A predicate function to pass to [subgraph_by].  The arguments are an initial node and the
/// current node while traversing the graph.
pub type SubGraphPred<N> =
    fn(start_node: &N, test_node: &N) -> bool;

/// Return a sub-graph of `graph` that includes the node given by `start_idx` and all connected
/// nodes and edges where `match_pred` returns true.  `match_pred` is passed the start node and
/// the current node while traversing the graph.
pub fn subgraph_by<N: Clone, E: Clone>(graph: &Graph<N, E>,
                                       start_idx: NodeIndex,
                                       match_pred: &SubGraphPred<N>)
    -> Result<Graph<N, E>>
{
    let mut stack = vec![];
    let mut seen_nodes = HashSet::new();

    let mut ret_graph = Graph::<N, E>::new();

    let old_start_node = graph.node_weight(start_idx)
        .ok_or(anyhow!("node not found: {:?}", start_idx))?;

    let new_start_idx = ret_graph.add_node(old_start_node.to_owned());

    stack.push((start_idx, new_start_idx));
    seen_nodes.insert(start_idx);

    loop {
        let Some((old_current_idx, new_current_idx)) = stack.pop()
        else {
            break;
        };

        let outgoing_iter = graph.edges_directed(old_current_idx, Direction::Outgoing)
            .map(|e| (e, Direction::Outgoing));
        let incoming_iter = graph.edges_directed(old_current_idx, Direction::Incoming)
            .map(|e| (e, Direction::Incoming));

        for (old_edge_ref, direction) in outgoing_iter.chain(incoming_iter) {
            let old_other_node_idx =
                if direction == Direction::Outgoing {
                    old_edge_ref.target()
                } else {
                    old_edge_ref.source()
                };
            let other_node = graph.node_weight(old_other_node_idx).unwrap();

            if seen_nodes.contains(&old_other_node_idx) {
                continue;
            } else {
                seen_nodes.insert(old_other_node_idx);
            }

            if match_pred(old_start_node, other_node) {
                let new_node_idx = ret_graph.add_node(other_node.to_owned());
                stack.push((old_other_node_idx, new_node_idx));
                let old_edge = old_edge_ref.weight().to_owned();

                if direction == Direction::Outgoing {
                    ret_graph.add_edge(new_current_idx, new_node_idx, old_edge.to_owned());
                } else {
                    ret_graph.add_edge(new_node_idx, new_current_idx, old_edge.to_owned());
                };
            }
        }
    }

    Ok(ret_graph)
}

#[cfg(test)]
mod tests {

    use petgraph::Graph;

    use super::*;

    #[test]
    fn parse_raw_test() {
        let match_pred: SubGraphPred<Option<u32>> =
            |a: &Option<u32>, b: &Option<u32>| -> bool {
                let Some(a) = a
                else {
                    return true;
                };
                let Some(b) = b
                else {
                    return true;
                };
                a == b
            };

        let mut graph = Graph::<Option<u32>, &'static str>::new();

        let n0 = graph.add_node(None);
        let n1 = graph.add_node(Some(10));
        let n2 = graph.add_node(None);
        let n3 = graph.add_node(Some(10));
        let n4 = graph.add_node(Some(10));
        let n5 = graph.add_node(None);
        let n6 = graph.add_node(Some(20));
        let n7 = graph.add_node(Some(30));

        graph.add_edge(n0, n1, "n0-n1");
        graph.add_edge(n1, n2, "n1-n2");
        graph.add_edge(n2, n3, "n2-n3");
        graph.add_edge(n3, n5, "n4-n5");
        graph.add_edge(n5, n6, "n5-n6");
        graph.add_edge(n2, n4, "n2-n4");
        graph.add_edge(n4, n7, "n4-n7");

        let sub10 = subgraph_by(&graph, n3, &match_pred).unwrap();
        let sub20 = subgraph_by(&graph, n6, &match_pred).unwrap();
        let sub30 = subgraph_by(&graph, n7, &match_pred).unwrap();
        let subnone = subgraph_by(&graph, n5, &match_pred).unwrap();

        assert_eq!(sub10.node_count(), 6);
        assert_eq!(sub20.node_count(), 2);
        assert_eq!(sub30.node_count(), 1);
        assert_eq!(subnone.node_count(), 8);
    }
}

