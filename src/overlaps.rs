use std::collections::{BTreeSet, HashMap, HashSet};

use petgraph::{graph::{EdgeReference, NodeIndex}, visit::EdgeRef, Direction};

use anyhow::Result;

use crate::{graph::{self, SubGraphPred}, raw::IndividualId, GoCamActivity, GoCamComponent, GoCamDirection, GoCamEdge, GoCamEnabledBy, GoCamGraph, GoCamModel, GoCamModelId, GoCamModelIdTitle, GoCamNode, GoCamNodeType, GoCamProcess};

/// An overlap returned by [GoCamModel::find_overlaps()]
#[derive(Serialize, Deserialize, Clone, Debug)]
#[serde(rename_all = "snake_case")]
pub struct GoCamNodeOverlap {
    pub node_id: String,
    pub node_label: String,
    pub node_type: GoCamNodeType,

    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub occurs_in: BTreeSet<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_process: Option<GoCamProcess>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub happens_during: Option<GoCamProcess>,
    pub overlapping_individual_ids: BTreeSet<IndividualId>,

    // the "home" model that this activity comes from
    // it will be added to other models to help with joining
    pub original_model_id: Option<GoCamModelId>,

    // a set of the model details for this overlap, with the direction of relations
    // into/out of the node in the given model
    pub models: BTreeSet<GoCamModelIdTitle>,
}

impl GoCamNodeOverlap {
    /// A unique ID for this overlap, created from the Individual IDs
    /// of the nodes/activities
    pub fn id(&self) -> String {
        self.overlapping_individual_ids.iter().cloned().collect::<Vec<_>>().join("-")
    }
}


/// Return the overlaps between models.  A [GoCamNodeOverlap] is
/// returned for each pair of models that have an activity in
/// common.  The pair of activities must have the same MF term, be
/// enabled by the same entity (gene, complex, modified protein or
/// chemical), have the same process and same the component
/// ("occurs in").  The process and component must be non-None.
pub fn find_overlaps(models: &[GoCamModel])
   -> Vec<GoCamNodeOverlap>
{
    let make_key = |node: &GoCamNode| {
        if node.occurs_in.is_empty() ||
            node.part_of_process.is_none() {
                return None;
            }

        let GoCamNodeType::Activity(GoCamActivity { enabler: ref enabled_by, .. }) = node.node_type
        else {
            return None;
        };

        Some((node.node_id.clone(),
              node.label.clone(),
              enabled_by.clone(),
              node.part_of_process.clone().unwrap(),
              node.happens_during.clone(),
              node.occurs_in.clone()))
    };

    let mut models_by_id = HashMap::new();

    let mut activities = HashMap::new();

    // put all activities in a HashMap
    for model in models {
        models_by_id.insert(model.id().to_owned(), model);

        for (node_idx, node) in model.node_iterator() {
            let Some(key) = make_key(node)
            else {
                continue;
            };

            activities
                .entry(key)
                .or_insert_with(HashMap::new)
                .entry(model.id())
                .or_insert_with(HashSet::new)
                .insert((node_idx, node));
        }
    }

    // From the HashMap of activities, find those that are in more than one model.  We ignore
    // activities that occur more than once in a model as we can't handle that
    let possible_overlapping_activities: HashMap<_, _> =
        activities
        .iter()
        .filter_map(|(key, node_details)| {
            if node_details.len() < 2 {
                // there is no overlap
                return None;
            }

            let mut models_and_nodes = vec![];

            for (model_id, nodes) in node_details.into_iter() {
                if nodes.len() == 1 {
                    let (node_idx, node) = nodes.iter().next().unwrap().clone();
                    models_and_nodes.push((model_id.to_owned(), node_idx, node));
                } else {
                    // for now ignore cases where an activity is duplicated in a model
                    return None;
                }
            }

            Some((key, models_and_nodes))
        })
        .collect();

    let mut ret = vec![];

    let mut overlapping_nodes_by_model = HashMap::new();

    for models_and_individual in possible_overlapping_activities.values() {
        for (model_id, _, node) in models_and_individual.iter() {
            overlapping_nodes_by_model
                .entry(model_id)
                .or_insert_with(HashSet::new)
                .insert(*node);
        }
    }

    // check that inputs or outputs match
    for (key, models_and_individuals) in possible_overlapping_activities.iter() {
        let (node_id, node_label, enabled_by, part_of_process, happens_during, occurs_in) = key;

        let activities_inputs_outputs: Vec<_> = models_and_individuals
            .iter()
            .map(|(_, _, node)| {
                let GoCamNodeType::Activity(GoCamActivity { ref enabler, ref inputs, ref outputs }) = node.node_type
                else {
                    panic!("internal error: expected an activity");
                };

                (enabler.clone(), inputs.clone(), outputs.clone())
            })
            .collect();



        let first_inputs = activities_inputs_outputs.first()
            .map(|(_, inputs, _)| inputs.clone()).unwrap();
        let first_outputs = activities_inputs_outputs.first()
            .map(|(_, _, outputs)| outputs.clone()).unwrap();

        let has_matching_inputs_or_outputs =
            activities_inputs_outputs.iter()
            .all(|(_, inputs, _)| *inputs == first_inputs)
            ||
            activities_inputs_outputs.iter()
            .all(|(_, _, outputs)| *outputs == first_outputs);

        if !has_matching_inputs_or_outputs {
            // nodes don't have all the same inputs and the nodes don't have all the same
            // outputs
            continue;
        }

        let mut overlapping_activity_nodes = HashSet::new();
        let mut models_with_complete_processes: Vec<(GoCamModelId, GoCamNode)> = vec![];

        let mut model_ids_and_titles = BTreeSet::new();

        let mut possible_input_output_overlaps = HashMap::new();

        let mut original_model_id = None;

        for (model_id, node_idx, node) in models_and_individuals.clone() {
            let &model = models_by_id.get(model_id).unwrap();

            let direction = node_rel_direction(&model.graph(), node_idx);

            model_ids_and_titles.insert((model_id.to_owned(), model.title().to_owned(),
                                         direction));

            overlapping_activity_nodes.insert(node);

            for input_output in inputs_outputs_of(model, node_idx) {
                let (ref input_output_edge, ref input_output_node, node_idx) = input_output;
                let key = (input_output_edge.id.clone(),
                           input_output_node.node_type.clone());

                let val = (model_id, model.title().to_owned(),
                           input_output_node.clone(), node_idx);

                possible_input_output_overlaps.entry(key)
                    .or_insert_with(Vec::new)
                    .push(val);
            }

            let this_model_overlaps_ids =
                overlapping_nodes_by_model.get(&model_id).unwrap()
                .iter().map(|n| n.individual_gocam_id.to_owned()).collect();

            if process_subgraph_in_overlap(model, node_idx, &this_model_overlaps_ids) {
                models_with_complete_processes.push((model.id().to_owned(), node.clone()));
            } else {
                original_model_id = Some(model.id().to_owned());
            }
        }

        if original_model_id.is_none() {
            let mut found_original_model_ids = HashSet::new();
            for (model_id, _, node) in models_and_individuals.clone() {
                let Some(ref process) = node.part_of_process
                else {
                    continue;
                };
                let Some(ref model) = models_by_id.get(model_id)
                else {
                    continue;
                };
                if model.title_process_term_ids.contains(&process.id) {
                    found_original_model_ids.insert(model_id.to_owned());
                }
                let Some(ref part_of_parent) = process.part_of_parent
                else {
                    continue;
                };
                if model.title_process_term_ids.contains(&part_of_parent.id) {
                    found_original_model_ids.insert(model_id.to_owned());
                }
            }
            if found_original_model_ids.len() == 1 {
                original_model_id = found_original_model_ids.iter().cloned().next();
            }
        }

        if models_with_complete_processes.is_empty() {
            continue;
        }

        let has_incoming = model_ids_and_titles.iter()
            .any(|(_, _, direction)| *direction == GoCamDirection::Incoming);
        let has_incoming_constitutively_upstream = model_ids_and_titles.iter()
            .any(|(_, _, direction)| *direction == GoCamDirection::IncomingConstitutivelyUpstream);
        let has_outgoing = model_ids_and_titles.iter()
            .any(|(_, _, direction)| *direction == GoCamDirection::Outgoing);
        let has_none = model_ids_and_titles.iter()
            .any(|(_, _, direction)| *direction == GoCamDirection::None);

        // See: https://github.com/pombase/pombase-gocam/issues/36#issuecomment-2982741903
        if !has_incoming_constitutively_upstream && // special case for complex assembly
            (!has_incoming && !has_outgoing ||
             has_incoming && !(has_outgoing || has_none) ||
             has_outgoing && !(has_incoming || has_none)) {
                // remove cases where there is no clear direction between models
                continue;
            }

        let enabled_by =
            if let GoCamEnabledBy::Complex(complex) = enabled_by {
                let mut complex = complex.to_owned();
                // the complex used in the key may not have all the has_part_genes
                for (_, _, node) in models_and_individuals {
                    if let GoCamNodeType::Activity(GoCamActivity {
                        enabler: GoCamEnabledBy::Complex(ref this_complex),
                        ..
                    }) = node.node_type {
                        for this_gene in &this_complex.has_part_genes {
                            complex.has_part_genes.insert(this_gene.to_owned());
                        }
                    }
                }
                GoCamEnabledBy::Complex(complex)

            } else {
                enabled_by.to_owned()
            };


        let overlapping_individual_ids = overlapping_activity_nodes.iter()
            .map(|n| n.individual_gocam_id.to_owned())
            .collect();
        let inputs = overlapping_activity_nodes.iter()
            .flat_map(|n| {
                let GoCamNodeType::Activity(GoCamActivity { enabler: _, ref inputs, outputs: _ }) = n.node_type
                else {
                    panic!("internal error: expected an activity");
                };
                inputs.to_owned()
            })
            .collect();
        let outputs = overlapping_activity_nodes.iter()
            .flat_map(|n| {
                let GoCamNodeType::Activity(GoCamActivity { enabler: _, inputs: _, ref outputs }) = n.node_type
                else {
                    panic!("internal error: expected an activity");
                };
                outputs.to_owned()
            })
            .collect();

        let node_overlap = GoCamNodeOverlap {
            node_id: node_id.to_owned(),
            node_label: node_label.to_owned(),
            node_type: GoCamNodeType::Activity(GoCamActivity { enabler: enabled_by, inputs, outputs, }),
            part_of_process: Some(part_of_process.to_owned()),
            occurs_in: occurs_in.to_owned(),
            happens_during: happens_during.to_owned(),
            overlapping_individual_ids,
            original_model_id,
            models: model_ids_and_titles,
        };

        ret.push(node_overlap);

        for (input_output_key, input_output_details) in possible_input_output_overlaps {
            if input_output_details.len() == 1 {
                // no overlap
                continue;
            }

            let (_, node_type) = input_output_key;
            let mut overlapping_individual_ids = BTreeSet::new();
            let mut model_ids_and_titles = BTreeSet::new();

            let (_, _, first_node, _) = input_output_details.get(0).unwrap().to_owned();

            let mut original_model_ids: BTreeSet<_> = BTreeSet::new();

            for (model_id, model_title, node, node_idx) in input_output_details {
                overlapping_individual_ids.insert(node.individual_gocam_id);
                model_ids_and_titles.insert((model_id.to_owned(), model_title, GoCamDirection::None));

                // Try to choose an original/home model ID for this chemical by using the
                // original_model_id the model where the chemical has an upstream and a
                // downstream activity.
                let model = models_by_id.get(model_id).unwrap();

                // in this case these are activities:
                let inputs_outputs = inputs_outputs_of(*model, node_idx);

                let has_upstream_activity = inputs_outputs
                    .iter()
                    .any(|(edge, _, _)| edge.id == "RO:0002234");

                let has_downstream_activity = inputs_outputs
                    .iter()
                    .any(|(edge, _, _)| edge.id == "RO:0002233");

                if has_upstream_activity && has_downstream_activity {
                    original_model_ids.insert(model_id.to_owned());
                }
            }

            let original_model_id =
                if original_model_ids.len() == 1 {
                    original_model_ids.first().map(String::to_owned)
                } else {
                    None
                };

            let input_output_node_overlap = GoCamNodeOverlap {
                node_id: node_type.id().to_owned(),
                node_label: node_type.label().to_owned(),
                node_type: first_node.node_type,
                part_of_process: None,
                occurs_in: BTreeSet::new(),
                happens_during: None,
                overlapping_individual_ids,
                original_model_id,
                models: model_ids_and_titles,
            };

            ret.push(input_output_node_overlap);
        }
    }

    ret.sort_by(|a, b| {
        a.models.cmp(&b.models)
            .then(a.node_label.cmp(&b.node_label))
    });

    ret
}

fn inputs_outputs_of(model: &GoCamModel, activity_index: NodeIndex)
      -> Vec<(GoCamEdge, GoCamNode, NodeIndex)>
{
    let mut ret = vec![];

    let graph = model.graph();

    let outgoing_iter = graph.edges_directed(activity_index, Direction::Outgoing);

    for edge_ref in outgoing_iter {
        let target_idx = edge_ref.target();
        let target_node = graph.node_weight(target_idx).unwrap();

        let edge = edge_ref.weight();

        if edge.id != "RO:0002234" && edge.id != "RO:0002233" {
            continue;
        }

        ret.push((edge.to_owned(), target_node.to_owned(), target_idx))
    }

    let incoming_iter = model.graph().edges_directed(activity_index, Direction::Incoming);

    for edge_ref in incoming_iter {
        let subject_idx = edge_ref.target();
        let subject_node = graph.node_weight(subject_idx).unwrap();

        let edge = edge_ref.weight();

        if edge.id != "RO:0002234" && edge.id != "RO:0002233" {
            continue;
        }

        ret.push((edge.to_owned(), subject_node.to_owned(), subject_idx))
    }

    ret
}

    // return the sub-graph of `graph` that includes `start_idx` and has the
    // same process as `start_idx`
    fn subgraph_by_process(graph: &GoCamGraph, start_idx: NodeIndex)
       -> Result<GoCamGraph>
    {
        let same_process: SubGraphPred<GoCamNode> =
        |a: &GoCamNode, b: &GoCamNode| -> bool {
            let Some(ref a_process) = a.part_of_process
            else {
                return true;
            };
            let Some(ref b_process) = b.part_of_process
            else {
                return true;
            };
            a_process == b_process
        };

        graph::subgraph_by(graph, start_idx, &same_process)
    }

fn process_subgraph_in_overlap(model: &GoCamModel, start_idx: NodeIndex,
                               id_overlaps: &HashSet<IndividualId>)
                               -> bool
{
    let sub_graph = match subgraph_by_process(&model.graph, start_idx) {
        Ok(sub_graph) => sub_graph,
        Err(_) => {
            return false;
        }
    };

    for node in sub_graph.node_weights() {
        if node.is_activity() && !id_overlaps.contains(&node.individual_gocam_id) {
            return false;
        }
    }

    true
}

// If the activity at `activity_idx` has only incoming relations (excluding inputs/outputs),
// return Incoming, unless all the incoming relations are "constitutively upstream" in which
// case return IncomingConstitutivelyUpstream.
// If the activity has only outgoing relations return Outgoing.
// Otherwise return None.
fn node_rel_direction(graph: &GoCamGraph, activity_idx: NodeIndex) -> GoCamDirection {
    let is_causal_edge = |edge_ref:  EdgeReference<GoCamEdge>| {
        let edge = edge_ref.weight();
        edge.id != "RO:0002234" && edge.id != "RO:0002233"
    };

    let is_constitutively_upstream_edge = |edge_ref:  EdgeReference<GoCamEdge>| {
        let edge = edge_ref.weight();
        edge.id == "RO:0012009"
    };

    let mut has_incoming = graph.edges_directed(activity_idx, Direction::Incoming)
        .into_iter()
        .any(is_causal_edge);


    let all_incoming_edge_constitutively_upstream = {
        // true if there are some incoming edges and the are all
        // "constitutively upstream" rels
        graph.edges_directed(activity_idx, Direction::Incoming).next().is_some() &&
            graph.edges_directed(activity_idx, Direction::Incoming).into_iter()
            .all(is_constitutively_upstream_edge)
    };

    let mut has_outgoing = graph.edges_directed(activity_idx, Direction::Outgoing)
        .into_iter()
        .any(is_causal_edge);

    if all_incoming_edge_constitutively_upstream && !has_outgoing {
        let node = graph.node_weight(activity_idx).unwrap();
        if node.type_string() == "enabled_by_complex" {
            // special case for activities enable by a complex that
            // only has "constitutively upstream" incoming relations
            return GoCamDirection::IncomingConstitutivelyUpstream;
        }
    }

    if !has_incoming || !has_outgoing {

        for edge in graph.edges_directed(activity_idx, Direction::Outgoing) {
            let is_input_edge =
                if edge.weight().id == "RO:0002233" {
                    true
                } else {
                    if edge.weight().id == "RO:0002234" {
                        false
                    } else {
                        continue;
                    }
                };
            let target_idx = edge.target();
            for incoming_edge in graph.edges_directed(target_idx, Direction::Incoming) {
                if is_input_edge && incoming_edge.weight().id == "RO:0002234" {
                    has_incoming = true;
                    break;
                }
                if !is_input_edge && incoming_edge.weight().id == "RO:0002233" {
                    has_outgoing = true;
                    break;
                }
            }
        }
    }

    match (has_incoming, has_outgoing) {
        (true, true) => GoCamDirection::None,
        (false, false) => GoCamDirection::None,
        (true, false) => GoCamDirection::Incoming,
        (false, true) => GoCamDirection::Outgoing,
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use crate::{overlaps::{find_overlaps, node_rel_direction}, parse_gocam_model, GoCamDirection};

    #[test]
    fn find_overlaps_test() {
        let mut source1 = File::open("tests/data/gomodel_66a3e0bb00001342.json").unwrap();
        let model1 = parse_gocam_model(&mut source1).unwrap();
        assert_eq!(model1.id(), "gomodel:66a3e0bb00001342");

        assert_eq!(model1.node_iterator().count(), 29);

        let mut source2 = File::open("tests/data/gomodel_665912ed00000015.json").unwrap();
        let model2 = parse_gocam_model(&mut source2).unwrap();
        assert_eq!(model2.id(), "gomodel:665912ed00000015");

        assert_eq!(model2.node_iterator().count(), 25);

        let mut source3 = File::open("tests/data/gomodel_678073a900003175.json").unwrap();
        let model3 = parse_gocam_model(&mut source3).unwrap();
        assert_eq!(model3.id(), "gomodel:678073a900003175");

        assert_eq!(model3.node_iterator().count(), 13);

        let overlaps = find_overlaps(&[model1, model2, model3]);

        assert_eq!(overlaps.len(), 2);

        let overlap_activity = &overlaps[0];

        assert_eq!(overlap_activity.node_label, "homoserine O-acetyltransferase activity");
        assert_eq!(overlap_activity.models.len(), 2);

        assert_eq!(overlap_activity.part_of_process.as_ref().unwrap().id, "GO:0071266");
        assert_eq!(overlap_activity.part_of_process.as_ref().unwrap().label,
                   "'de novo' L-methionine biosynthetic process");
        let first_occurs_in = overlap_activity.occurs_in.iter().next().unwrap();
        assert_eq!(first_occurs_in.id(), "GO:0005829");
        assert_eq!(first_occurs_in.label(), "cytosol");

        let first_overlapping_individual =
            overlap_activity.overlapping_individual_ids.iter().next().unwrap();
        assert_eq!(first_overlapping_individual,
                   "gomodel:66a3e0bb00001342/678073a900003752");

        let overlap_chemical = &overlaps[1];
        assert_eq!(overlap_chemical.node_label, "O-acetyl-L-homoserine");
    }

    #[test]
    fn test_node_rel_direction() {
        let mut source = File::open("tests/data/gomodel_67f85f2b00002766.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:67f85f2b00002766");

        let Some ((test_out_node_idx, _)) = model.node_iterator()
            .find(|(_, node)| {
                node.enabler_id() == "GO:0061671"
            })
        else {
            panic!()
        };

        // this a complex at the bottom of an assembly model
        let Some ((test_in_node_idx, _)) = model.node_iterator()
            .find(|(_, node)| {
                node.enabler_id() == "GO:0045275"
            })
        else {
            panic!()
        };

        let test_out_direction = node_rel_direction(model.graph(), test_out_node_idx);
        assert_eq!(test_out_direction, GoCamDirection::Outgoing);

        let test_in_direction = node_rel_direction(model.graph(), test_in_node_idx);
        assert_eq!(test_in_direction, GoCamDirection::IncomingConstitutivelyUpstream);
    }

}
