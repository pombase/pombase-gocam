//! A library for parsing and processing [GO-CAM](https://geneontology.org/docs/gocam-overview)
//! JSON format model files
//!
//! [GoCamRawModel] is a low level representation which closely
//! matches the JSON data, containing Fact, Individual and Annotation
//! structs.
//!
//! [GoCamModel] is a higher level representation, implemented as a
//! graph of nodes (activities, chemical, complexes etc.) and edges
//! (mostly causal relations).
//!
//! This high level representation is somewhat similar to the
//! [GO CAM Data Model - gocam-py](https://github.com/geneontology/gocam-py).
//!
//! ## Example
//!
//! ```
//! use std::fs::File;
//! use pombase_gocam::raw::gocam_parse_raw;
//!
//! let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
//!
//! // Low level representation:
//! let raw_model = gocam_parse_raw(&mut source).unwrap();
//! assert_eq!(raw_model.id(), "gomodel:66187e4700001744");
//!
//! for fact in raw_model.facts() {
//!     let subject_id = &fact.subject;
//!     println!("subject_id: {}", subject_id);
//!     let subject_individual = raw_model.get_individual(subject_id);
//!     let first_type = &subject_individual.types[0];
//!     if let Some(ref label) = first_type.label {
//!         println!("type label: {}", label);
//!     }
//! }
//!
//! // Higher level representation:
//! use pombase_gocam::{GoCamModel, GoCamNodeType};
//! let model = GoCamModel::new(raw_model);
//!
//! for (_, node) in model.node_iterator() {
//!     println!("node: {}", node);
//!     if let GoCamNodeType::Activity(ref enabler) = node.node_type {
//!         println!("enabler ID: {}", enabler.id());
//!     }
//! }
//!```

use std::{collections::{BTreeMap, BTreeSet, HashMap, HashSet},
          fmt::{self, Display},
          io::Read};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

pub mod raw;
pub mod graph;

use graph::SubGraphPred;
use raw::{gocam_parse_raw, FactId, GoCamRawModel, Individual, IndividualId, IndividualType};

use phf::phf_map;

use anyhow::{Result, anyhow};

use petgraph::{graph::{NodeIndex, NodeReferences}, visit::{EdgeRef, IntoNodeReferences}, Direction, Graph};

/// A map of edge relation term IDs to term names.  Example:
/// "RO:0002211" => "regulates",
pub static REL_NAMES: phf::Map<&'static str, &'static str> = phf_map! {
    "BFO:0000050" => "part of",
    "BFO:0000051" => "has part",
    "RO:0002233" => "has input",
    "RO:0002234" => "has output",
    "RO:0002333" => "enabled by",
    "RO:0002413" => "provides input for",
    "BFO:0000066" => "occurs in",
    "RO:0001025" => "located in",
    "RO:0002092" => "happens during",
    "RO:0002407" => "indirectly positively regulates",
    "RO:0002409" => "indirectly negatively regulates",
    "RO:0002629" => "directly positively regulates",
    "RO:0002630" => "directly negatively regulates",
    "RO:0002304" => "causally upstream of, positive effect",
    "RO:0002305" => "causally upstream of, negative effect",
    "RO:0012005" => "is small molecule activator of",
    "RO:0012006" => "is small molecule inhibitor of",
    "RO:0012009" => "constitutively upstream of",
    "RO:0000057" => "has participant",
    "RO:0001015" => "location of",
    "RO:0002020" => "transports",
    "RO:0002131" => "overlaps",
    "RO:0002162" => "in taxon",
    "RO:0002211" => "regulates",
    "RO:0002212" => "negatively regulates",
    "RO:0002213" => "positively regulates",
    "RO:0002220" => "adjacent to",
    "RO:0002221" => "surrounds",
    "RO:0002224" => "starts with",
    "RO:0002230" => "ends with",
    "RO:0002232" => "has end location",
    "RO:0002264" => "acts upstream of or within",
    "RO:0002295" => "results in developmental progression of",
    "RO:0002296" => "results in development of",
    "RO:0002297" => "results in formation of anatomical entity",
    "RO:0002298" => "results in morphogenesis of",
    "RO:0002299" => "results in maturation of",
    "RO:0002313" => "transports or maintains localization of",
    "RO:0002315" => "results in acquisition of features of",
    "RO:0002325" => "colocalizes with",
    "RO:0002326" => "contributes to",
    "RO:0002331" => "involved in",
    "RO:0002338" => "has target start location",
    "RO:0002339" => "has target end location",
    "RO:0002348" => "results in commitment to",
    "RO:0002349" => "results in determination of",
    "RO:0002356" => "results in specification of",
    "RO:0002400" => "has direct input",
    "RO:0002406" => "obsolete directly activates",
    "RO:0002411" => "causally upstream of",
    "RO:0002412" => "immediately causally upstream of",
    "RO:0002418" => "causally upstream of or within",
    "RO:0002432" => "is active in",
    "RO:0002490" => "existence overlaps",
    "RO:0002491" => "existence starts and ends during",
    "RO:0002502" => "depends on",
    "RO:0002565" => "results in movement of",
    "RO:0002578" => "directly regulates",
    "RO:0002588" => "results in assembly of",
    "RO:0002590" => "results in disassembly of",
    "RO:0002592" => "results in organization of",
    "RO:0004008" => "has primary output",
    "RO:0004009" => "has primary input",
    "RO:0004033" => "acts upstream of or within, negative effect",
    "RO:0004035" => "acts upstream of, negative effect",
    "RO:0004046" => "causally upstream of or within, negative effect",
    "RO:0004047" => "causally upstream of or within, positive effect",
    "RO:0012001" => "has small molecule activator",
    "RO:0012002" => "has small molecule inhibitor",
    "RO:0012003" => "acts on population of",
    "RO:0012010" => "removes input for",
};

pub type ModelId = String;
pub type ModelTitle = String;

/// The Graph representation of a GO-CAM model.  See: [GoCamModel::graph()]
pub type GoCamGraph = Graph::<GoCamNode, GoCamEdge>;

/// A gene ID with DB prefix, like "PomBase:SPAC9E9.05"
pub type GoCamGeneIdentifier = String;

type GoCamNodeMap = BTreeMap<IndividualId, GoCamNode>;

/// A high level representation of the model with nodes for
/// activities, chemicals, complexes etc. and edges for causal
/// dependencies between nodes/activities.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamModel {
    id: String,
    title: String,
    taxon: String,
    date: String,
    contributors: BTreeSet<String>,
    graph: GoCamGraph,
}

fn check_model_taxons(models: &[GoCamModel]) -> Result<String> {
    let mut model_iter = models.iter();

    let Some(first_model) = model_iter.next()
    else {
        return Err(anyhow!("no models passed to check_model_taxons()"));
    };

    for this_model in model_iter {
        // a hack to cope with models with more than one taxon:
        // "NCBITaxon:4896,NCBITaxon:559292"
        if !first_model.taxon.contains(&this_model.taxon) &&
           !this_model.taxon.contains(&first_model.taxon) {
            return Err(anyhow!("mismatched taxons: {} ({}) != {} ({})",
                               first_model.taxon, first_model.id,
                               this_model.taxon, this_model.id));
        }
    }

    Ok(first_model.taxon().to_owned())
}

impl GoCamModel {
    /// Create a [GoCamModel] from a [GoCamRawModel]
    ///
    /// ## Example
    /// ```
    /// use pombase_gocam::GoCamModel;
    /// let mut source = std::fs::File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
    /// let raw_model = pombase_gocam::raw::gocam_parse_raw(&mut source).unwrap();
    /// let model = pombase_gocam::GoCamModel::new(raw_model);
    /// ```
    ///
    /// See also [parse_gocam_model()].
    pub fn new(raw_model: GoCamRawModel) -> GoCamModel {
        let graph = make_graph(&raw_model);

        GoCamModel {
            id: raw_model.id().to_owned(),
            title: raw_model.title().to_owned(),
            taxon: raw_model.taxon().to_owned(),
            date: raw_model.date().to_owned(),
            contributors: raw_model.contributors(),
            graph,
        }
    }

    /// Return the [petgraph::Graph] representation of the model
    pub fn graph(&self) -> &GoCamGraph {
        &self.graph
    }

    /// Return an iterator over the nodes ([GoCamNode]) of the model
    pub fn node_iterator(&self) -> NodeIterator {
        NodeIterator {
            node_refs: self.graph().node_references(),
        }
    }

    /// Number of nodes in the graph
    pub fn node_count(&self) -> usize {
        self.graph().node_count()
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn title(&self) -> &str {
        &self.title
    }

    /// The taxon ID in the format "NCBITaxon:4896" or possible a
    /// comma separated list like: "NCBITaxon:4896,NCBITaxon:559292"
    pub fn taxon(&self) -> &str {
        &self.taxon
    }

    /// The date the model last changed
    pub fn date(&self) -> &str {
        &self.date
    }

    #[allow(rustdoc::bare_urls)]
    /// A set of contributor ORCIDs like:
    /// "https://orcid.org/0000-0001-6330-7526"
    pub fn contributors(&self) -> &BTreeSet<String> {
        &self.contributors
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
        let make_key = |node: &GoCamNode, input_output: &str| {
            if node.occurs_in.is_none() ||
                node.part_of_process.is_none() {
                return None;
            }

            let GoCamNodeType::Activity(ref enabled_by) = node.node_type
            else {
                return None;
            };

            let input =
                if input_output == "input" {
                    node.has_input.clone()
                } else {
                    BTreeSet::new()
                };

            let output =
                if input_output == "output" {
                    node.has_output.clone()
                } else {
                    BTreeSet::new()
                };

            Some((node.node_id.clone(),
                  node.label.clone(),
                  enabled_by.clone(),
                  node.part_of_process.clone().unwrap(),
                  node.occurs_in.clone().unwrap(),
                  input, output))
        };

        let mut models_by_id = HashMap::new();

        let mut activities = HashMap::new();

        for model in models {
            models_by_id.insert(model.id().to_owned(), model);

            for (node_idx, node) in model.node_iterator() {
                let Some(key) = make_key(node, "input")
                else {
                    continue;
                };

                activities
                    .entry(key)
                    .or_insert_with(HashMap::new)
                    .entry(model.id())
                    .or_insert_with(HashSet::new)
                    .insert((node_idx, node));

                let Some(key) = make_key(node, "output")
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

        let possible_overlapping_activities: HashMap<_, _> =
            activities
            .into_iter()
            .filter_map(|(key, node_details)| {
                if node_details.len() < 2 {
                    // there is no overlap
                    None
                } else {
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
                }
            })
            .collect();

        let mut ret = vec![];

        let mut overlapping_nodes_by_model = HashMap::new();

        for models_and_individual in possible_overlapping_activities.values() {
            for (model_id, _, node) in models_and_individual.iter() {
                overlapping_nodes_by_model
                    .entry(model_id.clone())
                    .or_insert_with(HashSet::new)
                    .insert(*node);
            }
        }

        for (key, models_and_individual) in possible_overlapping_activities.into_iter() {
            let (node_id, node_label, enabled_by,
                 part_of_process, occurs_in, _, _) = key;

            let mut overlapping_individual_ids = BTreeSet::new();
            let mut found_complete_process = false;

            let mut model_ids_and_titles = BTreeSet::new();

            let mut possible_chemical_overlaps = HashMap::new();

            for (model_id, node_idx, node) in models_and_individual.clone() {
                let &model = models_by_id.get(&model_id).unwrap();

                model_ids_and_titles.insert((model_id.clone(), model.title().to_owned()));

                overlapping_individual_ids.insert(node.individual_gocam_id.to_owned());

                for chemical_neighbour in Self::chemical_neighbours_of(model, node_idx) {
                    let (ref chemical_neighbour_edge, ref chemical_neighbour_node) = chemical_neighbour;
                    let key = (chemical_neighbour_edge.id.clone(),
                               chemical_neighbour_node.node_id.clone(),
                               chemical_neighbour_node.label.clone(),
                               chemical_neighbour_node.located_in.clone());

                    let val = (model_id.clone(), model.title().to_owned(),
                               chemical_neighbour_node.clone());

                    possible_chemical_overlaps.entry(key)
                        .or_insert_with(Vec::new)
                        .push(val);
                }

                let this_model_overlaps_ids =
                    overlapping_nodes_by_model.get(&model_id).unwrap()
                    .iter().map(|n| n.individual_gocam_id.to_owned()).collect();

                if Self::is_process_sub_graph(model, node_idx, &this_model_overlaps_ids) {
                    found_complete_process = true;
                }
            }

            if !found_complete_process {
                continue;
            }

            let node_overlap = GoCamNodeOverlap {
                node_id: node_id.to_owned(),
                node_label,
                node_type: GoCamNodeType::Activity(enabled_by),
                has_input: vec![],
                has_output: vec![],
                part_of_process: Some(part_of_process),
                occurs_in: Some(occurs_in),
                located_in: None,
                overlapping_individual_ids,
                models: model_ids_and_titles,
            };

            ret.push(node_overlap);

            for (chemical_key, chemical_details) in possible_chemical_overlaps {
                if chemical_details.len() == 1 {
                    // no overlap
                    continue;
                }
                let (_, node_id, node_label, located_in) = chemical_key;
                let mut overlapping_individual_ids = BTreeSet::new();
                let mut model_ids_and_titles = BTreeSet::new();
                for (model_id, model_title, chemical_node) in chemical_details {
                    overlapping_individual_ids.insert(chemical_node.individual_gocam_id);
                    model_ids_and_titles.insert((model_id, model_title));
                }

                let node_overlap = GoCamNodeOverlap {
                        node_id,
                        node_label,
                        node_type: GoCamNodeType::Chemical,
                        has_input: vec![],
                        has_output: vec![],
                        part_of_process: None,
                        occurs_in: None,
                        located_in: located_in,
                        overlapping_individual_ids,
                        models: model_ids_and_titles,
                    };

                    ret.push(node_overlap);
            }
        }

        ret.sort_by(|a, b| {
            a.models.cmp(&b.models)
                .then(a.node_label.cmp(&b.node_label))
        });

        ret
    }

    fn chemical_neighbours_of(model: &GoCamModel, activity_index: NodeIndex)
       -> Vec<(GoCamEdge, GoCamNode)>
    {
        let mut ret = vec![];

        let graph = model.graph();

        let outgoing_iter = graph.edges_directed(activity_index, Direction::Outgoing);

        for edge_ref in outgoing_iter {
            let target_node = graph.node_weight(edge_ref.target()).unwrap();

            if target_node.node_type != GoCamNodeType::Chemical {
                continue;
            }

            ret.push((edge_ref.weight().to_owned(), target_node.to_owned()))
        }


        let incoming_iter = model.graph().edges_directed(activity_index, Direction::Incoming);

        for edge_ref in incoming_iter {
            let subject_node = graph.node_weight(edge_ref.target()).unwrap();

            if subject_node.node_type != GoCamNodeType::Chemical {
                continue;
            }

            ret.push((edge_ref.weight().to_owned(), subject_node.to_owned()))
        }

        ret
    }

    fn is_process_sub_graph(model: &GoCamModel, start_idx: NodeIndex, id_overlaps: &HashSet<IndividualId>)
       -> bool
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

        let sub_graph = match graph::subgraph_by(&model.graph, start_idx, &same_process) {
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


    /// Merge the `models` that have nodes in common, returning a new
    /// [GoCamModel] with the `new_id` as the ID and `new_title` as
    /// the title.
    ///
    /// We use the result of calling [Self::find_overlaps()] to find
    /// nodes in common between all the `models`.
    pub fn merge_models(new_id: &str, new_title: &str, models: &[GoCamModel])
        -> Result<GoCamModel>
    {
        let mut merged_graph = GoCamGraph::new();

        let overlaps = Self::find_overlaps(models);

        let mut overlap_map = HashMap::new();

        for overlap in overlaps.into_iter() {
            let overlap_id = overlap.id();

            let overlap_node = GoCamNode {
                individual_gocam_id: overlap_id,
                node_id: overlap.node_id,
                label: overlap.node_label,
                node_type: overlap.node_type,
                has_input: BTreeSet::new(),
                has_output: BTreeSet::new(),
                occurs_in: overlap.occurs_in,
                part_of_process: overlap.part_of_process,
                located_in: overlap.located_in,
                source_ids: overlap.overlapping_individual_ids.clone(),
                models: overlap.models.clone(),
            };

            let overlap_node_idx = merged_graph.add_node(overlap_node);

            for overlapping_individual in overlap.overlapping_individual_ids.iter() {
                overlap_map.insert(overlapping_individual.clone(), overlap_node_idx);
            }
        }

        let mut idx_map = HashMap::new();
        let taxon = check_model_taxons(models)?;

        let mut contributors = BTreeSet::new();

        let mut seen_edge = HashSet::new();

        for model in models {
            for (old_idx, node) in model.graph().node_references() {
                if let Some(overlap_node_idx) = overlap_map.get(&node.individual_gocam_id) {
                    idx_map.insert(old_idx, *overlap_node_idx);
                } else {
                    let new_idx = merged_graph.add_node(node.to_owned());
                    idx_map.insert(old_idx, new_idx);
                }
            }

            for edge_ref in model.graph().edge_references() {
                let edge = edge_ref.weight();

                let old_source_idx = edge_ref.source();
                let old_target_idx = edge_ref.target();

                let new_source_idx = idx_map.get(&old_source_idx).unwrap().to_owned();
                let new_target_idx = idx_map.get(&old_target_idx).unwrap().to_owned();

                if !seen_edge.contains(&(new_source_idx, new_target_idx, edge.label.to_owned())) {
                    seen_edge.insert((new_source_idx, new_target_idx, edge.label.to_owned()));
                    merged_graph.add_edge(new_source_idx, new_target_idx, edge.to_owned());
                }
            }

            contributors.extend(model.contributors().iter().cloned());
        }

        Ok(GoCamModel {
            id: new_id.to_owned(),
            title: new_title.to_owned(),
            taxon,
            date: "now".to_owned(),
            contributors,
            graph: merged_graph,
        })
    }

    /// Remove a copy of the model with chemicals removed.  Where a
    /// chemical is between two activities, replace the chemical with
    /// a "provides input for" edge.
    pub fn remove_chemicals(&self) -> GoCamModel {
        let mut chemical_nodes = vec![];
        let mut new_model = self.clone();

        for (node_idx, node) in new_model.graph().node_references() {
            if node.node_type == GoCamNodeType::Chemical {
                chemical_nodes.push(node_idx);
            }
        }

        for chemical_node_idx in &chemical_nodes {
            let edges: Vec<_> = new_model.graph()
                .edges_directed(*chemical_node_idx, Direction::Incoming)
                .collect();

            let mut sources = vec![];
            let mut targets = vec![];

            for edge_ref in edges {
                let edge = edge_ref.weight();

                let edge_source_idx = edge_ref.source();

                if edge.id == "RO:0002234" {
                    // has output
                    sources.push(edge_source_idx);
                } else {
                    // has input
                    targets.push(edge_source_idx);
                }
            }


            for source_idx in &sources {
                let source_individual_gocam_id = {
                    new_model.graph.node_weight(*source_idx).unwrap().individual_gocam_id.clone()
                };

                for target_idx in &targets {
                    let existing_edges = new_model.graph.edges_connecting(*source_idx, *target_idx);

                    if existing_edges.into_iter().count() > 0 {
                        continue;
                    }

                    let target_node = new_model.graph.node_weight(*target_idx).unwrap();

                    let fact_gocam_id = format!("RO:0002413-{}-{}", source_individual_gocam_id,
                                                target_node.individual_gocam_id);

                    let edge_value = GoCamEdge {
                        fact_gocam_id,
                        id: "RO:0002413".to_owned(),
                        label: "provides input for".to_owned(),
                    };

                    new_model.graph.add_edge(*source_idx, *target_idx, edge_value);
                }
            }
        }

        new_model.graph.retain_nodes(|g, idx| {
            g.node_weight(idx).unwrap().node_type != GoCamNodeType::Chemical
        });

        new_model
    }
}

/// An overlap returned by [GoCamModel::find_overlaps()]
#[derive(Serialize, Deserialize, Clone, Debug)]
#[serde(rename_all = "snake_case")]
pub struct GoCamNodeOverlap {
    pub node_id: String,
    pub node_label: String,
    pub node_type: GoCamNodeType,
    pub has_input: Vec<GoCamInput>,
    pub has_output: Vec<GoCamOutput>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub occurs_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_process: Option<GoCamProcess>,
    pub overlapping_individual_ids: BTreeSet<IndividualId>,
    pub models: BTreeSet<(ModelId, String)>,
}

impl GoCamNodeOverlap {
    /// A unique ID for this overlap, created from the Individual IDs
    /// of the nodes/activities
    pub fn id(&self) -> String {
        self.overlapping_individual_ids.iter().cloned().collect::<Vec<_>>().join("-")
    }
}

/// An iterator over [GoCamNode], returned by [GoCamModel::node_iterator()]
pub struct NodeIterator<'a> {
    node_refs: NodeReferences<'a, GoCamNode>,
}

impl<'a> Iterator for NodeIterator<'a> {
    type Item = (NodeIndex, &'a GoCamNode);

    fn next(&mut self) -> Option<Self::Item> {
        self.node_refs.next().map(|(idx, node)| (idx, node))
    }
}

macro_rules! from_individual_type {
    ($type_name:ident, $doc:expr) => {

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[doc=$doc]
pub struct $type_name {
    pub id: String,
    pub label: String,
}

impl $type_name {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn label(&self) -> &str {
        &self.label
    }

    pub fn label_or_id(&self) -> &str {
        if self.label().len() > 0 {
            self.label()
        } else {
            self.id()
        }
    }
}

impl From<&IndividualType> for $type_name {
    fn from(individual_gene: &IndividualType) -> $type_name {
        $type_name {
            id: individual_gene.id().to_owned(),
            label: individual_gene.label().to_owned(),
        }
    }
}

impl Display for $type_name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label, self.id)?;
        Ok(())
    }
}

    };
}

from_individual_type!{GoCamGene, "A gene in a [GoCamNode], possibly enabling an activity"}

from_individual_type!{GoCamMRNA, "An mRNA - either an input/output or the enabler of an activity"}

from_individual_type!{GoCamChemical, "A chemical in a node, possibly enabling an activity"}

from_individual_type!{GoCamModifiedProtein, "A PRO modified protein in a node, possibly enabling an activity"}

from_individual_type!{GoCamInput, "The `has_input` of an activity"}

from_individual_type!{GoCamOutput, "The `has_output` of an activity"}

from_individual_type!{GoCamComplexComponent, "The GO component for the ComplexComponent variant of GoCamComponent"}

from_individual_type!{GoCamOtherComponent, "The GO component for the OtherComponent variant of GoCamComponent"}



/// A complex can have a GO complex ID (from the CC GO aspect) or a
/// Complex Portal ID
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct GoCamComplex {
    pub id: String,
    pub label: String,
    pub has_part_genes: Vec<GoCamGeneIdentifier>,
}

impl GoCamComplex {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn label(&self) -> &str {
        &self.label
    }
}

impl From<&IndividualType> for GoCamComplex {
    fn from(individual_complex: &IndividualType) -> GoCamComplex {
        GoCamComplex {
            id: individual_complex.id().to_owned(),
            label: individual_complex.label().to_owned(),
            has_part_genes: vec![],
        }
    }
}



/// A GO biological process
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, Hash)]
pub struct GoCamProcess {
    pub id: String,
    pub label: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_parent: Option<Box<GoCamProcess>>,
}

fn is_gene_id(identifier: &str) -> bool {
    ["PomBase:", "FB:", "UniProtKB:", "MGI:", "WB:", "RGD:", "RefSeq:",
     "Xenbase:", "SGD:", "ZFIN:", "RNAcentral:", "EMAPA:", "AGI_LocusCode:"]
        .iter().any(|s| identifier.starts_with(*s))
}

impl GoCamProcess {
    pub fn label_or_id(&self) -> &str {
        if self.label.len() > 0 {
            self.label.as_str()
        } else {
            self.id.as_str()
        }
    }
}

impl From<&IndividualType> for GoCamProcess {
    fn from(individual_process: &IndividualType) -> GoCamProcess {
        GoCamProcess {
            id: individual_process.id().to_owned(),
            label: individual_process.label().to_owned(),
            part_of_parent: None,
        }
    }
}

impl Display for GoCamProcess {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label, self.id)?;
        Ok(())
    }
}

/// A GO cellular component
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum GoCamComponent {
    ComplexComponent(GoCamComplexComponent),
    OtherComponent(GoCamOtherComponent),
}

impl GoCamComponent {
    /// The component ID, which will be a GO component term ID
    pub fn id(&self) -> &str {
        match self {
            GoCamComponent::ComplexComponent(complex) => &complex.id,
            GoCamComponent::OtherComponent(complex) => &complex.id,
        }
    }

    /// The component name, which will be a GO component term name
    pub fn label(&self) -> &str {
        match self {
            GoCamComponent::ComplexComponent(complex) => &complex.label,
            GoCamComponent::OtherComponent(complex) => &complex.label,
        }
    }

    /// The label (if set) otherwise the ID
    pub fn label_or_id(&self) -> &str {
        let label = self.label();

        if label.is_empty() {
            self.id()
        } else {
            self.label()
        }
    }
}

impl Display for GoCamComponent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "component: {} ({})", self.label(), self.id())?;
        Ok(())
    }
}

/// An enabler of an activity
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum GoCamEnabledBy {
    Complex(GoCamComplex),
    Gene(GoCamGene),
    Chemical(GoCamChemical),
    ModifiedProtein(GoCamModifiedProtein),
}

impl Display for GoCamEnabledBy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "enabled by: {} ({})", self.label(), self.id())?;
        Ok(())
    }
}

impl GoCamEnabledBy {
    /// The ID of the variant
    pub fn id(&self) -> &str {
        match self {
            GoCamEnabledBy::Complex(complex) => &complex.id,
            GoCamEnabledBy::Gene(gene) => &gene.id,
            GoCamEnabledBy::Chemical(chemical) => &chemical.id,
            GoCamEnabledBy::ModifiedProtein(modified_protein) => &modified_protein.id,
        }
    }

    /// The label of the variant
    pub fn label(&self) -> &str {
        match self {
            GoCamEnabledBy::Complex(complex) => &complex.label,
            GoCamEnabledBy::Gene(gene) => &gene.label,
            GoCamEnabledBy::Chemical(chemical) => &chemical.label,
            GoCamEnabledBy::ModifiedProtein(modified_protein) => &modified_protein.label,
        }
    }
}

/// The type of a node in a GoCamModel
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum GoCamNodeType {
    Unknown,
    Chemical,
    Gene(GoCamGene),
    MRNA(GoCamMRNA),
    ModifiedProtein(GoCamModifiedProtein),
    UnknownMRNA,
    Activity(GoCamEnabledBy),
}

impl GoCamNodeType {
    pub fn is_activity(&self) -> bool {
        if let GoCamNodeType::Activity(_) = self {
            return true;
        } else {
            return false;
        }
    }
}

impl Display for GoCamNodeType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GoCamNodeType::Unknown => write!(f, "unknown")?,
            GoCamNodeType::Chemical => write!(f, "chemical")?,
            GoCamNodeType::Gene(gene) => write!(f, "gene: {} ({})", gene.label, gene.id)?,
            GoCamNodeType::MRNA(mrna) => write!(f, "mrna: {} ({})", mrna.label, mrna.id)?,
            GoCamNodeType::ModifiedProtein(modified_protein) => {
                write!(f, "modified protein: {} ({})",
                       modified_protein.label, modified_protein.id)?;
            },
            GoCamNodeType::UnknownMRNA => write!(f, "unknown mRNA")?,
            GoCamNodeType::Activity(activity) => {
                write!(f, "{}", activity)?;
            },
        }
        Ok(())
    }
}

/// A gene, chemical, complex or modified protein OR an activity
/// (enabled by gene, chemical, complex or modified protein).  These
/// fields more or less match Figure 1 in
/// [the GO-CAM paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7012280/)
/// except for:
///  - `individual_gocam_id` which is the ID of the corresponding
///  Individual in the [GoCamRawModel]
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
pub struct GoCamNode {
    /// Individual ID in raw model from the JSON file
    pub individual_gocam_id: IndividualId,
    pub node_id: String,
    pub label: String,
    pub node_type: GoCamNodeType,
    pub has_input: BTreeSet<GoCamInput>,
    pub has_output: BTreeSet<GoCamOutput>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub occurs_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_process: Option<GoCamProcess>,
    pub source_ids: BTreeSet<IndividualId>,
    pub models: BTreeSet<(ModelId, ModelTitle)>,
}

impl GoCamNode {
    /// Return true iff this node is an activity
    pub fn is_activity(&self) -> bool {
        self.node_type.is_activity()
    }

    /// The type of this node, e.g. "chemical" or "enabled_by_gene"
    pub fn type_string(&self) -> &str {
        match &self.node_type {
            GoCamNodeType::Unknown => "unknown",
            GoCamNodeType::Chemical => "chemical",
            GoCamNodeType::UnknownMRNA => "unknown_mrna",
            GoCamNodeType::Gene(_) => "gene",
            GoCamNodeType::MRNA(_) => "mrna",
            GoCamNodeType::ModifiedProtein(_) => "modified_protein",
            GoCamNodeType::Activity(activity) => match activity {
                GoCamEnabledBy::Chemical(_) => "enabled_by_chemical",
                GoCamEnabledBy::Gene(_) => "enabled_by_gene",
                GoCamEnabledBy::ModifiedProtein(_) => "enabled_by_modified_protein",
                GoCamEnabledBy::Complex(_) => "enabled_by_complex",
            }
        }
    }

    /// Returns "X [enabled by] Y" for activities, otherwise returns
    /// the node label
    pub fn description(&self) -> String {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            format!("{} [enabled by] {}", self.label, enabler.label())
        } else {
            self.label.to_owned()
        }
    }

    /// The label of the enabler, otherwise ""
    pub fn enabler_label(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.label()
        } else {
            ""
        }
    }

    /// The label of the enabler, otherwise ""
    pub fn enabler_id(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.id()
        } else {
            ""
        }
    }

    /// If this node is an activity, return the ID of the enabler.
    /// Otherwise return the ID of the node (i.e. the chemical,
    /// complex, modified protein or gene ID).
    pub fn db_id(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.id()
        } else {
            &self.node_id
        }
    }
}

impl Display for GoCamNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {} {}", self.node_id,
               self.label, self.node_type)?;
        for has_input in &self.has_input {
            write!(f, " [has input] {}", has_input)?;
        }
        for has_output in &self.has_output {
            write!(f, " [has output] {}", has_output)?;
        }
        if let Some(ref located_in) = self.located_in {
            write!(f, " [located in] {}", located_in)?;
        }
        if let Some(ref occurs_in) = self.occurs_in {
            write!(f, " [occurs in] {}", occurs_in)?;
        }
        if let Some(ref part_of_process) = self.part_of_process {
            write!(f, " [part of] {}", part_of_process)?;
        }

        Ok(())
    }
}

/// An edge in the model - a causal relation between two activities.
///
///  - `id` - the term ID of the relation connecting two nodes, for
///  example "RO:0002629"
///  - `label` - the term name of the relation, e.g. "directly
///  positively regulates",
///
/// see [REL_NAMES] for a list of possible relations
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamEdge {
    pub fact_gocam_id: FactId,
    pub id: String,
    pub label: String,
}

impl Display for GoCamEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.label)?;
        Ok(())
    }
}

fn process_from_individual(individual: &Individual, model: &GoCamRawModel) -> GoCamProcess {
    let individual_type = individual.get_individual_type().unwrap();
    let mut process: GoCamProcess = individual_type.into();

    let subject_facts = model.facts_of_subject(&individual.id);

    for fact in subject_facts.iter() {
        if fact.property != "BFO:0000050" {
            continue;
        }
        let f_object = model.get_individual(&fact.object);
        let Some(f_object_type) = f_object.get_individual_type()
        else {
            continue;
        };
        if f_object_type.id().starts_with("GO:") {
           process.part_of_parent = Some(Box::new(f_object_type.into()));
        }
    }

    process
}

fn make_nodes(model: &GoCamRawModel) -> GoCamNodeMap {
    // genes, modified proteins and mRNAs that are the object
    // of a has_input or has_output relation
    let mut interesting_inputs_and_outputs = HashSet::new();

    for fact in model.facts() {
        let object = model.fact_object(fact);

        if !object.individual_is_gene() && !object.individual_is_modified_protein() &&
           !object.individual_is_unknown_mrna() && !object.individual_is_mrna() {
            continue;
        };

        if fact.property_label == "has input" || fact.property_label == "has output" {
            interesting_inputs_and_outputs.insert(fact.object.clone());
        }
    }

    let mut node_map = BTreeMap::new();

    for individual in model.individuals() {
        if individual.individual_is_activity(model) ||
            interesting_inputs_and_outputs.contains(&individual.id) ||
            individual.individual_is_chemical() &&
            !individual.individual_is_unknown_protein()
        {
            let model_id = model.id().to_owned();
            let model_title = model.title().to_owned();

            let Some(individual_type) = individual.get_individual_type()
            else {
                continue;
            };
            let detail =
                if individual.individual_is_chemical() {
                    GoCamNodeType::Chemical
                } else {
                    if interesting_inputs_and_outputs.contains(&individual.id) {
                        if individual.individual_is_unknown_mrna() {
                            GoCamNodeType::UnknownMRNA
                        } else {
                            if individual.individual_is_modified_protein() {
                                GoCamNodeType::ModifiedProtein(individual_type.into())
                            } else {
                                if individual.individual_is_mrna() {
                                    GoCamNodeType::MRNA(individual_type.into())
                                } else {
                                    GoCamNodeType::Gene(individual_type.into())
                                }
                            }
                        }
                    } else {
                        GoCamNodeType::Unknown
                    }
                };
            let mut source_ids = BTreeSet::new();
            source_ids.insert(individual.id.clone());
            let mut models = BTreeSet::new();
            models.insert((model_id, model_title));
            let gocam_node = GoCamNode {
                individual_gocam_id: individual.id.clone(),
                node_id: individual_type.id.clone().unwrap_or_else(|| "NO_ID".to_owned()),
                label: individual_type.label.clone().unwrap_or_else(|| "NO_LABEL".to_owned()),
                node_type: detail,
                has_input: BTreeSet::new(),
                has_output: BTreeSet::new(),
                located_in: None,
                occurs_in: None,
                part_of_process: None,
                source_ids,
                models,
            };

            node_map.insert(individual.id.clone(), gocam_node);
        }
    }

    let mut complex_map = HashMap::new();

    for individual in model.individuals() {
        if individual.individual_is_complex() {
            let Some(complex_type) = individual.types.get(0)
            else {
                continue;
            };

            let complex: GoCamComplex = complex_type.into();

            complex_map.insert(individual.id.clone(), complex);
        }
    }

    for fact in model.facts() {
        if fact.property_label == "has part" {
            let Some(complex) = complex_map.get_mut(&fact.subject)
            else {
                continue;
            };

            let object_individual = model.fact_object(fact);

            let Some(complex_part_type) = object_individual.types.get(0)
            else {
                continue;
            };

            let complex_gene = complex_part_type.id().to_owned();
            //            eprintln!("{}", complex_gene);
            complex.has_part_genes.push(complex_gene);
        }
    }

    for fact in model.facts() {
        let Some(subject_node) = node_map.get_mut(&fact.subject)
        else {
            continue;
        };

        let object_individual = model.fact_object(fact);
        let Some(object_type) = object_individual.get_individual_type()
        else {
            continue;
        };

        match fact.property_label.as_str() {
            "enabled by" => {
                if let Some(ref object_type_id) = object_type.id {
                    if is_gene_id(object_type_id) {
                        let gene_enabler = GoCamEnabledBy::Gene(object_type.into());
                        subject_node.node_type = GoCamNodeType::Activity(gene_enabler);
                    }
                    else if object_type_id.starts_with("CHEBI:") {
                        let chemical_enabler = GoCamEnabledBy::Chemical(object_type.into());
                        subject_node.node_type = GoCamNodeType::Activity(chemical_enabler);
                    }
                    else if object_type_id.starts_with("GO:") || object_type_id.starts_with("ComplexPortal:") {
                        let complex = complex_map.get(&fact.object)
                            .expect(&format!("expected complex {}", fact.object))
                            .to_owned();
                        let complex_enabler = GoCamEnabledBy::Complex(complex);
                        subject_node.node_type = GoCamNodeType::Activity(complex_enabler);
                    }
                    else if object_type_id.starts_with("PR:") {
                        let modified_protein_enabler = GoCamEnabledBy::ModifiedProtein(object_type.into());
                        subject_node.node_type = GoCamNodeType::Activity(modified_protein_enabler);
                    }
                    else  {
                        eprintln!("can't handle enabled by object: {} - {}", object_type_id, object_individual.id);
                    }
                }
            },
            "has input" => {
                subject_node.has_input.insert(object_type.into());
            },
            "has output" => {
                subject_node.has_output.insert(object_type.into());
            },
            "located in" => {
                if subject_node.located_in.is_some() {
                    panic!("{}: {} is located in multiple components", model.id(),
                           subject_node.description());
                }
                let located_in =
                    if object_individual.individual_is_complex() {
                        GoCamComponent::ComplexComponent(object_type.into())
                    } else {
                        GoCamComponent::OtherComponent(object_type.into())
                    };
                subject_node.located_in = Some(located_in);
            },
            "occurs in" => {
                if subject_node.occurs_in.is_some() {
                    eprintln!("{}: {} occurs in multiple components", model.id(),
                              subject_node.description());
                }
                let occurs_in =
                    if object_individual.individual_is_complex() {
                        GoCamComponent::ComplexComponent(object_type.into())
                    } else {
                        GoCamComponent::OtherComponent(object_type.into())
                    };
                subject_node.occurs_in = Some(occurs_in);
            },
            "part of" => {
                let process = process_from_individual(object_individual, model);

                subject_node.part_of_process = Some(process);
            },
            &_ => {
                // eprintln!("ignoring rel from fact: {} {}", fact.property_label, fact.id());
            }
        }
    }

    node_map
}

/// Read from a JSON source and return a [GoCamModel].
///
/// ## Example
/// ```
/// let mut source = std::fs::File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
/// let model = pombase_gocam::parse_gocam_model(&mut source).unwrap();
/// ```
pub fn parse_gocam_model(source: &mut dyn Read) -> Result<GoCamModel> {
    let raw_model = gocam_parse_raw(source)?;

    let model = GoCamModel::new(raw_model);

    Ok(model)
}

fn make_graph(model: &GoCamRawModel) -> GoCamGraph {
    let mut graph = GoCamGraph::new();

    let node_map = make_nodes(model);

    let mut id_map = HashMap::new();

    for node in node_map.values() {
        let idx = graph.add_node(node.to_owned());
        id_map.insert(node.individual_gocam_id.clone(), idx);
    }

    for fact in model.facts() {
        let subject_id = &fact.subject;
        let object_id = &fact.object;

        if let (Some(__), Some(_)) =
            (node_map.get(subject_id), node_map.get(object_id))
        {
            let subject_idx = id_map.get(subject_id).unwrap();
            let object_idx = id_map.get(object_id).unwrap();

            let edge = GoCamEdge {
                fact_gocam_id: fact.id(),
                id: fact.property.clone(),
                label: fact.property_label.clone(),
            };

            graph.add_edge(*subject_idx, *object_idx, edge);
        }
    }

    graph
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn parse_test() {
        let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:66187e4700001744");

        assert_eq!(model.title(), "meiotic cohesion protection in anaphase I (GO:1990813)");
        assert_eq!(model.taxon(), "NCBITaxon:4896");

        let (first_node_idx, first_node) = model.node_iterator().next().unwrap();

        assert_eq!(first_node_idx, 0.into());
        assert_eq!(first_node.node_id, "GO:0140483");
        let (first_node_model_id, first_node_model_title) = first_node.models.iter().next().unwrap();
        assert_eq!(first_node_model_id, "gomodel:66187e4700001744");
        assert_eq!(first_node_model_title, "meiotic cohesion protection in anaphase I (GO:1990813)");

        assert_eq!(model.node_iterator().count(), 12);

        assert_eq!(first_node.to_string(),
                   "GO:0140483 kinetochore adaptor activity enabled by: moa1 Spom (PomBase:SPAC15E1.07c) [occurs in] component: kinetochore (GO:0000776) [part of] meiotic centromeric cohesion protection in anaphase I (GO:1990813)");

    }

    #[test]
    fn overlap_test() {
        let mut source1 = File::open("tests/data/gomodel_665912ed00000192.json").unwrap();
        let model1 = parse_gocam_model(&mut source1).unwrap();
        assert_eq!(model1.id(), "gomodel:665912ed00000192");
        assert_eq!(model1.node_iterator().count(), 34);

        let mut source2 = File::open("tests/data/gomodel_663d668500002178.json").unwrap();
        let model2 = parse_gocam_model(&mut source2).unwrap();
        assert_eq!(model2.id(), "gomodel:663d668500002178");
        assert_eq!(model2.node_iterator().count(), 14);

        let overlaps = GoCamModel::find_overlaps(&[model1, model2]);

        assert_eq!(overlaps.len(), 2);

        let mut expected_node_ids = HashSet::new();
        expected_node_ids.insert("CHEBI:16749".to_owned());
        expected_node_ids.insert("GO:0003881".to_owned());

        let activity_node_ids: HashSet<_> = overlaps.iter().map(|o| o.node_id.clone()).collect();
        assert_eq!(activity_node_ids, expected_node_ids);

        let mut expected_enabled_by_ids = HashSet::new();
        expected_enabled_by_ids.insert("PomBase:SPAC1D4.08".to_owned());
        let mut expected_chemical_ids = HashSet::new();
        expected_chemical_ids.insert("CHEBI:16749".to_owned());

        let mut activity_enabled_by_ids = HashSet::new();
        let mut chemical_ids = HashSet::new();

        for overlap in overlaps {
            match overlap.node_type {
                GoCamNodeType::Activity(ref enabled_by) => {
                    activity_enabled_by_ids.insert(enabled_by.id().to_owned());
                },
                GoCamNodeType::Chemical => {
                    chemical_ids.insert(overlap.node_id);
                },
                _ => panic!(),
            }
        }

        assert_eq!(activity_enabled_by_ids, expected_enabled_by_ids);
        assert_eq!(chemical_ids, expected_chemical_ids);
    }

    #[test]
    fn merge_test() {
        let mut source1 = File::open("tests/data/gomodel_663d668500002178.json").unwrap();
        let model1 = parse_gocam_model(&mut source1).unwrap();

        let mut source2 = File::open("tests/data/gomodel_665912ed00000192.json").unwrap();
        let model2 = parse_gocam_model(&mut source2).unwrap();

        let merged = GoCamModel::merge_models("new_id", "new_title",
                                              &[model1, model2]).unwrap();

        assert_eq!(merged.node_iterator().count(), 46);

        let merged_ids: HashSet<_> =
            merged.node_iterator().filter_map(|(_, node)| if node.models.len() >= 2 {
                Some((node.label.clone(), node.db_id().to_owned()))
            } else {
                None
            })
            .collect();

        let mut expected_ids = HashSet::new();
        expected_ids.insert(("1-phosphatidyl-1D-myo-inositol".to_owned(),
                             "CHEBI:16749".to_owned()));
        expected_ids.insert(("CDP-diacylglycerol-inositol 3-phosphatidyltransferase activity".to_owned(),
                                         "PomBase:SPAC1D4.08".to_owned()));

        assert_eq!(merged_ids, expected_ids);
    }

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

        let overlaps = GoCamModel::find_overlaps(&[model1, model2, model3]);

        assert_eq!(overlaps.len(), 2);

        let overlap_chemical = &overlaps[0];
        assert_eq!(overlap_chemical.node_label, "O-acetyl-L-homoserine");

        let overlap_activity = &overlaps[1];

        assert_eq!(overlap_activity.node_label, "homoserine O-acetyltransferase activity");
        assert_eq!(overlap_activity.models.len(), 2);

        assert_eq!(overlap_activity.part_of_process.as_ref().unwrap().id, "GO:0071266");
        assert_eq!(overlap_activity.part_of_process.as_ref().unwrap().label,
                   "'de novo' L-methionine biosynthetic process");
        assert_eq!(overlap_activity.occurs_in.as_ref().unwrap().id(), "GO:0005829");
        assert_eq!(overlap_activity.occurs_in.as_ref().unwrap().label(), "cytosol");

        let first_overlapping_individual =
            overlap_activity.overlapping_individual_ids.iter().next().unwrap();
        assert_eq!(first_overlapping_individual,
                   "gomodel:66a3e0bb00001342/678073a900003752");
    }

    #[test]
    fn remove_chemicals_test() {
        let mut source1 = File::open("tests/data/gomodel_66a3e0bb00001342.json").unwrap();
        let model = parse_gocam_model(&mut source1).unwrap();

        let new_model = model.remove_chemicals();

        assert_eq!(new_model.node_count(), 18);
    }
}
