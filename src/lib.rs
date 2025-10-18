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
//! use pombase_gocam::{GoCamModel, GoCamNodeType, GoCamActivity};
//! let model = GoCamModel::new(raw_model);
//!
//! for (_, node) in model.node_iterator() {
//!     println!("node: {}", node);
//!     if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = node.node_type {
//!         println!("enabler ID: {}", enabler.id());
//!     }
//! }
//!```

use std::{cmp::Ordering, collections::{BTreeMap, BTreeSet, HashMap, HashSet}, fmt::{self, Display}, io::Read, sync::LazyLock};
use std::hash::{Hash, Hasher};

extern crate serde_json;
extern crate serde_yaml;
#[macro_use] extern crate serde_derive;

pub mod raw;
pub mod graph;
pub mod overlaps;
pub mod gocam_py;

use raw::{gocam_parse_raw, FactId, GoCamRawModel, Individual, IndividualId, IndividualType};

use phf::phf_map;

use anyhow::{Result, anyhow};

use petgraph::{graph::{EdgeReference, NodeIndex, NodeReferences},
               visit::{Bfs, EdgeRef, IntoNodeReferences, UndirectedAdaptor},
               Direction, Graph};
use regex::Regex;

use crate::overlaps::find_overlaps;

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

pub type GoCamModelId = String;
pub type GoCamModelTitle = String;

/// The Graph representation of a GO-CAM model.  See: [GoCamModel::graph()]
pub type GoCamGraph = Graph::<GoCamNode, GoCamEdge>;

/// A gene ID with DB prefix, like "PomBase:SPAC9E9.05"
pub type GoCamGeneIdentifier = String;

/// A gene name like "cdc2"
pub type GoCamGeneName = String;

type GoCamNodeMap = BTreeMap<IndividualId, GoCamNode>;

static TITLE_GO_TERM_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\((\s*GO:\d+\s*)\)").unwrap());

/// A high level representation of the model with nodes for
/// activities, chemicals, complexes etc. and edges for causal
/// dependencies between nodes/activities.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamModel {
    id: String,
    title: String,
    title_process_term_ids: HashSet<String>,
    taxon: String,
    date: String,
    contributors: BTreeSet<String>,
    graph: GoCamGraph,
    gene_name_map: HashMap<String, String>,
    pro_term_to_gene_map: HashMap<String, String>,
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

#[derive(PartialEq, Eq, Hash)]
pub enum RemoveType {
    Chemicals,
    Targets,
}

fn process_term_ids_from_title(title: &str) -> HashSet<String> {
    TITLE_GO_TERM_RE.captures_iter(title).map(|c| c[1].to_owned()).collect()
}

impl GoCamModel {
    /// Create a [GoCamModel] from a [GoCamRawModel]
    /// `gene_name_map` is a map from gene indentifiers to gene names (used
    /// to fill gene names of complex members)
    /// `pro_term_to_gene_map` is a map from PRO ID to gene identifier (needed
    /// because this information isn't in the raw JSON file)
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

        let title_process_term_ids = process_term_ids_from_title(&raw_model.title());

        let model = GoCamModel {
            id: raw_model.id().to_owned(),
            title: raw_model.title().to_owned(),
            title_process_term_ids,
            taxon: raw_model.taxon().to_owned(),
            date: raw_model.date().to_owned(),
            contributors: raw_model.contributors(),
            graph,
            gene_name_map: HashMap::new(),
            pro_term_to_gene_map: HashMap::new(),
        };

        model
    }

    /// Add a map from gene ID to gene name to the model. This map is used to
    /// make the return value of genes_enabling_activities() more useful.
    pub fn add_gene_name_map(&mut self, gene_name_map: &HashMap<String, String>) {
        let mut enabler_gene_name_map = HashMap::new();

        for (orig_gene_id, _) in self.genes_enabling_activities() {
            let gene_id =
                if orig_gene_id.contains(":") {
                    orig_gene_id.split(":").last().unwrap().to_owned()
                } else {
                    orig_gene_id.to_owned()
                };

            if let Some(name) = gene_name_map.get(&gene_id) {
                enabler_gene_name_map.insert(orig_gene_id, name.to_owned());
                enabler_gene_name_map.insert(gene_id.clone(), name.to_owned());
            }
        }

        self.gene_name_map = enabler_gene_name_map;
    }

    /// Add a map from PRO ID to gene ID to the model.  This is needed so that
    /// we can accurately count and list enabling genes in models where
    /// activities are enabled by modified proteins.
    pub fn add_pro_term_to_gene_map(&mut self, pro_term_to_gene_map: &HashMap<String, String>) {
        self.pro_term_to_gene_map = pro_term_to_gene_map.to_owned();
    }

    /// Return the [petgraph::Graph] representation of the model
    pub fn graph(&self) -> &GoCamGraph {
        &self.graph
    }

    /// Return an iterator over the nodes ([GoCamNode]) of the model
    pub fn node_iterator(&self) -> NodeIterator<'_> {
        NodeIterator {
            node_refs: self.graph().node_references(),
        }
    }

    /// Number of nodes in the graph
    pub fn node_count(&self) -> usize {
        self.graph().node_count()
    }

    /// Return the IDs of all the genes in the model, including input and output genes, genes
    /// of mRNAs and genes in complexes.
    pub fn genes_in_model(&self)
         -> BTreeSet<GoCamGeneIdentifier>
    {
        let mut ret_genes = BTreeSet::new();

        for (_, node) in self.node_iterator() {
            match &node.node_type {
                GoCamNodeType::Gene(gene) => {
                    ret_genes.insert(gene.id().to_owned());
                },
                GoCamNodeType::MRNA(mrna) => {
                    // temporary hack
                    if let Some((prefix, suffix)) = mrna.id().rsplit_once('.') {
                        if suffix.len() == 1 && suffix.chars().next().unwrap().is_numeric() {
                            ret_genes.insert(prefix.to_owned());
                        }
                    }
                },
                GoCamNodeType::ModifiedProtein(modified_protein) => {
                    let pro_id = modified_protein.id();
                    if let Some(gene) = self.pro_term_to_gene_map.get(pro_id) {
                        ret_genes.insert(gene.to_owned());
                    }
                },
                GoCamNodeType::Activity(GoCamActivity { enabler: enabled_by, .. }) => {
                    match enabled_by {
                        GoCamEnabledBy::Gene(gene) => {
                            ret_genes.insert(gene.id().to_owned());
                        },
                        GoCamEnabledBy::Complex(complex) => {
                            for gene in &complex.has_part_genes {
                                ret_genes.insert(gene.to_owned());
                            }
                        },
                        GoCamEnabledBy::ModifiedProtein(modified_protein) => {
                            let pro_id = modified_protein.id();
                            if let Some(gene) = self.pro_term_to_gene_map.get(pro_id) {
                                ret_genes.insert(gene.to_owned());
                            }
                        },
                        _ => (),
                    }
                },
                _ => (),
            }
        }

        ret_genes
    }

    /// Return the IDs and name of the genes that enable an activity, including genes in complexes.
    /// The pro_term_to_gene_map is a map of protein ontology term IDs to gene IDs that
    /// allows counting modified genes.
    pub fn genes_enabling_activities(&self)
          -> HashMap<GoCamGeneIdentifier, Option<GoCamGeneName>>
    {
        let mut ret_genes = HashMap::new();

        let get_gene_and_name = |gene_id: &GoCamGeneIdentifier| {
            let gene_name =
                if gene_id.contains(":") {
                    let id = gene_id.split(":").last().unwrap();
                    self.gene_name_map.get(id)
                } else {
                    self.gene_name_map.get(gene_id)
                };
            (gene_id.to_owned(), gene_name.map(String::to_owned))
        };

        for (_, node) in self.node_iterator() {
            match &node.node_type {
                GoCamNodeType::Activity(GoCamActivity { enabler: enabled_by, .. }) => {
                    match enabled_by {
                        GoCamEnabledBy::Gene(gene) => {
                            let (gene_id,_name) = get_gene_and_name(&gene.id);
                            ret_genes.insert(gene_id,_name);
                        },
                        GoCamEnabledBy::Complex(complex) => {
                            for gene_id in &complex.has_part_genes {
                                let (gene_id,_name) = get_gene_and_name(gene_id);
                                ret_genes.insert(gene_id,_name);
                            }
                        },
                        GoCamEnabledBy::ModifiedProtein(modified_protein) => {
                            let pro_id = modified_protein.id();
                            if let Some(gene_id) = self.pro_term_to_gene_map.get(pro_id) {
                                let (gene_id,_name) = get_gene_and_name(gene_id);
                                ret_genes.insert(gene_id,_name);
                            }
                        },
                        _ => (),
                    }
                },
                _ => (),
            }
        }

        ret_genes
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

    /// Merge the `models` that have nodes in common, returning a new
    /// [GoCamModel] with the `new_id` as the ID and `new_title` as
    /// the title.
    ///
    /// We use the result of calling [find_overlaps()] to find
    /// nodes in common between all the `models`.
    pub fn merge_models(new_id: &str, new_title: &str, models: &[GoCamModel])
        -> Result<GoCamModel>
    {
        let mut merged_graph = GoCamGraph::new();

        let overlaps = find_overlaps(models);

        let mut overlap_map = HashMap::new();

        for overlap in overlaps.into_iter() {
            let overlap_id = overlap.id();

            let models = overlap.models.iter()
                .map(|(model_id, model_title, _)| (model_id.to_owned(), model_title.to_owned()))
                .collect();

            let overlap_node = GoCamNode {
                individual_gocam_id: overlap_id,
                node_id: overlap.node_id,
                label: overlap.node_label,
                node_type: overlap.node_type,
                occurs_in: overlap.occurs_in,
                part_of_process: overlap.part_of_process,
                happens_during: overlap.happens_during,
                source_ids: overlap.overlapping_individual_ids.clone(),
                original_model_id: overlap.original_model_id,
                models,
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

        let mut gene_name_map = HashMap::new();

        for model in models {
            gene_name_map.extend(model.gene_name_map.clone().into_iter());
        }

        let first_model = models.first().unwrap();

        Ok(GoCamModel {
            id: new_id.to_owned(),
            title: new_title.to_owned(),
            title_process_term_ids: HashSet::new(),
            taxon,
            date: "now".to_owned(),
            contributors,
            graph: merged_graph,
            gene_name_map,
            pro_term_to_gene_map: first_model.pro_term_to_gene_map.clone(),
        })
    }

    /// Return a copy of the model with chemicals or all inputs and outputs removed.
    /// Where a chemical is between two activities, replace the chemical with
    /// a "provides input for" edge.
    pub fn remove_nodes(&self, remove_types: HashSet<RemoveType>) -> GoCamModel {
        let mut nodes_to_remove = BTreeSet::new();
        let mut new_model = self.clone();

        for (node_idx, node) in new_model.graph().node_references() {
            if remove_types.contains(&RemoveType::Chemicals) &&
                node.node_type.is_chemical() {
                    nodes_to_remove.insert(node_idx);
            }
            if remove_types.contains(&RemoveType::Targets) &&
                !node.is_activity() && !node.node_type.is_chemical() {
                nodes_to_remove.insert(node_idx);
            }
        }

        let nodes_to_remove: BTreeSet<_> = nodes_to_remove.into_iter()
            .filter(|node_idx| {
                let in_iter = new_model.graph().edges_directed(*node_idx, Direction::Incoming);
                let out_iter = new_model.graph().edges_directed(*node_idx, Direction::Outgoing);
                let edge_iter = in_iter.chain(out_iter);

                for in_edge_ref in edge_iter {
                    let in_edge = in_edge_ref.weight();
                    if in_edge.id != "RO:0002233" && in_edge.id != "RO:0002234" {
                        return false;
                    }
                }

                true
            })
            .collect();

        for remove_node_idx in &nodes_to_remove {
            let edges: Vec<_> = new_model.graph()
                .edges_directed(*remove_node_idx, Direction::Incoming)
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

        new_model.graph.retain_nodes(|_, idx| {
            !nodes_to_remove.contains(&idx)
        });

        new_model
    }

    // If the activity at `activity_idx` has only incoming relations (excluding inputs/outputs),
    // return Incoming, unless all the incoming relations are "constitutively upstream" in which
    // case return IncomingConstitutivelyUpstream.
    // If the activity has only outgoing relations return Outgoing.
    // Otherwise return None.
    fn node_rel_direction(&self, activity_idx: NodeIndex) -> GoCamDirection {
        let is_causal_edge = |edge_ref:  EdgeReference<GoCamEdge>| {
            let edge = edge_ref.weight();
            edge.id != "RO:0002234" && edge.id != "RO:0002233"
        };

        let is_constitutively_upstream_edge = |edge_ref:  EdgeReference<GoCamEdge>| {
            let edge = edge_ref.weight();
            edge.id == "RO:0012009"
        };

        let graph = self.graph();

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
            let node = self.graph().node_weight(activity_idx).unwrap();
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

    /// Return a new GoCamModel after removing all but the largest connected subgraph
    pub fn retain_largest_subgraph(&self) -> GoCamModel {
        let mut subgraphs: BTreeMap<NodeIndex, BTreeSet<NodeIndex>> = BTreeMap::new();

        let graph = UndirectedAdaptor(self.graph());

        for (node_idx, _) in self.node_iterator() {
            for subgraph_node_indexes in subgraphs.values() {
                if subgraph_node_indexes.contains(&node_idx) {
                    continue;
                }
            }

            let mut subgraph_nodes = BTreeSet::new();

            let mut bfs = Bfs::new(graph, node_idx);

            while let Some(nx) = bfs.next(&graph) {
                subgraph_nodes.insert(nx);
            }

            subgraphs.insert(node_idx, subgraph_nodes);
        }

        let (mut current_max_len, mut current_max_node_idx) = {
            let first_entry = subgraphs.first_entry().unwrap();
            (first_entry.get().len(), first_entry.key().to_owned())
        };

        for (node_idx, subgraph_nodes) in subgraphs.iter() {
            if subgraph_nodes.len() > current_max_len {
                current_max_len = subgraph_nodes.len();
                current_max_node_idx = *node_idx;
            }
        }

        let max_subgraph_nodes = subgraphs.get(&current_max_node_idx).unwrap();

        let mut graph = self.graph.clone();

        graph.retain_nodes(|_, node_idx| {
            max_subgraph_nodes.contains(&node_idx)
        });

        let title_process_term_ids = process_term_ids_from_title(&self.title);

        GoCamModel {
            id: self.id.clone(),
            title: self.title.clone(),
            title_process_term_ids,
            taxon: self.taxon.clone(),
            contributors: self.contributors.clone(),
            date: self.date.clone(),
            graph,
            gene_name_map: self.gene_name_map.clone(),
            pro_term_to_gene_map: self.pro_term_to_gene_map.clone(),
        }
    }

    /// Given a model with the inputs/outputs removed (using remove_nodes()), return a clone
    /// containing only those activities enabled by a gene in the `retain_genes` set.
    pub fn retain_enabling_genes(&self, retain_genes: &BTreeSet<GoCamGeneIdentifier>)
       -> GoCamModel
    {
        let mut graph = self.graph.clone();

        graph.retain_nodes(|_, node_idx| {
            let node = self.graph.node_weight(node_idx).unwrap();
            if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = node.node_type {
                if let GoCamEnabledBy::Gene(gene) = enabler {
                    let split = gene.id().split(":").last().unwrap();
                    if retain_genes.contains(split) {
                       return true;
                   }
                }
                if let GoCamEnabledBy::Complex(complex) = enabler {
                    for gene in &complex.has_part_genes {
                        let split = gene.split(":").last().unwrap();
                        if retain_genes.contains(split) {
                           return true;
                        }
                    }
                }
            }
            false
        });

        GoCamModel {
            id: self.id.clone(),
            title: self.title.clone(),
            title_process_term_ids: self.title_process_term_ids.clone(),
            taxon: self.taxon.clone(),
            contributors: self.contributors.clone(),
            date: self.date.clone(),
            graph,
            gene_name_map: self.gene_name_map.clone(),
            pro_term_to_gene_map: self.pro_term_to_gene_map.clone(),
        }
    }

    /// Return true iff any of the genes in the set enable an activity in this model,
    /// or are part of a complex that enables an activity in the model
    pub fn model_activity_enabled_by(&self, genes: &BTreeSet<GoCamGeneIdentifier>)
        -> bool
    {
        for (_, node) in self.node_iterator() {
            if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = node.node_type {
                if let GoCamEnabledBy::Gene(gene) = enabler {
                    let split = gene.id().split(":").last().unwrap();
                    if genes.contains(split) {
                       return true;
                   }
                }
                if let GoCamEnabledBy::Complex(complex) = enabler {
                    for gene in &complex.has_part_genes {
                        let split = gene.split(":").last().unwrap();
                        if genes.contains(split) {
                           return true;
                        }
                    }
                }
            }
        }

        false
    }
}

pub type GoCamModelIdTitle = (GoCamModelId, GoCamModelTitle, GoCamDirection);

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum GoCamDirection {
    Incoming,
    IncomingConstitutivelyUpstream,
    Outgoing,
    None,
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

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GoCamActivity {
    pub enabler: GoCamEnabledBy,
    pub inputs: BTreeSet<GoCamInput>,
    pub outputs: BTreeSet<GoCamOutput>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[doc="A gene in a [GoCamNode], possibly enabling an activity"]
pub struct GoCamGene {
    pub id: String,
    pub label: String,
    pub part_of_complex: Option<GoCamComplex>,
}

impl GoCamGene {
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn label(&self) -> String {
        if let Some(ref complex) = self.part_of_complex {
            format!("{} part of complex {}", self.label, complex.label())
        } else {
            self.label.to_owned()
        }
    }

    pub fn label_or_id(&self) -> String {
        if self.label().len() > 0 {
            self.label()
        } else {
            self.id().to_owned()
        }
    }
}

impl GoCamGene {
    fn new(individual_gene: &IndividualType, part_of_complex: &Option<GoCamComplex>) -> GoCamGene {
        GoCamGene {
            id: individual_gene.id().to_owned(),
            label: individual_gene.label().to_owned(),
            part_of_complex: part_of_complex.to_owned(),
        }
    }
}

impl Display for GoCamGene {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label, self.id)?;
        Ok(())
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[doc = "The `has_input` of an activity"]
pub struct GoCamInput {
    pub id: String,
    pub label: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub occurs_in: BTreeSet<GoCamComponent>,
}
impl GoCamInput {
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
impl From<&IndividualType> for GoCamInput {
    fn from(individual_gene: &IndividualType) -> GoCamInput {
        GoCamInput {
            id: individual_gene.id().to_owned(),
            label: individual_gene.label().to_owned(),
            located_in: None,
            occurs_in: BTreeSet::new(),
        }
    }
}
impl Display for GoCamInput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label, self.id)?;
        if let Some(ref located_in) = self.located_in {
            write!(f, " located in {}", located_in)?;
        }
        if !self.occurs_in.is_empty() {
            let occurs_in = self.occurs_in.iter()
                .map(|o| o.to_string())
                .collect::<Vec<_>>()
                .join(",");
            write!(f, " occurs in {}", occurs_in)?;
        }
        Ok(())
    }
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[doc = "The `has_output` of an activity"]
pub struct GoCamOutput {
    pub id: String,
    pub label: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub occurs_in: BTreeSet<GoCamComponent>,
}

impl GoCamOutput {
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
impl From<&IndividualType> for GoCamOutput {
    fn from(individual_gene: &IndividualType) -> GoCamOutput {
        GoCamOutput {
            id: individual_gene.id().to_owned(),
            label: individual_gene.label().to_owned(),
            located_in: None,
            occurs_in: BTreeSet::new(),
        }
    }
}
impl Display for GoCamOutput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label, self.id)?;
        if let Some(ref located_in) = self.located_in {
            write!(f, " located in {}", located_in)?;
        }
        if !self.occurs_in.is_empty() {
            let occurs_in = self.occurs_in.iter()
                .map(|o| o.to_string())
                .collect::<Vec<_>>()
                .join(",");
            write!(f, " occurs in {}", occurs_in)?;
        }
        Ok(())
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[doc = "A chemical which will be the input, output or enabler of an activity "]
pub struct GoCamChemical {
    pub id: String,
    pub label: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<GoCamComponent>,
}
impl GoCamChemical {
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

from_individual_type!{GoCamMRNA, "An mRNA - either an input/output or the enabler of an activity"}

from_individual_type!{GoCamModifiedProtein, "A PRO modified protein in a node, possibly enabling an activity"}

from_individual_type!{GoCamComplexComponent, "The component for the ComplexComponent variant of GoCamComponent"}

from_individual_type!{GoCamOtherComponent, "The component for the OtherComponent variant of GoCamComponent"}



/// A complex can have a GO complex ID (from the CC GO aspect) or a
/// Complex Portal ID
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamComplex {
    pub id: String,
    pub label: String,
    pub has_part_genes: BTreeSet<GoCamGeneIdentifier>,
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
            has_part_genes: BTreeSet::new(),
        }
    }
}

impl PartialEq for GoCamComplex {
    fn eq(&self, other: &GoCamComplex) -> bool {
        self.id == other.id
    }
}
impl Eq for GoCamComplex { }

impl Ord for GoCamComplex {
    fn cmp(&self, other: &GoCamComplex) -> Ordering {
        self.id.cmp(&other.id)
    }
}
impl PartialOrd for GoCamComplex {
    fn partial_cmp(&self, other: &GoCamComplex) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


impl Hash for GoCamComplex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// A GO biological process
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GoCamProcess {
    pub id: String,
    pub label: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_parent: Option<Box<GoCamProcess>>,
}

fn is_gene_id(identifier: &str) -> bool {
    ["PomBase:", "FB:", "UniProtKB:", "MGI:", "WB:", "RGD:", "RefSeq:",
     "Xenbase:", "SGD:", "ZFIN:", "RNAcentral:", "EMAPA:", "AGI_LocusCode:",
     "EcoCyc:"]
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

/// A complex or a GO cellular component
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum GoCamComponent {
    /// A complex from Complex Portal
    ComplexComponent(GoCamComplexComponent),
    /// a GO component
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
        write!(f, "{} ({})", self.label(), self.id())?;
        Ok(())
    }
}

/// An enabler of an activity
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
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

impl From<&IndividualType> for GoCamChemical {
    fn from(individual_process: &IndividualType) -> GoCamChemical {
        GoCamChemical {
            id: individual_process.id().to_owned(),
            label: individual_process.label().to_owned(),
            located_in: None,
        }
    }
}

/// The type of a node in a GoCamModel
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
pub enum GoCamNodeType {
    Unknown,
    Chemical(GoCamChemical),
    Gene(GoCamGene),
    MRNA(GoCamMRNA),
    ModifiedProtein(GoCamModifiedProtein),
    UnknownMRNA,
    Activity(GoCamActivity),
}

impl GoCamNodeType {
    pub fn is_activity(&self) -> bool {
        if let GoCamNodeType::Activity { .. } = self {
            return true;
        } else {
            return false;
        }
    }

    pub fn is_chemical(&self) -> bool {
        if let GoCamNodeType::Chemical(_) = self {
            return true;
        } else {
            return false;
        }
    }

    pub fn id(&self) -> &str {
        match self {
            GoCamNodeType::Unknown => "unknown",
            GoCamNodeType::Chemical(chemical) => &chemical.id,
            GoCamNodeType::Gene(gene) => &gene.id,
            GoCamNodeType::MRNA(mrna) => &mrna.id,
            GoCamNodeType::ModifiedProtein(modified_protein) => &modified_protein.id,
            GoCamNodeType::UnknownMRNA => "unknown mRNA",
            GoCamNodeType::Activity(GoCamActivity { enabler, inputs: _, outputs: _ }) => enabler.id(),
        }
    }

    pub fn label(&self) -> &str {
        match self {
            GoCamNodeType::Unknown => "unknown",
            GoCamNodeType::Chemical(chemical) => &chemical.label,
            GoCamNodeType::Gene(gene) => &gene.label,
            GoCamNodeType::MRNA(mrna) => &mrna.label,
            GoCamNodeType::ModifiedProtein(modified_protein) => &modified_protein.label,
            GoCamNodeType::UnknownMRNA => "unknown mRNA",
            GoCamNodeType::Activity(GoCamActivity { enabler, inputs: _, outputs: _ }) => enabler.label(),
        }
    }
}

impl Display for GoCamNodeType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GoCamNodeType::Unknown => write!(f, "unknown")?,
            GoCamNodeType::Chemical(chemical) => {
                write!(f, "chemical: {} ({})", chemical.label, chemical.id)?
            },
            GoCamNodeType::Gene(gene) => write!(f, "gene: {} ({})", gene.label, gene.id)?,
            GoCamNodeType::MRNA(mrna) => write!(f, "mrna: {} ({})", mrna.label, mrna.id)?,
            GoCamNodeType::ModifiedProtein(modified_protein) => {
                write!(f, "modified protein: {} ({})",
                       modified_protein.label, modified_protein.id)?;
            },
            GoCamNodeType::UnknownMRNA => write!(f, "unknown mRNA")?,
            GoCamNodeType::Activity(GoCamActivity { enabler, inputs: _, outputs: _ }) => {
                write!(f, "{}", enabler)?;
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

    #[serde(skip_serializing_if="BTreeSet::is_empty", default)]
    pub occurs_in: BTreeSet<GoCamComponent>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of_process: Option<GoCamProcess>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub happens_during: Option<GoCamProcess>,
    pub source_ids: BTreeSet<IndividualId>,
    pub original_model_id: Option<GoCamModelId>,
    pub models: BTreeSet<(GoCamModelId, GoCamModelTitle)>,
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
            GoCamNodeType::Chemical(_) => "chemical",
            GoCamNodeType::UnknownMRNA => "unknown_mrna",
            GoCamNodeType::Gene(_) => "gene",
            GoCamNodeType::MRNA(_) => "mrna",
            GoCamNodeType::ModifiedProtein(_) => "modified_protein",
            GoCamNodeType::Activity(GoCamActivity { enabler, inputs: _, outputs: _ }) => match enabler {
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
        if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = self.node_type {
            format!("{} [enabled by] {}", self.label, enabler.label())
        } else {
            self.label.to_owned()
        }
    }

    /// The label of the enabler, otherwise ""
    pub fn enabler_label(&self) -> &str {
        if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = self.node_type {
            enabler.label()
        } else {
            ""
        }
    }

    /// The label of the enabler, otherwise ""
    pub fn enabler_id(&self) -> &str {
        if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = self.node_type {
            enabler.id()
        } else {
            ""
        }
    }

    /// If this node is an activity, return the ID of the enabler.
    /// Otherwise return the ID of the node (i.e. the chemical,
    /// complex, modified protein or gene ID).
    pub fn db_id(&self) -> &str {
        if let GoCamNodeType::Activity(GoCamActivity { ref enabler, inputs: _, outputs: _ }) = self.node_type {
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
        if let GoCamNodeType::Activity(GoCamActivity { enabler: _, ref inputs, ref outputs }) = self.node_type {
            for has_input in inputs {
                write!(f, " [has input] {}", has_input)?;
            }
            for has_output in outputs {
                write!(f, " [has output] {}", has_output)?;
            }
        }
        for occurs_in in &self.occurs_in {
            write!(f, " [occurs in] {}", occurs_in)?;
        }
        if let Some(ref part_of_process) = self.part_of_process {
            write!(f, " [part of] {}", part_of_process)?;
        }
        if let Some(ref happens_during) = self.happens_during {
            write!(f, " [happens during] {}", happens_during)?;
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
            let node_type =
                if individual.individual_is_chemical() {
                    let chemical = GoCamChemical {
                         id: individual_type.id().to_owned(),
                         label: individual_type.label().to_owned(),
                         located_in: None,
                    };
                    GoCamNodeType::Chemical(chemical)
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
                                    let gene = GoCamGene::new(individual_type, &None);
                                    GoCamNodeType::Gene(gene)
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
            models.insert((model_id.clone(), model_title));
            let gocam_node = GoCamNode {
                individual_gocam_id: individual.id.clone(),
                node_id: individual_type.id.clone().unwrap_or_else(|| "NO_ID".to_owned()),
                label: individual_type.label.clone().unwrap_or_else(|| "NO_LABEL".to_owned()),
                node_type,
                occurs_in: BTreeSet::new(),
                part_of_process: None,
                happens_during: None,
                source_ids,
                original_model_id: Some(model_id),
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
            complex.has_part_genes.insert(complex_gene);
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
                        let facts = model.facts_of_subject(&object_individual.id);
                        let complex = facts.iter()
                            .find(|f| f.property_label == "part of")
                            .map(|f| complex_map.get(&f.object))
                            .flatten()
                            .cloned();

                        let gene = GoCamGene::new(object_type, &complex);
                        let enabler = GoCamEnabledBy::Gene(gene);

                        subject_node.node_type = GoCamNodeType::Activity(GoCamActivity {
                             enabler, inputs: BTreeSet::new(), outputs: BTreeSet::new(),
                        });
                    }
                    else if object_type_id.starts_with("CHEBI:") {
                        let enabler = GoCamEnabledBy::Chemical(object_type.into());
                        subject_node.node_type = GoCamNodeType::Activity(GoCamActivity {
                            enabler, inputs: BTreeSet::new(), outputs: BTreeSet::new(),
                        });
                    }
                    else if object_type_id.starts_with("GO:") || object_type_id.starts_with("ComplexPortal:") {
                        let complex = complex_map.get(&fact.object)
                            .expect(&format!("expected complex {}", fact.object))
                            .to_owned();
                        let enabler = GoCamEnabledBy::Complex(complex);
                        subject_node.node_type = GoCamNodeType::Activity(GoCamActivity {
                            enabler, inputs: BTreeSet::new(), outputs: BTreeSet::new(),
                        });
                    }
                    else if object_type_id.starts_with("PR:") {
                        let enabler = GoCamEnabledBy::ModifiedProtein(object_type.into());
                        subject_node.node_type = GoCamNodeType::Activity(GoCamActivity {
                            enabler, inputs: BTreeSet::new(), outputs: BTreeSet::new(),
                        });
                    }
                    else  {
                        eprintln!("{}: can't handle enabled by object: {} - {}",
                                  model.id(), object_type_id, object_individual.id);
                    }
                }
            },
            "located in" => {
                let GoCamNodeType::Chemical(ref mut chemical) = subject_node.node_type
                else {
                    panic!("located_in relation for non-chemical");
                };
                if chemical.located_in.is_some() {
                    panic!("{}: {} is located in multiple components", model.id(),
                           subject_node.description());
                }
                let located_in =
                    if object_individual.individual_is_complex() {
                        GoCamComponent::ComplexComponent(object_type.into())
                    } else {
                        GoCamComponent::OtherComponent(object_type.into())
                    };
                chemical.located_in = Some(located_in);
            },
            "occurs in" => {
                let occurs_in =
                    if object_individual.individual_is_complex() {
                        GoCamComponent::ComplexComponent(object_type.into())
                    } else {
                        GoCamComponent::OtherComponent(object_type.into())
                    };
                subject_node.occurs_in.insert(occurs_in);
            },
            "part of" => {
                let process = process_from_individual(object_individual, model);

                subject_node.part_of_process = Some(process);
            },
            "happens during" => {
                let during = process_from_individual(object_individual, model);
                subject_node.happens_during = Some(during);
            },
            &_ => (),
        }
    }

    for fact in model.facts() {
        let object_individual = model.fact_object(fact);

        let Some(object_type) = object_individual.get_individual_type()
        else {
            continue;
        };

        let object_node_located_in =
             if let Some(object_node) = node_map.get(&object_individual.id) {
                if let GoCamNodeType::Chemical(ref chemical) = object_node.node_type {
                    chemical.located_in.clone()
                } else {
                    None
                }
            } else {
                continue;
            };

        let object_node_occurs_in = {
            if let Some(ref object_node) = node_map.get(&object_individual.id) {
                object_node.occurs_in.clone()
            } else {
                panic!("internal error: can't find node for {}", object_individual.id);
            }
        };

        let Some(subject_node) = node_map.get_mut(&fact.subject)
        else {
            continue;
        };

        match fact.property_label.as_str() {
            "has input" => {
                let mut input: GoCamInput = object_type.into();
                input.occurs_in = object_node_occurs_in;
                input.located_in = object_node_located_in;
                match subject_node.node_type {
                    GoCamNodeType::Activity(GoCamActivity { enabler: _, ref mut inputs, outputs: _ }) => {
                        inputs.insert(input);
                    },
                    _ => (),
                }
            },
            "has output" => {
                let mut output: GoCamOutput = object_type.into();
                output.occurs_in = object_node_occurs_in;
                output.located_in = object_node_located_in;
                match subject_node.node_type {
                    GoCamNodeType::Activity(GoCamActivity { enabler: _, inputs: _, ref mut outputs }) => {
                        outputs.insert(output);
                    },
                    _ => (),
                }
            },
            &_ => (),
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
        let subject_id = fact.subject.as_str();
        let object_id = fact.object.as_str();

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
                   "GO:0140483 kinetochore adaptor activity enabled by: moa1 Spom (PomBase:SPAC15E1.07c) [occurs in] kinetochore (GO:0000776) [part of] meiotic centromeric cohesion protection in anaphase I (GO:1990813)");

    }

    #[test]
    fn test_genes_in_model() {
        let mut source = File::open("tests/data/gomodel_67f85f2b00003383.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:67f85f2b00003383");

        let expected_ids = vec![
            "PomBase:SPAC12B10.06c", "PomBase:SPAC140.01",
            "PomBase:SPAC1556.02c", "PomBase:SPAC664.12c",
            "PomBase:SPBC26H8.16", "PomBase:SPBP23A10.03c"
        ];

        let expected_genes_in_model: BTreeSet<String> =
            expected_ids.into_iter().map(String::from).collect();

        assert_eq!(expected_genes_in_model, model.genes_in_model());
    }

    #[test]
    fn test_modified_genes_in_model() {
        let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
        let mut model = parse_gocam_model(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:66187e4700001744");

        let expected_ids = vec![
            "PomBase:SPAC15E1.07c", "PomBase:SPAC23C11.16",
            "PomBase:SPAC664.01c", "PomBase:SPBC106.01",
            "PomBase:SPBC16H5.07c", "PomBase:SPBP35G2.03c",
            "PomBase:SPCC1322.12c", "PomBase:SPCC1020.02",
            "PomBase:SPCC622.08c", "PomBase:SPAC19G12.06c",
            "PomBase:SPBC29A10.14"
        ];

        let expected_genes_in_model: BTreeSet<String> =
            expected_ids.into_iter().map(String::from).collect();

        let mut pro_term_to_gene_map = HashMap::new();

        pro_term_to_gene_map.insert("PR:000059631".to_owned(),
                                    "PomBase:SPCC1020.02".to_owned());
        pro_term_to_gene_map.insert("PR:000027566".to_owned(),
                                    "PomBase:SPCC622.08c".to_owned());
        pro_term_to_gene_map.insert("PR:000027557".to_owned(),
                                    "PomBase:SPAC19G12.06c".to_owned());
        pro_term_to_gene_map.insert("PR:000050512".to_owned(),
                                    "PomBase:SPBC29A10.14".to_owned());

        model.add_pro_term_to_gene_map(&pro_term_to_gene_map);

        let actual_genes = model.genes_in_model();

        assert_eq!(expected_genes_in_model, actual_genes);
    }

    #[test]
    fn test_genes_enabling_activities() {
        let mut source = File::open("tests/data/gomodel_67f85f2b00003383.json").unwrap();
        let mut gene_name_map = HashMap::new();
        gene_name_map.insert("SPAC12B10.06c".to_owned(), "sdh5".to_owned());
        let mut model = parse_gocam_model(&mut source).unwrap();
        model.add_gene_name_map(&gene_name_map);

        assert_eq!(model.id(), "gomodel:67f85f2b00003383");

        let expected = vec![
            ("PomBase:SPAC12B10.06c", Some("sdh5".to_owned())),
            ("PomBase:SPAC664.12c", None),
            ("PomBase:SPBC26H8.16", None),
            ("PomBase:SPBP23A10.03c", None)
        ];
        let expected_genes_in_model: HashMap<String, Option<String>> =
            expected.into_iter()
            .map(|(id, name)| {
                (id.to_owned(), name)
            }).collect();

        assert_eq!(expected_genes_in_model,
                   model.genes_enabling_activities());
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

        let test_out_direction = model.node_rel_direction(test_out_node_idx);
        assert_eq!(test_out_direction, GoCamDirection::Outgoing);

        let test_in_direction = model.node_rel_direction(test_in_node_idx);
        assert_eq!(test_in_direction, GoCamDirection::IncomingConstitutivelyUpstream);
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

        let overlaps = find_overlaps(&[model1, model2]);

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
                GoCamNodeType::Activity(GoCamActivity { enabler: ref enabled_by, .. }) => {
                    activity_enabled_by_ids.insert(enabled_by.id().to_owned());
                },
                GoCamNodeType::Chemical(_) => {
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
    fn merge_test_mss51() {
        let mut source1 = File::open("tests/data/gomodel_67086be200000519.json").unwrap();
        let model1 = parse_gocam_model(&mut source1).unwrap();

        let mut source2 = File::open("tests/data/gomodel_67e5e74400003073.json").unwrap();
        let model2 = parse_gocam_model(&mut source2).unwrap();

        let merged = GoCamModel::merge_models("new_id", "new_title",
                                              &[model1, model2]).unwrap();

        assert_eq!(merged.node_iterator().count(), 89);

        // find IDs in merged nodes
        let merged_ids: HashSet<_> =
            merged.node_iterator().filter_map(|(_, node)| if node.models.len() >= 2 {
                Some((node.label.clone(), node.db_id().to_owned()))
            } else {
                None
            })
            .collect();

        let mut expected_ids = HashSet::new();
        expected_ids.insert(("gene product or complex activity".to_owned(), "PomBase:SPAC25B8.04c".to_owned()));
        expected_ids.insert(("cox1 Spom".to_owned(), "PomBase:SPMIT.01".to_owned()));

        assert_eq!(merged_ids, expected_ids);
    }

    #[test]
    fn remove_chemicals_test() {
        let mut source1 = File::open("tests/data/gomodel_66a3e0bb00001342.json").unwrap();
        let model = parse_gocam_model(&mut source1).unwrap();

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Chemicals);
        let new_model = model.remove_nodes(remove_types);

        assert_eq!(new_model.node_count(), 18);
    }

    #[test]
    fn remove_chemicals_inputs_outputs_only_test() {
        // test removing chemicals but remove those that have
        // relations that aren't input or outputs
        let mut source1 = File::open("tests/data/gomodel_66187e4700003150.json").unwrap();
        let model = parse_gocam_model(&mut source1).unwrap();

        assert_eq!(model.node_count(), 20);

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Chemicals);
        let new_model = model.remove_nodes(remove_types);

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Targets);
        assert_eq!(new_model.node_count(), 18);

        let new_model2 = model.remove_nodes(remove_types);
        assert_eq!(new_model2.node_count(), 20);

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Chemicals);
        remove_types.insert(RemoveType::Targets);
        let new_model3 = model.remove_nodes(remove_types);
        assert_eq!(new_model3.node_count(), 18);
    }

    #[test]
    fn merge_test_mss51_remove_chemicals() {
        let mut source1 = File::open("tests/data/gomodel_67086be200000519.json").unwrap();
        let model1 = parse_gocam_model(&mut source1).unwrap();

        let mut source2 = File::open("tests/data/gomodel_67e5e74400003073.json").unwrap();
        let model2 = parse_gocam_model(&mut source2).unwrap();

        let merged = GoCamModel::merge_models("new_id", "new_title",
                                              &[model1, model2]).unwrap();

        let chemical_nodes_iter = merged.node_iterator();
        let chemical_count = chemical_nodes_iter
            .filter(|(_, n)| n.node_type.is_chemical())
            .count();
        assert_eq!(chemical_count, 6);

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Chemicals);
        let new_model = merged.remove_nodes(remove_types);

        let chemical_nodes_iter = new_model.node_iterator();
        let chemical_count = chemical_nodes_iter
            .filter(|(_, n)| n.node_type.is_chemical())
            .count();
        assert_eq!(chemical_count, 0);

        assert_eq!(new_model.node_iterator().count(), 83);

        // find IDs in merged nodes
        let merged_ids: HashSet<_> =
            new_model.node_iterator().filter_map(|(_, node)| if node.models.len() >= 2 {
                Some((node.label.clone(), node.db_id().to_owned()))
            } else {
                None
            })
            .collect();

        let mut expected_ids = HashSet::new();
        expected_ids.insert(("gene product or complex activity".to_owned(), "PomBase:SPAC25B8.04c".to_owned()));
        expected_ids.insert(("cox1 Spom".to_owned(), "PomBase:SPMIT.01".to_owned()));

        assert_eq!(merged_ids, expected_ids);

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Targets);
        let no_inputs_model = merged.remove_nodes(remove_types);

        let input_nodes_iter = no_inputs_model.node_iterator();
        let inputs_count = input_nodes_iter
            .filter(|(_, n)| !n.is_activity())
            .count();
        assert_eq!(inputs_count, 6);

        let mut remove_types = HashSet::new();
        remove_types.insert(RemoveType::Chemicals);
        remove_types.insert(RemoveType::Targets);
        let no_inputs_model = merged.remove_nodes(remove_types);

        let input_nodes_iter = no_inputs_model.node_iterator();
        let inputs_count = input_nodes_iter
            .filter(|(_, n)| !n.is_activity())
            .count();
        assert_eq!(inputs_count, 0);
    }

    #[test]
    fn merge_test_three_models() {
        let mut source1 = File::open("tests/data/gomodel_665912ed00000459.json").unwrap();
        let model1 = parse_gocam_model(&mut source1).unwrap();

        let mut source2 = File::open("tests/data/gomodel_671ae02600003596.json").unwrap();
        let model2 = parse_gocam_model(&mut source2).unwrap();

        let mut source3 = File::open("tests/data/gomodel_665912ed00000192.json").unwrap();
        let model3 = parse_gocam_model(&mut source3).unwrap();

        let merged = GoCamModel::merge_models("new_id", "new_title",
                                              &[model1, model2, model3]).unwrap();

        use petgraph::visit::NodeRef;

        for node in merged.node_iterator() {
            if node.weight().node_id == "GO:0004582" {
                assert_eq!(node.weight().models.iter().count(), 3);
            }
        }
    }

    #[test]
    fn retain_largest_subgraph_test() {
        let mut source = File::open("tests/data/gomodel_671ae02600003596.json").unwrap();
        let mut model = parse_gocam_model(&mut source).unwrap();

        let chemical_node_indexes: BTreeSet<_> = model.node_iterator()
            .filter(|(_, node)| {
                !node.is_activity()
            })
            .map(|(node_idx, _)| node_idx)
            .collect();

        assert_eq!(model.graph.node_count(), 9);

        model.graph.retain_nodes(|_, node_idx| {
            !chemical_node_indexes.contains(&node_idx)
        });

        assert_eq!(model.graph.node_count(), 6);

        let largest_subgraph_model = model.retain_largest_subgraph();

        assert_eq!(largest_subgraph_model.graph.node_count(), 3);
    }

    #[test]
    fn retain_enabling_genes_test() {
        let mut source = File::open("tests/data/gomodel_671ae02600003596.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();

        let chemical_node_count = model.node_iterator()
            .filter(|(_, node)| {
                node.type_string() == "chemical"
            })
            .count();
        assert_eq!(chemical_node_count, 3);

        assert_eq!(model.graph.node_count(), 9);

        let mut remove_list = BTreeSet::new();
        remove_list.insert("SPBC21B10.11".to_owned());
        remove_list.insert("SPBC1677.02".to_owned());
        remove_list.insert("SPAC31G5.16c".to_owned());

        let model = model.retain_enabling_genes(&remove_list);

        assert_eq!(model.graph.node_count(), 3);

        let chemical_node_count = model.node_iterator()
            .filter(|(_, node)| {
                node.type_string() == "chemical"
            })
            .count();
        assert_eq!(chemical_node_count, 0);
    }

    #[test]
    fn model_activity_enabled_by_test() {
        let mut source = File::open("tests/data/gomodel_671ae02600003596.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();

        assert_eq!(model.graph.node_count(), 9);

        let mut test_genes = BTreeSet::new();
        test_genes.insert("SPAC11D3.09".to_owned());
        test_genes.insert("SPBC1677.02".to_owned());
        test_genes.insert("SPAC31G5.16c".to_owned());

        assert!(model.model_activity_enabled_by(&test_genes));

        test_genes.clear();
        test_genes.insert("SPBC18H10.20c".to_owned());

        assert!(!model.model_activity_enabled_by(&test_genes));
    }

    #[test]
    fn location_of_chemical() {
        let mut source = File::open("tests/data/gomodel_66187e4700003150.json").unwrap();
        let model = parse_gocam_model(&mut source).unwrap();

        let (_, cgs2_node) =
            model.node_iterator().filter(|(_, node)| {
                node.node_id == "GO:0004016"
            }).next().unwrap();
        assert!(cgs2_node.is_activity());

        let GoCamNodeType::Activity(GoCamActivity { enabler: _, ref inputs, outputs: _ }) = cgs2_node.node_type
        else {
            panic!();
        };

        assert_eq!(inputs.len(), 1);

        let input = inputs.iter().next().unwrap();
        assert_eq!(input.located_in.clone().unwrap().to_string(), "cytosol (GO:0005829)");
    }
}
