use std::{cmp::Ordering, collections::{BTreeMap, BTreeSet, HashMap, HashSet}, fmt::{self, Display}, io::{BufReader, Read}};

extern crate serde_json;
#[macro_use] extern crate serde_derive;
use phf::phf_map;

use anyhow::{Result, anyhow};

use petgraph::{graph::NodeReferences, visit::{EdgeRef, IntoNodeReferences}, Graph};

static REL_NAMES: phf::Map<&'static str, &'static str> = phf_map! {
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
pub type FactId = String;
pub type IndividualId = String;
pub type PropertyId = String;
pub type AnnotationKey = String;
pub type AnnotationValue = String;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Annotation {
    pub key: AnnotationKey,
    pub value: AnnotationValue,
    #[serde(skip_serializing_if="Option::is_none")]
    #[serde(rename = "value-type")]
    pub value_type: Option<String>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Fact {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub annotations: Vec<Annotation>,
    pub subject: IndividualId,
    pub object: IndividualId,
    pub property: PropertyId,
    #[serde(rename = "property-label")]
    pub property_label: String,
}

impl Fact {
    pub fn id(&self) -> FactId {
        format!("{}-{}-{}", self.subject, self.property, self.object)
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct IndividualType {
    pub id: Option<String>,
    pub label: Option<String>,
    #[serde(rename = "type")]
    pub type_string: String,
}

const MOLECULAR_FUNCTION_ID: &str = "GO:0003674";
/*
const CELLULAR_COMPONENT_ID: &str = "GO:0032991";
const BIOLOGICAL_PROCESS_ID: &str = "GO:0008150";
*/
const PROTEIN_CONTAINING_COMPLEX_ID: &str = "GO:0032991";
const CHEBI_PROTEIN_ID: &str = "CHEBI:36080";
const CHEBI_CHEMICAL_ENTITY_ID: &str = "CHEBI:24431";

impl IndividualType {
    pub fn id(&self) -> &str {
        self.id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_ID")
    }

    pub fn label(&self) -> &str {
        self.label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_LABEL")
    }

    pub fn label_or_id(&self) -> &str {
        self.label.as_ref().map(|s| s.as_str())
            .or(self.id.as_ref().map(|s| s.as_str()))
            .unwrap_or("UNKNOWN")
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Individual {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub annotations: Vec<Annotation>,
    pub id: IndividualId,
    #[serde(rename = "type")]
    pub types: Vec<IndividualType>,
    #[serde(rename = "root-type")]
    pub root_types: Vec<IndividualType>,
}

impl Individual {
    pub fn has_root_term(&self, term_id: &str) -> bool {
        for individual_type in &self.root_types {
            if let Some(ref individual_type_id) = individual_type.id {
                if individual_type_id == term_id {
                    return true;
                }
            }
        }

        false
    }

    pub fn individual_is_activity(&self) -> bool {
        self.has_root_term(MOLECULAR_FUNCTION_ID)
    }

    /*
    fn individual_is_component(individual: &Individual) -> bool {
    has_root_term(individual, CELLULAR_COMPONENT_ID)
    }

    fn individual_is_process(individual: &Individual) -> bool {
    has_root_term(individual, BIOLOGICAL_PROCESS_ID)
    }
    */

    pub fn individual_is_complex(&self) -> bool {
        self.has_root_term(PROTEIN_CONTAINING_COMPLEX_ID)
    }

    pub fn individual_is_gene(&self) -> bool {
        if !self.has_root_term(CHEBI_CHEMICAL_ENTITY_ID) {
            return false;
        }

        if self.has_root_term(CHEBI_PROTEIN_ID) {
            return false;
        }

        let Some(individual_type) = self.get_individual_type()
        else {
            return false;
        };

        if let Some(ref id) = individual_type.id {
            if id.starts_with("CHEBI:") {
                return false;
            }
        }

        true
    }

    pub fn individual_is_chemical(&self) -> bool {
        if !self.has_root_term(CHEBI_CHEMICAL_ENTITY_ID) {
            return false;
        }

        let Some(individual_type) = self.get_individual_type()
        else {
            return false;
        };

        if let Some(ref id) = individual_type.id {
            if id.starts_with("CHEBI:") {
                return true;
            }
        }

        false
    }

    pub fn get_individual_type(&self) -> Option<&IndividualType> {
        self.types.get(0)
    }

    pub fn individual_is_unknown_protein(&self) -> bool {
        let Some(individual_type) = self.get_individual_type()
        else {
            return false;
        };

        if let Some(ref individual_type_id) = individual_type.id {
            if individual_type_id == CHEBI_PROTEIN_ID {
                return true;
            }
        }

        false
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct SerdeModel {
    annotations: Vec<Annotation>,
    id: ModelId,
    facts: Vec<Fact>,
    individuals: Vec<Individual>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct GoCamRawModel {
    _annotations: Vec<Annotation>,
    _id: ModelId,
    _facts: BTreeMap<FactId, Fact>,
    _individuals: BTreeMap<IndividualId, Individual>,
}

fn annotation_values(annotations: Box<dyn Iterator<Item = &Annotation> + '_>, key: &str)
                     -> Vec<String>
{
    annotations.filter(|annotation| annotation.key == key)
        .map(|annotation| &annotation.value)
        .map(|s| s.replace(&['\n', '\t'], " ").trim_matches(' ').to_owned())
        .collect()
}

impl GoCamRawModel {
    pub fn id(&self) -> &ModelId {
        &self._id
    }

    pub fn annotations(&self) -> Box<dyn Iterator<Item = &Annotation> + '_> {
        Box::new(self._annotations.iter())
    }

    pub fn title(&self) -> String {
        annotation_values(self.annotations(), "title")
            .join(",")
    }

    pub fn date(&self) -> String {
        annotation_values(self.annotations(), "date")
            .join(",")
    }

    pub fn taxon(&self) -> String {
        annotation_values(self.annotations(), "https://w3id.org/biolink/vocab/in_taxon")
            .join(",")
    }

    pub fn contributors(&self) -> BTreeSet<String> {
        BTreeSet::from_iter(annotation_values(self.annotations(), "contributor"))
    }

    pub fn facts(&self) -> Box<dyn Iterator<Item = &Fact> + '_>  {
        Box::new(self._facts.values())
    }

    pub fn individuals(&self) -> Box<dyn Iterator<Item = &Individual> + '_> {
        Box::new(self._individuals.values())
    }

    pub fn fact_subject(&self, fact: &Fact) -> &Individual {
        self.get_individual(&fact.subject)
    }

    pub fn fact_object(&self, fact: &Fact) -> &Individual {
        self.get_individual(&fact.object)
    }

    pub fn get_individual(&self, individual_id: &IndividualId)
        -> &Individual
    {
        self._individuals.get(individual_id)
            .expect(&format!("can't find individual: {}",
                             individual_id))
    }
}

pub type GoCamGraph = Graph::<GoCamNode, GoCamEdge>;

pub type GoCamGene = IndividualType;
pub type GoCamChemical = IndividualType;
pub type GoCamModifiedProtein = IndividualType;
pub type GoCamProcess = IndividualType;
pub type GoCamInput = IndividualType;
pub type GoCamOutput = IndividualType;
pub type GoCamGeneIdentifier = String;

pub type GoCamNodeMap = BTreeMap<IndividualId, GoCamNode>;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum GoCamComponent {
    ComplexComponent(IndividualType),
    OtherComponent(IndividualType),
}

impl GoCamComponent {
    pub fn id(&self) -> &str {
        match self {
            GoCamComponent::ComplexComponent(individual_type) |
            GoCamComponent::OtherComponent(individual_type) => {
                individual_type.id()
            }
        }
    }

    pub fn label(&self) -> &str {
        match self {
            GoCamComponent::ComplexComponent(individual_type) |
            GoCamComponent::OtherComponent(individual_type) => {
                individual_type.label()
            }
        }
    }

    pub fn label_or_id(&self) -> &str {
        match self {
            GoCamComponent::ComplexComponent(individual_type) |
            GoCamComponent::OtherComponent(individual_type) => {
                individual_type.label_or_id()
            }
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamModel {
    id: String,
    title: String,
    taxon: String,
    date: String,
    contributors: BTreeSet<String>,
    graph: GoCamGraph,
}

fn check_model_taxons(models: &[&GoCamModel]) -> Result<String> {
    let mut model_iter = models.iter();

    let Some(first_model) = model_iter.next()
    else {
        return Err(anyhow!("no models passed to check_model_taxons()"));
    };

    for this_model in model_iter {
        if first_model.taxon != this_model.taxon {
            return Err(anyhow!("mismatched taxons: {} ({}) != {} ({})",
                               first_model.taxon, first_model.id,
                               this_model.taxon, this_model.id));
        }
    }

    Ok(first_model.taxon().to_owned())
}

impl GoCamModel {
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

    pub fn graph(&self) -> &GoCamGraph {
        &self.graph
    }

    pub fn node_iterator(&self) -> NodeIterator {
        NodeIterator {
            node_refs: self.graph().node_references(),
        }
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn title(&self) -> &str {
        &self.title
    }

    pub fn taxon(&self) -> &str {
        &self.taxon
    }

    pub fn date(&self) -> &str {
        &self.date
    }

    pub fn contributors(&self) -> &BTreeSet<String> {
        &self.contributors
    }

    pub fn find_activity_overlaps(models: &[&GoCamModel])
        -> Vec<GoCamNodeOverlap>
    {
        let mut seen_activities = HashMap::new();

        let make_key = |node: &GoCamNode| {
            let Some(ref occurs_in_type) = node.occurs_in
            else {
                return None;
            };
            let Some(ref part_of_process_type) = node.part_of_process
            else {
                return None;
            };

            let (node_type, enabled_by_type, enabled_by_id, _) =
                node.node_type_summary_strings();

            Some((node.node_id.to_owned(),
                  node.description(),
                  node_type.to_owned(),
                  enabled_by_type.to_owned(), enabled_by_id.to_owned(),
                  part_of_process_type.id().to_owned(),
                  part_of_process_type.label().to_owned(),
                  occurs_in_type.id().to_owned(),
                  occurs_in_type.label().to_owned()))
        };

        for model in models {
            for node in model.node_iterator() {

                let Some(key) = make_key(node)
                else {
                    continue;
                };

                seen_activities
                    .entry(key)
                    .or_insert_with(Vec::new)
                    .push(model.id());
            }
        }

        let mut ret = vec![];

        for (key, model_ids) in seen_activities {
            if model_ids.len() > 1 {
                let (node_id, node_description, node_type, enabled_by_type, enabled_by_id,
                     part_of_process_id, part_of_process_label,
                     occurs_in_id, occurs_in_label) = key;
                let node_overlap = GoCamNodeOverlap {
                    node_id,
                    node_description,
                    node_type,
                    enabled_by_type,
                    enabled_by_id,
                    part_of_process_id,
                    part_of_process_label,
                    occurs_in_id,
                    occurs_in_label,
                    model_ids: model_ids.iter().map(|&s| s.to_owned()).collect(),
                };
                ret.push(node_overlap);
            }
        }

        ret.sort_by(|a, b| {
            let ord = a.model_ids.cmp(&b.model_ids);

            if ord == Ordering::Equal {
                a.node_description.cmp(&b.node_description)
            } else {
                ord
            }
        });

        ret
    }

    pub fn merge_models(new_id: &str, new_title: &str, models: &[&GoCamModel])
        -> Result<GoCamModel>
    {
        let mut merged_graph = GoCamGraph::new();

        let mut idx_map = HashMap::new();
        let taxon = check_model_taxons(models)?;

        let mut contributors = BTreeSet::new();

        for model in models {
            for (old_idx, node) in model.graph().node_references() {
                let new_idx = merged_graph.add_node(node.to_owned());
                idx_map.insert(old_idx, new_idx);
            }

            for edge_ref in model.graph().edge_references() {
                let edge = edge_ref.weight();

                let old_source_idx = edge_ref.source();
                let old_target_idx = edge_ref.target();

                let new_source_idx = idx_map.get(&old_source_idx).unwrap();
                let new_target_idx = idx_map.get(&old_target_idx).unwrap();

                merged_graph.add_edge(*new_source_idx, *new_target_idx, edge.to_owned());
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
}

#[derive(Clone, Debug)]
pub struct GoCamNodeOverlap {
    pub node_id: String,
    pub node_description: String,
    pub node_type: String,
    pub enabled_by_type: String,
    pub enabled_by_id: String,
    pub part_of_process_id: String,
    pub part_of_process_label: String,
    pub occurs_in_id: String,
    pub occurs_in_label: String,
    pub model_ids: BTreeSet<ModelId>,
}

pub struct NodeIterator<'a> {
    node_refs: NodeReferences<'a, GoCamNode>,
}

impl<'a> Iterator for NodeIterator<'a> {
    type Item = &'a GoCamNode;

    fn next(&mut self) -> Option<Self::Item> {
        self.node_refs.next().map(|(_, node)| node)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamComplex {
    pub id: Option<String>,
    pub label: Option<String>,
    pub has_part_genes: Vec<GoCamGeneIdentifier>,
}

impl GoCamComplex {
    pub fn id(&self) -> &str {
        self.id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_ID")
    }

    pub fn label(&self) -> &str {
        self.label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_LABEL")
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum GoCamEnabledBy {
    Complex(GoCamComplex),
    Gene(GoCamGene),
    Chemical(GoCamChemical),
    ModifiedProtein(GoCamModifiedProtein),
}

impl GoCamEnabledBy {
    pub fn id(&self) -> &str {
        let maybe_id = match self {
            GoCamEnabledBy::Complex(complex) => &complex.id,
            GoCamEnabledBy::Gene(gene) => &gene.id,
            GoCamEnabledBy::Chemical(chemical) => &chemical.id,
            GoCamEnabledBy::ModifiedProtein(modified_protein) => &modified_protein.id,
        };
        maybe_id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }

    pub fn label(&self) -> &str {
        let maybe_label = match self {
            GoCamEnabledBy::Complex(complex) => &complex.label,
            GoCamEnabledBy::Gene(gene) => &gene.label,
            GoCamEnabledBy::Chemical(chemical) => &chemical.label,
            GoCamEnabledBy::ModifiedProtein(modified_protein) => &modified_protein.label,
        };
        maybe_label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN")
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum GoCamNodeType {
    Unknown,
    Chemical,
    Gene(GoCamGene),
    ModifiedProtein(GoCamModifiedProtein),
    Activity(GoCamEnabledBy),
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct GoCamNode {
    pub individual_gocam_id: IndividualId,
    pub node_id: String,
    pub label: String,
    pub node_type: GoCamNodeType,
    pub has_input: Vec<GoCamInput>,
    pub has_output: Vec<GoCamOutput>,
    pub located_in: Option<GoCamComponent>,
    pub occurs_in: Option<GoCamComponent>,
    pub part_of_process: Option<GoCamProcess>,
}

impl GoCamNode {
    fn node_type_summary_strings(&self) -> (&str, &str, &str, &str) {
        match &self.node_type {
            GoCamNodeType::Unknown => ("unknown", "unknown", "unknown", "unknown"),
            GoCamNodeType::Chemical => ("chemical", "", "", ""),
            GoCamNodeType::Gene(_) => ("gene", "", "", ""),
            GoCamNodeType::ModifiedProtein(_) => ("modified_protein", "", "", ""),
            GoCamNodeType::Activity(enabled_by) => match enabled_by {
                GoCamEnabledBy::Chemical(chem) => ("activity", "chemical", chem.id(), chem.label()),
                GoCamEnabledBy::Gene(gene) => ("activity", "gene", gene.id(), gene.label()),
                GoCamEnabledBy::ModifiedProtein(prot) => ("activity", "modified_protein", prot.id(), prot.label()),
                GoCamEnabledBy::Complex(complex) => ("activity", "complex", complex.id(), complex.label()),
            }
        }
    }
}

impl Display for GoCamNode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t", self.node_id)?;
        write!(f, "{}\t", self.label)?;

        let (node_type, enabled_by_type, enabled_by_id, enabled_by_label) =
            self.node_type_summary_strings();

        write!(f, "{}\t{}\t{}\t{}\t", node_type, enabled_by_type, enabled_by_id, enabled_by_label)?;

        if let Some(ref part_of_process) = self.part_of_process {
            write!(f, "{}\t", part_of_process.label_or_id())?;
        } else {
            write!(f, "\t")?;
        }
        let has_input_string =
            self.has_input.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if has_input_string.len() > 0 {
            write!(f, "{}\t", has_input_string)?;
        } else {
            write!(f, "\t")?;
        }
        let has_output_string =
            self.has_output.iter().map(|l| l.label_or_id()).collect::<Vec<_>>().join(",");
        if has_output_string.len() > 0 {
            write!(f, "{}\t", has_output_string)?;
        } else {
            write!(f, "\t")?;
        }

        if let Some(ref occurs_in) = self.occurs_in {
            write!(f, "{}\t", occurs_in.label_or_id())?;
        } else {
            write!(f, "\t")?;
        }
        if let Some(ref located_in) = self.located_in {
            write!(f, "{}", located_in.label_or_id())?;
        } else {
            write!(f, "\t")?;
        }

        Ok(())
    }
}

impl GoCamNode {
    pub fn type_string(&self) -> &str {
        match &self.node_type {
            GoCamNodeType::Unknown => "unknown",
            GoCamNodeType::Chemical => "chemical",
            GoCamNodeType::Gene(_) => "gene",
            GoCamNodeType::ModifiedProtein(_) => "modified_protein",
            GoCamNodeType::Activity(activity) => match activity {
                GoCamEnabledBy::Chemical(_) => "enabled_by_chemical",
                GoCamEnabledBy::Gene(_) => "enabled_by_gene",
                GoCamEnabledBy::ModifiedProtein(_) => "enabled_by_modified_protein",
                GoCamEnabledBy::Complex(_) => "enabled_by_complex",
            }
        }
    }

    pub fn description(&self) -> String {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            format!("{} [enabled by] {}", self.label, enabler.label())
        } else {
            self.label.to_owned()
        }
    }

    pub fn enabler_label(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.label()
        } else {
            ""
        }
    }

    pub fn enabler_id(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.id()
        } else {
            ""
        }
    }

    pub fn db_id(&self) -> &str {
        if let GoCamNodeType::Activity(ref enabler) = self.node_type {
            enabler.id()
        } else {
            &self.node_id
        }
    }
}

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


fn is_gene_id(identifier: &str) -> bool {
    ["PomBase:", "FB:", "UniProtKB:", "MGI:", "WB:", "RGD:", "RefSeq:",
     "Xenbase:", "SGD:", "ZFIN:", "RNAcentral:", "EMAPA:"]
        .iter().any(|s| identifier.starts_with(*s))
}

fn is_connecting_fact(rel_name: &str) -> bool {
    ["causally upstream of, negative effect",
     "causally upstream of, positive effect",
     "provides input for",
     "directly provides input for",
     "removes input for",
     "constitutively upstream of",
     "directly negatively regulates",
     "directly positively regulates",
     "indirectly negatively regulates",
     "indirectly positively regulates",
     "has small molecular activator",
     "has small molecular inhibitor",
     "is small molecule activator of",
     "is small molecule inhibitor of",
     "input of",
     "has output"]
        .iter().any(|s| rel_name == *s)
}

fn gocam_parse_raw(source: &mut dyn Read) -> Result<SerdeModel> {
    let reader = BufReader::new(source);

    let raw_model: SerdeModel = serde_json::from_reader(reader)?;

    Ok(raw_model)
}


/// Parses a GO-CAM model from a stream
///
/// # Example:
///
/// ```
/// use std::fs::File;
/// use pombase_gocam::gocam_parse;
///
/// let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
/// let model = gocam_parse(&mut source).unwrap();
/// assert_eq!(model.id(), "gomodel:66187e4700001744");
///
/// for fact in model.facts() {
///   let subject_id = &fact.subject;
///   println!("subject_id: {}", subject_id);
///   let subject_individual = model.get_individual(subject_id);
///   let first_type = &subject_individual.types[0];
///   if let Some(ref label) = first_type.label {
///     println!("type label: {}", label);
///   }
/// }
/// ```
pub fn gocam_parse(source: &mut dyn Read) -> Result<GoCamRawModel> {
    let raw_model = gocam_parse_raw(source)?;

    let mut fact_map = BTreeMap::new();
    let mut individual_map = BTreeMap::new();

    for mut fact in raw_model.facts.into_iter() {
        if let Some(&rel_name) = REL_NAMES.get(&fact.property) {
            fact.property_label = rel_name.to_owned();
        }
        fact_map.insert(fact.id(), fact);
    }

    for individual in raw_model.individuals.into_iter() {
        individual_map.insert(individual.id.clone(), individual);
    }

    Ok(GoCamRawModel {
        _annotations: raw_model.annotations,
        _id: raw_model.id,
        _facts: fact_map,
        _individuals: individual_map,
    })
}

fn make_nodes(model: &GoCamRawModel) -> GoCamNodeMap {
    // genes that are the object of a has_input or has_output relation
    let mut bare_genes_and_modified_proteins = HashSet::new();

    for fact in model.facts() {
        if fact.property_label == "has input" ||
            fact.property_label == "has output" {
                bare_genes_and_modified_proteins.insert(fact.object.clone());
            }
    }

    let mut node_map = BTreeMap::new();

    for individual in model.individuals() {
        if individual.individual_is_activity() ||
            bare_genes_and_modified_proteins.contains(&individual.id) ||
            individual.individual_is_chemical() &&
            !individual.individual_is_unknown_protein()
        {
            let Some(individual_type) = individual.get_individual_type()
            else {
                continue;
            };
            let detail =
                if individual.individual_is_chemical() {
                    GoCamNodeType::Chemical
                } else {
                    if bare_genes_and_modified_proteins.contains(&individual.id) {
                        if individual_type.id().starts_with("PR:") {
                            GoCamNodeType::ModifiedProtein(individual_type.clone())
                        } else {
                            GoCamNodeType::Gene(individual_type.clone())
                        }
                    } else {
                        GoCamNodeType::Unknown
                    }
                };
            let gocam_node = GoCamNode {
                individual_gocam_id: individual.id.clone(),
                node_id: individual_type.id.clone().unwrap_or_else(|| "NO_ID".to_owned()),
                label: individual_type.label.clone().unwrap_or_else(|| "NO_LABEL".to_owned()),
                node_type: detail,
                has_input: vec![],
                has_output: vec![],
                located_in: None,
                occurs_in: None,
                part_of_process: None,
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

            let complex = GoCamComplex {
                id: complex_type.id.clone(),
                label: complex_type.label.clone(),
                has_part_genes: vec![],
            };

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
                        let gene_enabler = GoCamEnabledBy::Gene(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(gene_enabler);
                    }
                    else if object_type_id.starts_with("CHEBI:") {
                        let chemical_enabler = GoCamEnabledBy::Chemical(object_type.clone());
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
                        let modified_protein_enabler = GoCamEnabledBy::ModifiedProtein(object_type.clone());
                        subject_node.node_type = GoCamNodeType::Activity(modified_protein_enabler);
                    }
                    else  {
                        eprintln!("can't handle enabled by object: {} - {}", object_type_id, object_individual.id);
                    }
                }
            },
            "has input" => {
                subject_node.has_input.push(object_type.clone());
            },
            "has output" => {
                subject_node.has_output.push(object_type.clone());
            },
            "located in" => {
                if subject_node.located_in.is_some() {
                    panic!("{}: {} is located in multiple components", model.id(),
                           subject_node.description());
                }
                let located_in =
                    if object_individual.individual_is_complex() {
                        GoCamComponent::ComplexComponent(object_type.clone())
                    } else {
                        GoCamComponent::OtherComponent(object_type.clone())
                    };
                subject_node.located_in = Some(located_in);
            },
            "occurs in" => {
                if subject_node.occurs_in.is_some() {
                    panic!("{}: {} is occurs in multiple components", model.id(),
                           subject_node.description());
                }
                let occurs_in =
                    if object_individual.individual_is_complex() {
                        GoCamComponent::ComplexComponent(object_type.clone())
                    } else {
                        GoCamComponent::OtherComponent(object_type.clone())
                    };
                subject_node.occurs_in = Some(occurs_in);
            },
            "part of" => {
                subject_node.part_of_process = Some(object_type.clone());
            },
            &_ => {
                // eprintln!("ignoring rel from fact: {} {}", fact.property_label, fact.id());
            }
        }
    }

    let mut connectons = vec![];

    for fact in model.facts() {
        if is_connecting_fact(fact.property_label.as_str()) {
            let connection = (fact.subject.clone(), fact.property_label.clone(),
                              fact.object.clone());
            connectons.push(connection);
        }
    }

    node_map
}

pub fn make_gocam_model(source: &mut dyn Read) -> Result<GoCamModel> {
    let raw_model = gocam_parse(source)?;

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
    fn parse_raw_test() {
        let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
        let model = gocam_parse(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:66187e4700001744");
        assert_eq!(model.facts().count(), 42);
        assert_eq!(model.individuals().count(), 82);

        assert_eq!(model.title(), "meiotic cohesion protection in anaphase I (GO:1990813)");
        assert_eq!(model.taxon(), "NCBITaxon:4896");

        let first_fact = model.facts().next().unwrap();
        assert_eq!(first_fact.property, "BFO:0000050");
        assert_eq!(first_fact.property_label, "part of");
        assert_eq!(first_fact.id(), "gomodel:66187e4700001744/66187e4700001758-BFO:0000050-gomodel:66187e4700001744/66187e4700001760");

        let fact_object = model.fact_object(&first_fact);
        assert_eq!(first_fact.object, fact_object.id);

        let object_first_type = &fact_object.types[0];
        assert_eq!(object_first_type.label.clone().unwrap(),
                   "meiotic centromeric cohesion protection in anaphase I");

        let object_first_root_type = &fact_object.root_types[0];
        assert_eq!(object_first_root_type.id.clone().unwrap(),
                   "GO:0008150");
    }


    #[test]
    fn parse_test() {
        let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
        let model = make_gocam_model(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:66187e4700001744");

        assert_eq!(model.title(), "meiotic cohesion protection in anaphase I (GO:1990813)");
        assert_eq!(model.taxon(), "NCBITaxon:4896");

        let first_node = model.node_iterator().next().unwrap();

        assert_eq!(first_node.node_id, "GO:0140483");

        assert_eq!(model.node_iterator().count(), 12);
    }

    #[test]
    fn merge_test() {
        let mut source1 = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
        let model1 = make_gocam_model(&mut source1).unwrap();
        assert_eq!(model1.id(), "gomodel:66187e4700001744");

        assert_eq!(model1.node_iterator().count(), 12);

        let mut source2 = File::open("tests/data/gomodel:665912ed00000015.json").unwrap();
        let model2 = make_gocam_model(&mut source2).unwrap();
        assert_eq!(model2.id(), "gomodel:665912ed00000015");

        assert_eq!(model2.node_iterator().count(), 25);

        let merged = GoCamModel::merge_models("new_id", "new_title",
                                              &[&model1, &model2]).unwrap();

        assert_eq!(merged.node_iterator().count(), 37);
    }

    #[test]
    fn find_activity_overlaps() {
        let mut source1 = File::open("tests/data/gomodel:66a3e0bb00001342.json").unwrap();
        let model1 = make_gocam_model(&mut source1).unwrap();
        assert_eq!(model1.id(), "gomodel:66a3e0bb00001342");

        assert_eq!(model1.node_iterator().count(), 29);

        let mut source2 = File::open("tests/data/gomodel:665912ed00000015.json").unwrap();
        let model2 = make_gocam_model(&mut source2).unwrap();
        assert_eq!(model2.id(), "gomodel:665912ed00000015");

        assert_eq!(model2.node_iterator().count(), 25);

        let mut source3 = File::open("tests/data/gomodel:678073a900003175.json").unwrap();
        let model3 = make_gocam_model(&mut source3).unwrap();
        assert_eq!(model3.id(), "gomodel:678073a900003175");

        assert_eq!(model3.node_iterator().count(), 13);

        let overlaps = GoCamModel::find_activity_overlaps(&[&model1, &model2, &model3]);

        assert_eq!(overlaps.len(), 1);

        let overlap = &overlaps[0];

        assert_eq!(overlap.node_description,
                   "homoserine O-acetyltransferase activity [enabled by] met6 Spom");
        assert_eq!(overlap.model_ids.len(), 2);
    }
}
