use std::{collections::HashMap, io::{BufReader, Read}};

use anyhow::Result;
use serde::{Deserialize, Deserializer};

pub fn gocam_py_parse(source: &mut dyn Read) -> Result<GoCamPyModel> {
    let reader = BufReader::new(source);

    let raw_model: GoCamPyModel = serde_json::from_reader(reader)?;

    Ok(raw_model)
}


fn deserialize_null_default<'de, D, T>(deserializer: D) -> Result<T, D::Error>
where
    T: Default + Deserialize<'de>,
    D: Deserializer<'de>,
{
    let opt = Option::deserialize(deserializer)?;
    Ok(opt.unwrap_or_default())
}


type UriOrCurie = String;

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum ModelStateEnum {
    Production,
    Development,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum InformationBiomacromoleculeCategory {
    GeneOrReferenceProtein,
    ProteinIsoform,
    MacromolecularComplex,
    Unknown,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum CausalPredicateEnum {
    CausallyUpstreamOfPositiveEffect,
    CausallyUpstreamOfNegativeEffect,
    CausallyUpstreamOf,
    ImmediatelyCausallyUpstreamOf,
    CausallyUpstreamOfOrWithin,
    CausallyUpstreamOfOrWithinNegativeEffect,
    CausallyUpstreamOfOrWithinPositiveEffect,
    Regulates,
    NegativelyRegulates,
    PositivelyRegulates,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GoCamPyModel {
    pub id: UriOrCurie,
    pub title: String,
    pub taxon: TaxonTermObject,
    pub status: String,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub comments: Vec<String>,
    pub activities: Vec<Activity>,
    pub objects: Vec<Object>,
    pub provenances: Vec<ProvenanceInfo>,
    pub query_index: QueryIndex
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Activity {
    pub id: UriOrCurie,
    pub enabled_by: EnabledByAssociation,
    pub molecular_function: MolecularFunctionAssociation,
    pub occurs_in: CellularAnatomicalEntityAssociation,
    pub part_of: BiologicalProcessAssociation,
    pub has_input: Vec<MoleculeAssociation>,
    pub has_primary_input: MoleculeAssociation,
    pub has_output: Vec<MoleculeAssociation>,
    pub has_primary_output: MoleculeAssociation,
    pub causal_associations: Vec<CausalAssociation>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EvidenceItem {
    pub term: EvidenceTermObject,
    pub reference: PublicationObject,
    pub with_objects: Vec<Object>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Association {
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByAssociation {
    pub term: InformationBiomacromoleculeTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByGeneProductAssociation {
    pub term: GeneProductTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByProteinComplexAssociation {
    pub members: Vec<GeneProductTermObject>,
    pub term: ProteinComplexTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CausalAssociation {
    pub predicate: PredicateTermObject,
    pub downstream_activity: Box<Activity>,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct TermAssociation {
    pub term: TermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MolecularFunctionAssociation {
    pub term: MolecularFunctionTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct BiologicalProcessAssociation {
    pub happens_during: PhaseTermObject,
    pub part_of: BiologicalProcessTermObject,
    pub term: BiologicalProcessTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellularAnatomicalEntityAssociation {
    pub part_of: CellTypeAssociation,
    pub term: CellularAnatomicalEntityTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellTypeAssociation {
    pub part_of: GrossAnatomyAssociation,
    pub term: CellTypeTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GrossAnatomyAssociation {
    pub part_of: Box<GrossAnatomyAssociation>,
    pub term: GrossAnatomicalStructureTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAssociation {
    pub term: MoleculeTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: HashMap<String, EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Object {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct TermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct PublicationObject {
    pub abstract_text: String,
    pub full_text: String,
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EvidenceTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MolecularFunctionTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct BiologicalProcessTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellularAnatomicalEntityTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MoleculeTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellTypeTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GrossAnatomicalStructureTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct PhaseTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct InformationBiomacromoleculeTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GeneProductTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ProteinComplexTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct TaxonTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct PredicateTermObject {
    pub id: UriOrCurie,
    pub label: String,
    pub type_: UriOrCurie,
    pub obsolete: bool
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ProvenanceInfo {
    pub contributor: Vec<String>,
    pub created: String,
    pub date: String,
    pub provided_by: Vec<String>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct QueryIndex {
    pub number_of_activities: isize,
    pub number_of_enabled_by_terms: isize,
    pub number_of_causal_associations: isize,
    pub length_of_longest_causal_association_path: isize,
    pub number_of_strongly_connected_components: isize,
    pub flattened_references: Vec<PublicationObject>,
    pub model_activity_molecular_function_terms: Vec<TermObject>,
    pub model_activity_molecular_function_closure: Vec<TermObject>,
    pub model_activity_occurs_in_terms: Vec<TermObject>,
    pub model_activity_occurs_in_closure: Vec<TermObject>,
    pub model_activity_part_of_terms: Vec<TermObject>,
    pub model_activity_part_of_closure: Vec<TermObject>,
    pub model_activity_has_input_terms: Vec<TermObject>,
    pub model_activity_has_input_closure: Vec<TermObject>,
    pub taxon_closure: Vec<TermObject>,
    pub annoton_terms: Vec<TermObject>
}
