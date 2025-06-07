use std::io::{BufReader, Read};

use anyhow::Result;
use serde::{Deserialize, Deserializer};

pub fn gocam_py_parse(source: &mut dyn Read) -> Result<GoCamPyModel> {
    let reader = BufReader::new(source);

    let raw_model: GoCamPyModel = serde_yaml::from_reader(reader)?;

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
    pub taxon: String,
    pub status: String,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub comments: Vec<String>,
    pub activities: Vec<Activity>,
    pub objects: Vec<Object>,
    pub provenances: Vec<ProvenanceInfo>,
    pub query_index: Option<QueryIndex>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Activity {
    pub id: UriOrCurie,
    pub enabled_by: EnabledByAssociation,
    pub molecular_function: MolecularFunctionAssociation,
    pub occurs_in: Option<CellularAnatomicalEntityAssociation>,
    pub part_of: Option<BiologicalProcessAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub has_input: Vec<MoleculeAssociation>,
    pub has_primary_input: Option<MoleculeAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub has_output: Vec<MoleculeAssociation>,
    pub has_primary_output: Option<MoleculeAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub causal_associations: Vec<CausalAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EvidenceItem {
    pub term: EvidenceTermObject,
    pub reference: Option<PublicationObject>,
#[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub with_objects: Vec<String>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Association {
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByAssociation {
    pub term: String,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByGeneProductAssociation {
    pub term: String,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByProteinComplexAssociation {
    pub members: Vec<String>,
    pub term: ProteinComplexTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CausalAssociation {
    pub predicate: PredicateTermObject,
    pub downstream_activity: String,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct TermAssociation {
    pub term: TermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MolecularFunctionAssociation {
    pub term: String,
#[serde(rename = "type")]
    pub type_: String,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct BiologicalProcessAssociation {
    pub happens_during: Option<PhaseTermObject>,
    pub part_of: Option<BiologicalProcessTermObject>,
    pub term: BiologicalProcessTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellularAnatomicalEntityAssociation {
    pub part_of: Option<CellTypeAssociation>,
    pub term: CellularAnatomicalEntityTermObject,
#[serde(rename = "type")]
    pub type_: String,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellTypeAssociation {
    pub part_of: GrossAnatomyAssociation,
    pub term: CellTypeTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GrossAnatomyAssociation {
    pub part_of: Box<GrossAnatomyAssociation>,
    pub term: GrossAnatomicalStructureTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAssociation {
    pub term: MoleculeTermObject,
#[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Object {
    pub id: UriOrCurie,
    pub label: Option<String>,
#[serde(rename = "type")]
    pub type_: UriOrCurie,
    #[serde(default)]
    pub obsolete: bool
}

pub type TermObject = String;
pub type PublicationObject = String;
pub type EvidenceTermObject = String;

pub type MolecularFunctionTermObject = String;


pub type BiologicalProcessTermObject = String;


pub type CellularAnatomicalEntityTermObject = String;


pub type MoleculeTermObject = String;


pub type CellTypeTermObject = String;


pub type GrossAnatomicalStructureTermObject = String;


pub type PhaseTermObject = String;


pub type InformationBiomacromoleculeTermObject = String;


pub type GeneProductTermObject = String;


pub type ProteinComplexTermObject = String;


pub type TaxonTermObject = String;


pub type PredicateTermObject = String;


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ProvenanceInfo {
    pub contributor: Vec<String>,
    pub created: Option<String>,
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


#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn parse_test() {
        let mut source = File::open("tests/data/6690711d00000331.yaml").unwrap();
        let model = gocam_py_parse(&mut source).unwrap();

        assert_eq!(model.id, "gomodel:6690711d00000331");
        assert_eq!(model.activities.len(), 19);
    }
}
