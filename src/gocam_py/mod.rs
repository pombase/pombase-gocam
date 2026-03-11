//! Code for parsing gocam-py YAML models
//! See: <https://github.com/geneontology/gocam-py>

use std::{collections::HashMap, io::{BufReader, Read}};

use serde::{Deserialize, Deserializer};

use crate::GoCamError;

/// Parse a model in gocam-py YAML format.
pub fn gocam_py_parse(source: &mut dyn Read) -> Result<GoCamPyModel, GoCamError> {
    let reader = BufReader::new(source);

    let mut raw_model: GoCamPyModel = serde_yaml::from_reader(reader)?;

    for obj in &raw_model.objects {
        raw_model.objects_by_id.insert(obj.id.clone(), obj.clone());
    }

    let model_title = raw_model.title.replace("\n", " ");

    raw_model.title = trim_whitespace(&model_title);

    Ok(raw_model)
}

fn trim_whitespace(s: &str) -> String {
    let mut new_str = s.trim().to_owned();
    let mut prev = ' '; // The initial value doesn't really matter
    new_str.retain(|ch| {
        let result = ch != ' ' || prev != ' ';
        prev = ch;
        result
    });
    new_str
}

/// A deserialiser for null values that treats them the same as
/// missing values.
fn deserialize_null_default<'de, D, T>(deserializer: D) -> Result<T, D::Error>
where
    T: Default + Deserialize<'de>,
    D: Deserializer<'de>,
{
    let opt = Option::deserialize(deserializer)?;
    Ok(opt.unwrap_or_default())
}

pub type UriOrCurie = String;
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

pub type GoCamPyObjectMap = HashMap<String, Object>;

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GoCamPyModel {
    pub id: UriOrCurie,
    pub title: String,
    pub taxon: TaxonTermObject,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub additional_taxa: Vec<TaxonTermObject>,
    pub status: String,
    #[serde(skip_serializing_if="Option::is_none")]
    pub date_modified: Option<String>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub comments: Vec<String>,
    pub activities: Vec<Activity>,
    pub molecules: Vec<MoleculeNode>,
    pub objects: Vec<Object>,
    #[serde(skip_serializing_if="HashMap::is_empty", default)]
    objects_by_id: HashMap<UriOrCurie, Object>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub query_index: Option<QueryIndex>
}

impl GoCamPyModel {
    pub fn get_object(&self, object_id: &UriOrCurie) -> Option<&Object> {
        self.objects_by_id.get(object_id)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Activity {
    pub id: UriOrCurie,
    pub enabled_by: EnabledByAssociation,
    pub molecular_function: MolecularFunctionAssociation,
    #[serde(skip_serializing_if="Option::is_none")]
    pub occurs_in: Option<CellularAnatomicalEntityAssociation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of: Option<BiologicalProcessAssociation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub happens_during: Option<PhaseAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub has_input: Vec<MoleculeAssociation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub has_primary_input: Option<MoleculeAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub has_output: Vec<MoleculeAssociation>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub has_primary_output: Option<MoleculeAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub causal_associations: Vec<CausalAssociation>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum ModelStateEnum {
    Development,
    Production,
    Delete,
    Review,
    InternalTest,
    Closed,
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
    ProvidesInputFor,
    RemovesInputFor,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum EvidenceCodeEnum {
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum CellularAnatomicalEntityEnum {
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum PhaseEnum {
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MoleculeNode {
    pub id: UriOrCurie,
    pub term: MoleculeTermObject,
    #[serde(skip_serializing_if="Option::is_none")]
    pub located_in: Option<CellularAnatomicalEntityAssociation>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EvidenceItem {
    pub term: EvidenceTermObject,
    #[serde(skip_serializing_if="Option::is_none")]
    pub reference: Option<PublicationObject>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub with_objects: Vec<String>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Association {
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByAssociation {
    pub term: String,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub members: Vec<ProteinComplexMemberAssociation>,  // if type is "EnabledByProteinComplexAssociation"
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ProteinComplexMemberAssociation {
    pub term: String,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    pub provenances: Vec<ProvenanceInfo>,
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub enum GoCamPyEnablerType {
    Complex,
    Gene,
    Chemical,
    ModifiedProtein,
}

impl EnabledByAssociation {
    pub fn enabler_type(&self) -> GoCamPyEnablerType {
        if self.term.starts_with("PR:") {
            return GoCamPyEnablerType::ModifiedProtein;
        }
        if self.term.starts_with("CHEBI:") || &self.term == "SO:0000234" ||
            &self.term == "SO:0000185" {
            return GoCamPyEnablerType::Chemical;
        }
        match self.type_.as_str() {
            "EnabledByGeneProductAssociation" => GoCamPyEnablerType::Gene,
            "EnabledByProteinComplexAssociation" => GoCamPyEnablerType::Complex,
            _ => panic!("{}", self.type_),
        }
    }
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByGeneProductAssociation {
    pub term: String,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct EnabledByProteinComplexAssociation {
    pub members: Vec<String>,
    pub term: ProteinComplexTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CausalAssociation {
    pub predicate: PredicateTermObject,
    pub downstream_activity: UriOrCurie,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct TermAssociation {
    pub term: TermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
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
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of: Option<Box<BiologicalProcessAssociation>>,
    pub term: BiologicalProcessTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct BiologicalProcessPhaseAssociation {
    pub term: PhaseTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct PhaseAssociation {
    pub term: PhaseTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellularAnatomicalEntityAssociation {
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of: Option<CellTypeAssociation>,
    pub term: CellularAnatomicalEntityTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct CellTypeAssociation {
    #[serde(skip_serializing_if="Option::is_none")]
    pub part_of: Option<GrossAnatomyAssociation>,
    pub term: CellTypeTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct GrossAnatomyAssociation {
    pub part_of: Box<GrossAnatomyAssociation>,
    pub term: GrossAnatomicalStructureTermObject,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}


#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct MoleculeAssociation {
    pub molecule: UriOrCurie,
    #[serde(rename = "type")]
    pub type_: String,
    pub evidence: Vec<EvidenceItem>,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provenances: Vec<ProvenanceInfo>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct Object {
    pub id: UriOrCurie,
    #[serde(skip_serializing_if="Option::is_none")]
    pub label: Option<String>,
    #[serde(rename = "type")]
    pub type_: UriOrCurie,
    #[serde(default)]
    pub obsolete: bool
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct ProvenanceInfo {
    pub contributor: Vec<String>,
    #[serde(skip_serializing_if="Option::is_none")]
    pub created: Option<String>,
    pub date: String,
    #[serde(skip_serializing_if="Vec::is_empty", default, deserialize_with = "deserialize_null_default")]
    pub provided_by: Vec<String>
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub struct QueryIndex {
    pub taxon_label: String,
    pub number_of_activities: isize,
    pub number_of_enabled_by_terms: isize,
    pub number_of_causal_associations: isize,
    pub length_of_longest_causal_association_path: isize,
    pub number_of_strongly_connected_components: isize,
    pub flattened_references: Vec<PublicationObject>,
    pub flattened_provided_by: Vec<Object>,
    pub flattened_contributors: Vec<String>,
    pub flattened_evidence_terms: Vec<EvidenceTermObject>,
    pub model_activity_molecular_function_terms: Vec<TermObject>,
    pub model_activity_molecular_function_closure: Vec<TermObject>,
    pub model_activity_molecular_function_rollup: Vec<TermObject>,
    pub model_activity_occurs_in_terms: Vec<TermObject>,
    pub model_activity_occurs_in_closure: Vec<TermObject>,
    pub model_activity_occurs_in_rollup: Vec<TermObject>,
    pub model_activity_enabled_by_terms: Vec<TermObject>,
    pub model_activity_enabled_by_closure: Vec<TermObject>,
    pub model_activity_enabled_by_rollup: Vec<TermObject>,
    pub model_activity_enabled_by_genes: Vec<TermObject>,
    pub model_activity_part_of_terms: Vec<TermObject>,
    pub model_activity_part_of_closure: Vec<TermObject>,
    pub model_activity_part_of_rollup: Vec<TermObject>,
    pub model_activity_has_input_terms: Vec<TermObject>,
    pub model_activity_has_input_closure: Vec<TermObject>,
    pub model_activity_has_input_rollup: Vec<TermObject>,
    pub model_taxon: Vec<TermObject>,
    pub model_taxon_closure: Vec<TermObject>,
    pub model_taxon_rollup: Vec<TermObject>,
    pub annoton_terms: Vec<TermObject>,
    pub start_activities: Vec<Activity>,
    pub end_activities: Vec<Activity>,
    pub intermediate_activities: Vec<Activity>,
    pub singleton_activities: Vec<Activity>,
    pub number_of_start_activities: isize,
    pub number_of_end_activities: isize,
    pub number_of_intermediate_activities: isize,
    pub number_of_singleton_activities: isize
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn parse_test() {
        let mut source = File::open("tests/data/67ae98b500000055.yaml").unwrap();
        let model = gocam_py_parse(&mut source).unwrap();

        assert_eq!(model.id, "gomodel:67ae98b500000055");
        assert_eq!(model.title, "iron import into cell (GO:0033212) / siderophore biosynthetic process (GO:0019290)");
        assert_eq!(model.activities.len(), 15);

        let first_activity = model.activities.first().unwrap();
        let enabled_by = &first_activity.enabled_by;
        let enabled_by_term_id = &enabled_by.term;
        let enabled_by_term_object = model.get_object(enabled_by_term_id).unwrap();
        assert_eq!(enabled_by_term_object.label.as_ref().unwrap(), "protein");
    }
}
