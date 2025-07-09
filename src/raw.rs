//! The module contains code for reading the low-level JSON
//! representation of [GO-CAM models](https://geneontology.org/docs/gocam-overview).
//!
//! [GoCamRawModel] closely matches the JSON data, containing Fact,
//! Individual and Annotation structs.
//!
//! ## Example
//! ```
//! use std::fs::File;
//! use pombase_gocam::raw::gocam_parse_raw;
//!
//! let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
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
//! ```

use std::{collections::{BTreeMap, BTreeSet, HashMap, HashSet}, fmt::{self, Display}, io::{BufReader, Read}};

use anyhow::Result;

use crate::{ModelId, REL_NAMES};

pub type FactId = String;
pub type IndividualId = String;
pub type PropertyId = String;
pub type AnnotationKey = String;
pub type AnnotationValue = String;

/// An item from the `annotations` collection of the model, a fact or
/// an Individual.
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, Hash)]
pub struct Annotation {
    /// A key like "date" or "contributor"
    pub key: AnnotationKey,
    /// The value, like "2025-03-07"
    pub value: AnnotationValue,
    #[serde(skip_serializing_if="Option::is_none")]
    #[serde(rename = "value-type")]
    pub value_type: Option<String>,
}

/// A fact/relation conecting two individuals in the model
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, Hash)]
pub struct Fact {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub annotations: Vec<Annotation>,
    /// The ID of the subject Individual
    pub subject: IndividualId,
    /// The ID of the object Individual
    pub object: IndividualId,
    /// The relation ID (from [REL_NAMES])
    pub property: PropertyId,
    /// The retation name
    #[serde(rename = "property-label")]
    pub property_label: String,
}

impl Fact {
    pub fn id(&self) -> FactId {
        format!("{}-{}-{}", self.subject, self.property, self.object)
    }
}
/// An ID and a label.  Used in the `filler` field.`
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, Hash)]
pub struct FillerType {
    pub id: Option<String>,
    pub label: Option<String>,
    #[serde(rename = "type")]
    pub type_string: String,
}

/// An ID and a label.  Used in the `type` and `root-type` fields.
#[derive(Deserialize, Serialize, Debug, Clone, PartialEq, Eq, Hash)]
pub struct IndividualType {
    pub id: Option<String>,
    pub label: Option<String>,
    #[serde(rename = "type")]
    pub type_string: String,

    pub filler: Option<FillerType>
}

const MOLECULAR_FUNCTION_ID: &str = "GO:0003674";
const PROTEIN_CONTAINING_COMPLEX_ID: &str = "GO:0032991";
const CHEBI_PROTEIN_ID: &str = "CHEBI:36080";
const CHEBI_CHEMICAL_ENTITY_ID: &str = "CHEBI:24431";
const SO_MRNA_ID: &str = "SO:0000234";

impl IndividualType {
    /// Return the ID or if the ID is None return "UNKNOWN_ID"
    pub fn id(&self) -> &str {
        self.id.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_ID")
    }

    /// Return the label or if the label is None return "UNKNOWN_ID"
    pub fn label(&self) -> &str {
        self.label.as_ref().map(|s| s.as_str()).unwrap_or("UNKNOWN_LABEL")
    }

    // Return the label, if set (not None).  Otherise return the ID.
    // If the label and ID are both unset, return "UNKNOWN"
    pub fn label_or_id(&self) -> &str {
        self.label.as_ref().map(|s| s.as_str())
            .or(self.id.as_ref().map(|s| s.as_str()))
            .unwrap_or("UNKNOWN")
    }
}

impl Display for IndividualType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ({})", self.label(), self.id())?;
        Ok(())
    }
}

/// A node in the raw GO-CAM model
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
    /// Return true if the term_id is in the root_terms of this
    /// Individual
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

    /// Return true if this Individual is an activity, by checking for
    /// "molecular_function" in the root_terms
    pub fn individual_is_activity(&self, model: &GoCamRawModel) -> bool {
        if let Some(individual_type) = self.get_individual_type() {
            // See: https://github.com/pombase/pombase-chado/issues/1262#issuecomment-2708083647
            if individual_type.label().starts_with("obsolete ") &&
                individual_type.id().starts_with("GO:") {
                    return model.facts_of_subject(individual_type.id()).len() == 0;
                }
        }

        self.has_root_term(MOLECULAR_FUNCTION_ID)
    }

    pub fn individual_is_complex(&self) -> bool {
        self.has_root_term(PROTEIN_CONTAINING_COMPLEX_ID)
    }

    pub fn individual_is_modified_protein(&self) -> bool {
        let Some(individual_type) = self.get_individual_type()
        else {
            return false;
        };

        if let Some(ref id) = individual_type.id {
            return id.starts_with("PR:") && id != "PR:000000001";
        }

        return false;
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
            if id.starts_with("CHEBI:") || id.starts_with("SO:") {
                return true;
            }
        }

        false
    }

    pub fn individual_is_unknown_mrna(&self) -> bool {
        let Some(individual_type) = self.get_individual_type()
        else {
            return false;
        };

        if let Some(ref id) = individual_type.id {
            return id == SO_MRNA_ID;
        }

        false
    }

    pub fn individual_is_mrna(&self) -> bool {
        if let Some(individual_type) = self.get_individual_type() {
           return is_mrna_id(individual_type.id());
        }
        false
    }

    /// Return the first element of the types collection
    pub fn get_individual_type(&self) -> Option<&IndividualType> {
        self.types.get(0)
    }

    /// Return true if the Individual is "protein" from ChEBI, but not
    /// a specific protein
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

pub(crate) fn is_mrna_id(id: &str) -> bool {
    if let Some(no_suffix) = id.strip_suffix(|c: char| c.is_numeric()) {
        return no_suffix.ends_with('.')
    }

    false
}

#[derive(Deserialize, Serialize, Debug, Clone)]
struct SerdeModel {
    annotations: Vec<Annotation>,
    id: ModelId,
    facts: Vec<Fact>,
    individuals: Vec<Individual>,
}

/// A container for the GO-CAM JSON format, containing annotations,
/// facts and individuals
///
/// ## Example
/// ```
/// use std::fs::File;
/// use pombase_gocam::raw::gocam_parse_raw;
///
/// let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
/// let raw_model = gocam_parse_raw(&mut source).unwrap();
/// assert_eq!(raw_model.id(), "gomodel:66187e4700001744");
///
/// for fact in raw_model.facts() {
///     let subject_id = &fact.subject;
///     println!("subject_id: {}", subject_id);
///     let subject_individual = raw_model.get_individual(subject_id);
///     let first_type = &subject_individual.types[0];
///     if let Some(ref label) = first_type.label {
///         println!("type label: {}", label);
///     }
/// }
/// ```
#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct GoCamRawModel {
    _annotations: Vec<Annotation>,
    _id: ModelId,
    _facts: BTreeMap<FactId, Fact>,
    _individuals: BTreeMap<IndividualId, Individual>,

    _facts_by_object: HashMap<IndividualId, HashSet<FactId>>,
    _facts_by_subject: HashMap<IndividualId, HashSet<FactId>>,
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

    /// An iterator over the top-level Annotations of this model
    pub fn annotations(&self) -> Box<dyn Iterator<Item = &Annotation> + '_> {
        Box::new(self._annotations.iter())
    }

    pub fn title(&self) -> String {
        annotation_values(self.annotations(), "title")
            .join(",")
    }

    /// The date the model last changed
    pub fn date(&self) -> String {
        annotation_values(self.annotations(), "date")
            .join(",")
    }

    /// The taxon ID in the format "NCBITaxon:4896" or possible a
    /// comma separated list like: "NCBITaxon:4896,NCBITaxon:559292"
    pub fn taxon(&self) -> String {
        annotation_values(self.annotations(), "https://w3id.org/biolink/vocab/in_taxon")
            .join(",")
    }

    #[allow(rustdoc::bare_urls)]
    /// A set of contributor ORCIDs like:
    /// "https://orcid.org/0000-0001-6330-7526"
    pub fn contributors(&self) -> BTreeSet<String> {
        BTreeSet::from_iter(annotation_values(self.annotations(), "contributor"))
    }

    /// An iterator over the facts
    pub fn facts(&self) -> Box<dyn Iterator<Item = &Fact> + '_>  {
        Box::new(self._facts.values())
    }

    /// An iterator over the individuals
    pub fn individuals(&self) -> Box<dyn Iterator<Item = &Individual> + '_> {
        Box::new(self._individuals.values())
    }

    /// Return the subject &Individual of a given Fact
    pub fn fact_subject(&self, fact: &Fact) -> &Individual {
        self.get_individual(&fact.subject)
    }

    /// Return the object &Individual of a given Fact
    pub fn fact_object(&self, fact: &Fact) -> &Individual {
        self.get_individual(&fact.object)
    }

    /// Return a copy of the subject Facts of an Individual - every
    /// Fact where the subject Individual has the given subject_id
    pub fn facts_of_subject(&self, subject_id: &str)
        -> HashSet<Fact>
    {
        if let Some(fact_ids) = self._facts_by_subject.get(subject_id) {
            fact_ids.iter().filter_map(|fact_id| self._facts.get(fact_id))
                    .cloned()
                    .collect()
        } else {
            HashSet::new()
        }
    }

    /// Return a copy of the object Facts of an Individual - every
    /// Fact where the object Individual has the given object_id
    pub fn facts_of_object(&self, object_id: &str)
        -> HashSet<Fact>
    {
        if let Some(fact_ids) = self._facts_by_object.get(object_id) {
            fact_ids.iter().filter_map(|fact_id| self._facts.get(fact_id))
                    .cloned()
                    .collect()
        } else {
            HashSet::new()
        }
    }

    /// Return the Individual with the given an IndividualId.  IDs
    /// have the format "gomodel:67c10cc400002026/67c10cc400002109",
    pub fn get_individual(&self, individual_id: &str)
        -> &Individual
    {
        self._individuals.get(individual_id)
            .expect(&format!("can't find individual: {}",
                             individual_id))
    }
}

fn gocam_parse_raw_helper(source: &mut dyn Read) -> Result<SerdeModel> {
    let reader = BufReader::new(source);

    let mut raw_model: SerdeModel = serde_json::from_reader(reader)?;

    for individual in &mut raw_model.individuals {
        for type_item in &mut individual.types {
            if let Some(ref mut filler) = type_item.filler.take() {
                type_item.id = filler.id.clone();
                type_item.label = filler.label.clone();
                type_item.type_string = "class".into();
            }
        }
    }

    Ok(raw_model)
}


/// Parses a GO-CAM model from a stream in a raw representation of
/// Individuals and Facts
///
/// ## Example
///
/// ```
/// use std::fs::File;
/// use pombase_gocam::raw::gocam_parse_raw;
///
/// let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
/// let model = gocam_parse_raw(&mut source).unwrap();
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
pub fn gocam_parse_raw(source: &mut dyn Read) -> Result<GoCamRawModel> {
    let raw_model = gocam_parse_raw_helper(source)?;

    let mut fact_map = BTreeMap::new();
    let mut individual_map = BTreeMap::new();

    let mut _facts_by_subject = HashMap::new();
    let mut _facts_by_object = HashMap::new();

    for mut fact in raw_model.facts.into_iter() {
        if let Some(&rel_name) = REL_NAMES.get(&fact.property) {
            fact.property_label = rel_name.to_owned();
        }

        _facts_by_subject.entry(fact.subject.clone())
                         .or_insert_with(HashSet::new)
                         .insert(fact.id());
        _facts_by_object.entry(fact.object.clone())
                        .or_insert_with(HashSet::new)
                        .insert(fact.id());

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

        _facts_by_subject,
        _facts_by_object,
    })
}


#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn parse_raw_test() {
        let mut source = File::open("tests/data/gomodel_66187e4700001744.json").unwrap();
        let model = gocam_parse_raw(&mut source).unwrap();
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
    fn parse_raw_test_with_complement() {
        // parse file containing:
        //   "type": "complement"
        let mut source = File::open("tests/data/gomodel_67369e7600002505.json").unwrap();
        let model = gocam_parse_raw(&mut source).unwrap();
        assert_eq!(model.id(), "gomodel:67369e7600002505");
        assert_eq!(model.individuals().count(), 198);

        let mut test_individual = None;

        for individual in model.individuals() {
            if individual.id == "gomodel:67369e7600002505/67369e7600002681" {
                test_individual = Some(individual);
            }
        }

        let Some(test_individual) = test_individual
        else {
            panic!();
        };

        let individual_type = test_individual.get_individual_type().unwrap();

        assert_eq!(individual_type.id(), "GO:0004674");
        assert_eq!(individual_type.label(), "protein serine/threonine kinase activity");
    }
}
