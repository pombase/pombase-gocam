use std::{collections::BTreeMap, io::{BufReader, Read}};

extern crate serde_json;
#[macro_use] extern crate serde_derive;
use phf::phf_map;

use anyhow::Result;

static REL_NAMES: phf::Map<&'static str, &'static str> = phf_map! {
    "BFO:0000050" => "part of",
    "BFO:0000051" => "has part",
    "RO:0002233" => "has input",
    "RO:0002234" => "has output",
    "RO:0002333" => "enabled by",
    "RO:0002413" => "directly provides input for",
    "BFO:0000066" => "occurs in",
    "RO:0001025" => "located in",
    "RO:0002092" => "happens during",
    "RO:0002407" => "indirectly activates",
    "RO:0002409" => "indirectly inhibits",
    "RO:0002629" => "directly positively regulates",
    "RO:0002630" => "directly negatively regulates",
    "RO:0002304" => "causally upstream of, positive effect",
    "RO:0002305" => "causally upstream of, negative effect",
    "RO:0012005" => "is small molecule activator of",
    "RO:0012006" => "is small molecule inhibitor of",
    "RO:0012009" => "constitutively upstream of",
};

type ModelId = String;
type FactId = String;
type IndividualId = String;
type PropertyId = String;
type AnnotationKey = String;
type AnnotationValue = String;

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

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Individual {
    #[serde(skip_serializing_if="Vec::is_empty", default)]
    pub annotations: Vec<Annotation>,
    pub id: IndividualId,
    #[serde(rename = "type")]
    pub types: Vec<IndividualType>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct RawModel {
    pub annotations: Vec<Annotation>,
    pub id: ModelId,
    pub facts: Vec<Fact>,
    pub individuals: Vec<Individual>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct GoCamModel {
    _annotations: Vec<Annotation>,
    _id: ModelId,
    _facts: BTreeMap<FactId, Fact>,
    _individuals: BTreeMap<IndividualId, Individual>,
}

impl GoCamModel {
    pub fn id(&self) -> &ModelId {
        &self._id
    }

    pub fn annotations(&self) -> Box<dyn Iterator<Item = &Annotation> + '_> {
        Box::new(self._annotations.iter())
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

    pub fn get_individual(&self, individual_id: &IndividualId)
        -> &Individual
    {
        self._individuals.get(individual_id)
            .expect(&format!("can't find individual: {}",
                            individual_id))
    }
}

pub fn parse_raw(source: &mut dyn Read) -> Result<RawModel> {
    let reader = BufReader::new(source);

    let raw_model: RawModel = serde_json::from_reader(reader)?;

    Ok(raw_model)
}


/// Parses a GO-CAM model from a stream
///
/// # Example:
///
/// ```
/// use std::fs::File;
/// use pombase_gocam::parse;
///
/// let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
/// let model = parse(&mut source).unwrap();
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
pub fn parse(source: &mut dyn Read) -> Result<GoCamModel> {
    let raw_model = parse_raw(source)?;

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

    Ok(GoCamModel {
        _annotations: raw_model.annotations,
        _id: raw_model.id,
        _facts: fact_map,
        _individuals: individual_map,
    })
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn parse_test() {
        let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
        let model = parse(&mut source).unwrap();
        assert!(model.id() == "gomodel:66187e4700001744");
        assert!(model.facts().count() == 42);
        assert!(model.individuals().count() == 82);

        let first_fact = model.facts().next().unwrap();
        assert_eq!(first_fact.property, "BFO:0000050");
        assert_eq!(first_fact.property_label, "part of");
        assert_eq!(first_fact.id(), "gomodel:66187e4700001744/66187e4700001758-BFO:0000050-gomodel:66187e4700001744/66187e4700001760");

        let individual1_id = &first_fact.object;

        let lookup_individual1 = model.get_individual(individual1_id);

        assert_eq!(*individual1_id, lookup_individual1.id);
    }
}
