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

fn annotation_values(annotations: Box<dyn Iterator<Item = &Annotation> + '_>,
                     key: &str)
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

    pub fn taxon(&self) -> String {
        annotation_values(self.annotations(),
                          "https://w3id.org/biolink/vocab/in_taxon")
            .join(",")
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

#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn parse_test() {
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
}
