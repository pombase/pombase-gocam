use std::{collections::HashMap, io::{BufReader, Read}};

extern crate serde_json;
#[macro_use] extern crate serde_derive;

use anyhow::Result;

type ModelId = String;
type FactId = String;
type IndividualId = String;
type PropertyId = String;

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Fact {
    pub subject: IndividualId,
    pub object: IndividualId,
    pub property: PropertyId,
}

impl Fact {
    pub fn id(&self) -> FactId {
        format!("{}-{}-{}", self.subject, self.property, self.object)
    }
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct IndividualType {
    pub id: String,
    pub label: String,
    #[serde(rename = "type")]
    pub type_string: String,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct Individual {
    pub id: IndividualId,
    #[serde(rename = "type")]
    pub types: Vec<IndividualType>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct RawModel {
    id: ModelId,
    facts: Vec<Fact>,
    individuals: Vec<Individual>,
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub struct NoctuaModel {
    _id: ModelId,
    _facts: HashMap<FactId, Fact>,
    _individuals: HashMap<IndividualId, Individual>,
}

impl NoctuaModel {
    pub fn id(&self) -> &ModelId {
        &self._id
    }

    pub fn all_facts(&self) -> Box<dyn Iterator<Item = &Fact> + '_>  {
        Box::new(self._facts.values())
    }

    pub fn all_individuals(&self) -> Box<dyn Iterator<Item = &Individual> + '_> {
        Box::new(self._individuals.values())
    }

    pub fn get_individual(&self, individual_id: &IndividualId)
        -> Option<&Individual>
    {
        self._individuals.get(individual_id)
    }
}

/// Parses a GO-CAM model from a stream
///
/// # Example:
///
/// ```
/// use std::fs::File;
/// use pombase_gocam::parse;
/// let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
/// let model = parse(&mut source).unwrap();
/// assert!(model.id() == "gomodel:66187e4700001744");
/// ```
pub fn parse(source: &mut dyn Read) -> Result<NoctuaModel> {
    let reader = BufReader::new(source);

    let raw_model: RawModel = serde_json::from_reader(reader)?;

    let mut fact_map = HashMap::new();
    let mut individual_map = HashMap::new();

    for fact in raw_model.facts.into_iter() {
        fact_map.insert(fact.id(), fact);
    }

    for individual in raw_model.individuals.into_iter() {
        individual_map.insert(individual.id.clone(), individual);
    }

    Ok(NoctuaModel {
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
        assert!(model.all_facts().count() == 42);
        assert!(model.all_individuals().count() == 82);

        let individual1_id =
            &model.all_facts().next().unwrap().object;

        let lookup_individual1 = model.get_individual(individual1_id).unwrap();

        assert!(*individual1_id == lookup_individual1.id);
    }
}
