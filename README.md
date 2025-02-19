# pombase-gocam

[![build status](https://github.com/pombase/pombase-gocam/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/pombase/pombase-gocam/actions)

Code for parsing [GO-CAM](https://geneontology.org/docs/gocam-overview)
JSON format model files.

# Example

```rust
use std::fs::File;
use pombase_gocam::parse;
let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
let model = parse(&mut source).unwrap();
assert!(model.id() == "gomodel:66187e4700001744");

for fact in model.facts() {
  let subject_id = &fact.subject;
  println!("subject_id: {}", subject_id);
  let subject_individual = model.get_individual(subject_id).unwrap();
  let first_type = &subject_individual.types[0];
  if let Some(ref label) = first_type.label {
    println!("first_type label: {}", label);
  }
}
```

# PomBase

The library was developed by the [PomBase](https://www.pombase.org/) project.
