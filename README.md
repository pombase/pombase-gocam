# pombase-gocam

Code for parsing [GO-CAM](https://geneontology.org/docs/gocam-overview)
JSON format model files.

# Example

```rust
use std::fs::File;
use pombase_gocam::parse;
let mut source = File::open("tests/data/gomodel:66187e4700001744.json").unwrap();
let model = parse(&mut source).unwrap();

for fact in model.all_facts() {
  let subject_id = &fact.subject;
  println!("subject_id: {}", subject_id);
  let subject_individual = model.get_individual(subject_id).unwrap();
  let first_type = &subject_individual.types[0];
  println!("first_type label: {}", first_type.label);
}
```

# PomBase

The library was developed by the [PomBase](https://www.pombase.org/) project.
