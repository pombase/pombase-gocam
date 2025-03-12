# pombase-gocam

[![build status](https://github.com/pombase/pombase-gocam/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/pombase/pombase-gocam/actions)

Code for parsing and processing [GO-CAM](https://geneontology.org/docs/gocam-overview)
JSON format model files.

# Example

```bash
curl -L https://live-go-cam.geneontology.io/product/json/low-level/665912ed00002626.json |
  jq . > gomodel:665912ed00002626.json

```

```rust
use std::fs::File;
use pombase_gocam::parse;

let mut source = File::open("gomodel_665912ed00002626.json").unwrap();
let model = gocam_parse(&mut source).unwrap();

for fact in model.facts() {
  let subject_id = &fact.subject;
  println!("subject_id: {}", subject_id);
  let subject_individual = model.get_individual(subject_id);
  let type = &subject_individual.types[0];
  if let Some(ref label) = type.label {
    println!("type label: {}", label);
  }
}
```

# Authors

The library was developed by the [PomBase](https://www.pombase.org/) project.
