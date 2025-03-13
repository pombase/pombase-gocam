# pombase-gocam

[![build status](https://github.com/pombase/pombase-gocam/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/pombase/pombase-gocam/actions)
[![Crates.io](https://img.shields.io/crates/v/pombase-gocam.svg)](https://crates.io/crates/pombase-gocam)
[![Documentation](https://docs.rs/pombase-gocam/badge.svg)](https://docs.rs/pombase-gocam)
[![GitHub](https://img.shields.io/badge/GitHub-pombase--gocam-blue)](https://github.com/pombase/pombase-gocam)

Code for parsing and processing [GO-CAM](https://geneontology.org/docs/gocam-overview)
JSON format model files.

The main struct is [GoCamModel] representating a graph of nodes
(activities, chemical, complexes etc.) and edges (most causal
relations).

This representation is similar to the
[GO CAM Data Model - gocam-py](https://github.com/geneontology/gocam-py)

See the [documentation](https://docs.rs/pombase-gocam/latest/pombase_gocam/)
for usage.

A lower level representation is available, [GoCamRawModel], which is
closely matches the GO-CAM JSON (with Fact, Individual and Annotation
structs).

# Example

```bash
curl -L https://live-go-cam.geneontology.io/product/json/low-level/665912ed00002626.json |
  jq . > gomodel_665912ed00002626.json

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

// Higher level representation:
use pombase_gocam::GoCamModel;
let model = GoCamModel::new(raw_model);

for node in model.node_iterator() {
    println!("node: {}", node);
}
```

# Authors

The library was developed by the [PomBase](https://www.pombase.org/) project.
