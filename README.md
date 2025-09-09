# pombase-gocam

[![build status](https://github.com/pombase/pombase-gocam/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/pombase/pombase-gocam/actions)
[![Crates.io](https://img.shields.io/crates/v/pombase-gocam.svg)](https://crates.io/crates/pombase-gocam)
[![Documentation](https://docs.rs/pombase-gocam/badge.svg)](https://docs.rs/pombase-gocam)
[![GitHub](https://img.shields.io/badge/GitHub-pombase--gocam-blue)](https://github.com/pombase/pombase-gocam)

Code for parsing and processing [GO-CAM](https://geneontology.org/docs/gocam-overview)
JSON format model files.

There is a low level representation which closely matches the JSON
data: [GoCamRawModel](https://docs.rs/pombase-gocam/latest/pombase_gocam/struct.GoCamRawModel.html)
(containing Fact, Individual and Annotation structs).

And a higher level representation, [GoCamModel](https://docs.rs/pombase-gocam/latest/pombase_gocam/struct.GoCamModel.html),
implemented as a graph of nodes (activities, chemical, complexes etc.)
and edges (mostly causal relations).

The high level representation is similar to the
[GO CAM Data Model - gocam-py](https://github.com/geneontology/gocam-py)

See the [documentation on docs.rs](https://docs.rs/pombase-gocam/latest/pombase_gocam/)
for usage.

# Example

```bash
curl -L https://live-go-cam.geneontology.io/product/json/low-level/665912ed00002626.json |
  jq . > gomodel_665912ed00002626.json

```

```rust
use std::fs::File;
use pombase_gocam::raw::gocam_parse_raw;

fn main() {
    let mut source = File::open("gomodel_665912ed00002626.json").unwrap();

    // Low level representation:
    let raw_model = gocam_parse_raw(&mut source).unwrap();
    assert_eq!(raw_model.id(), "gomodel:665912ed00002626");

    for fact in raw_model.facts() {
        let subject_id = &fact.subject;
        println!("subject_id: {}", subject_id);
        let subject_individual = raw_model.get_individual(subject_id);
        let first_type = &subject_individual.types[0];
        if let Some(ref label) = first_type.label {
            println!("type label: {}", label);
        }
    }

    // Higher level representation:
    use pombase_gocam::{GoCamModel, GoCamNodeType};
    let model = GoCamModel::new(raw_model);

    for (_, node) in model.node_iterator() {
        println!("node: {}", node);
        if let GoCamNodeType::Activity(ref enabler) = node.node_type {
            println!("enabler ID: {}", enabler.id());
        }
    }
}
```

```
cargo run
```

# Authors

The library was developed by the [PomBase](https://www.pombase.org/) project.
