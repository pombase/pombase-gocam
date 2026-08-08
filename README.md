# pombase-gocam

[![build status](https://github.com/pombase/pombase-gocam/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/pombase/pombase-gocam/actions)
[![Crates.io](https://img.shields.io/crates/v/pombase-gocam.svg)](https://crates.io/crates/pombase-gocam)
[![Documentation](https://docs.rs/pombase-gocam/badge.svg)](https://docs.rs/pombase-gocam)
[![GitHub](https://img.shields.io/badge/GitHub-pombase--gocam-blue)](https://github.com/pombase/pombase-gocam)

Code for parsing and processing [GO-CAM](https://geneontology.org/docs/gocam-overview)
gocam-py YAML format model files.

The main entry points are the high level representation,
[GoCamModel](https://docs.rs/pombase-gocam/latest/pombase_gocam/struct.GoCamModel.html),
which is implemented internally as a graph of nodes (activities,
chemical, complexes etc.) and edges (mostly causal relations).

There is a low level representation which closely matches the
[gocam-py YAML data](https://github.com/geneontology/gocam-py):
[GoCamPyModel](https://docs.rs/pombase-gocam/latest/pombase_gocam/gocam_py/struct.GoCamPyModel.html).

See the [documentation on docs.rs](https://docs.rs/pombase-gocam/latest/pombase_gocam/)
for usage.

# Example

```bash
curl -s -L https://live-go-cam.geneontology.io/product/yaml/go-cam/6994852c00004096.yaml > 6994852c00004096.yaml
```

```rust
use std::fs::File;
use pombase_gocam::gocam_py::gocam_py_parse;
use pombase_gocam::{GoCamModel, GoCamNodeType};

fn main() {
    let mut source = File::open("6994852c00004096.yaml").unwrap();

    // Low level representation:
    let gocam_py_model = gocam_py_parse(&mut source).unwrap();
    assert_eq!(gocam_py_model.id, "gomodel:6994852c00004096");

    // Higher level representation:
    let model = GoCamModel::new_from_gocam_py(gocam_py_model);

    for (_, node) in model.node_iterator() {
        println!("node: {}", node);
        if let GoCamNodeType::Activity(ref activity) = node.node_type {
            println!("enabler ID: {}", activity.enabler.id());
        }
    }
}
```

```
cargo run
```

# Authors

The library was developed by the [PomBase](https://www.pombase.org/) project.
