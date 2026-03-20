To build for browser:
wasm-pack build --target web --out-dir docs/pkg

To serve locally:
python -m http.server 8000

To edit configs:
http://localhost:8080/docs/editor.html 

To run configs:
http://localhost:8080/docs/index.html 

To benchmark
cargo flamegraph --bin bench --release -- .\config\checkerboard.json --steps 5000
