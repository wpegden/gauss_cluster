# gauss_cluster

A small C++17 simulator that builds a growing Gaussian cluster on the plane, records its construction tree, and reports the length of both the construction tree and a minimum spanning tree (MST). The program can also emit TikZ code to plot the points along with optional overlays for the construction tree and MST.

## Building

The repository includes a simple `Makefile` that produces a `gauss_cluster` binary in the repository root.

```bash
make        # builds gauss_cluster
make clean  # removes the binary
```

If you prefer to compile manually:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -o gauss_cluster src/main.cpp
```

## Usage

Run the simulator by providing an `\alpha > 0` value and the number of time steps `n` (the simulation produces `n+1` points, starting from the origin):

```bash
./gauss_cluster --alpha 1.2 --steps 200 [options]
```

### Options

- `--alpha <double>`: Required positive parameter controlling the variance decay.
- `--steps <int>`: Number of time steps `n` (default: `100`).
- `--seed <int>`: Optional non-negative RNG seed for reproducibility.
- `--[no-]plot`: Enable/disable TikZ output (default: on).
- `--[no-]tree`: Include/exclude the construction tree overlay in TikZ output (default: on).
- `--[no-]mst`: Include/exclude the MST overlay in TikZ output (default: on).
- `--[no-]lengths`: Print lengths for the construction tree and MST (default: on).
- `--output <path>`: Write TikZ output to the given file instead of stdout.
- `--help`: Show usage information.

### Output

- **Lengths**: When enabled (default), the program prints the total Euclidean length of the construction tree and the MST.
- **TikZ**: When plotting is enabled (default), TikZ code is printed to stdout or the path supplied with `--output`. Points are black dots, construction-tree edges are thin black lines, and MST edges are blue lines. The origin is centered and the aspect ratio is square.

### Example

Produce a 250-point simulation with both overlays drawn and TikZ saved to `cluster.tex`:

```bash
./gauss_cluster --alpha 0.8 --steps 249 --seed 42 --output cluster.tex
```

To generate TikZ without the MST overlay and suppress length reporting:

```bash
./gauss_cluster --alpha 1.0 --steps 150 --no-mst --no-lengths
```
