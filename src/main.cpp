#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct Options {
    double alpha = 1.0;
    std::size_t steps = 100; // Total time steps n (produces n+1 points)
    std::optional<std::uint64_t> seed;
    bool plot = true;
    double plot_width_cm = 10.0;
    bool show_tree = true;
    bool show_mst = true;
    bool print_lengths = true;
    double dot_size_pt = 1.0;
    double tree_line_width_pt = 0.25;
    double mst_line_width_pt = 0.4;
    std::string output_path;
};

struct Point {
    double x = 0.0;
    double y = 0.0;
    std::size_t parent = 0; // Index of the parent point
};

struct MstResult {
    double length = 0.0;
    std::vector<std::pair<std::size_t, std::size_t>> edges;
};

namespace {

constexpr double kInfinity = 1e300;
constexpr double kMinimumExtent = 1e-6;

void print_usage(const std::string &program) {
    std::cerr << "Usage: " << program << " --alpha <alpha> [options]\n"
              << "\n"
              << "Simulate Gaussian clustering with perturbations of variance t^{-alpha}.\n"
              << "\n"
              << "Required:\n"
              << "  --alpha <double>         Positive alpha parameter for variance decay.\n"
              << "\n"
              << "Options:\n"
              << "  --steps <int>            Number of steps n (default: 100, produces n+1 points).\n"
              << "  --seed <int>             Optional RNG seed for reproducibility.\n"
              << "  --[no-]plot              Enable/disable TikZ output (default: on).\n"
              << "  --width <double>         Target plot width in cm (default: 10.0).\n"
              << "  --[no-]tree              Show construction tree in TikZ (default: on).\n"
              << "  --[no-]mst               Show minimum spanning tree in TikZ (default: on).\n"
              << "  --[no-]lengths           Print lengths of the construction tree and MST (default: on).\n"
              << "  --dot-size <double>      Radius of plotted points in pt (default: 1.0).\n"
              << "  --tree-width <double>    Line width of construction tree in pt (default: 0.25).\n"
              << "  --mst-width <double>     Line width of MST in pt (default: 0.4).\n"
              << "  --output <path>          Write TikZ to a file instead of stdout.\n"
              << "  --help                   Show this message.\n";
}

Options parse_arguments(int argc, char **argv) {
    Options opts;
    bool alpha_set = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--alpha") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--alpha requires a value");
            }
            opts.alpha = std::stod(argv[++i]);
            alpha_set = true;
        } else if (arg == "--steps") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--steps requires a value");
            }
            long long value = std::stoll(argv[++i]);
            if (value < 0) {
                throw std::invalid_argument("--steps must be non-negative");
            }
            opts.steps = static_cast<std::size_t>(value);
        } else if (arg == "--width") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--width requires a value");
            }
            opts.plot_width_cm = std::stod(argv[++i]);
            if (!(opts.plot_width_cm > 0.0)) {
                throw std::invalid_argument("--width must be positive");
            }
        } else if (arg == "--seed") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--seed requires a value");
            }
            long long value = std::stoll(argv[++i]);
            if (value < 0) {
                throw std::invalid_argument("--seed must be non-negative");
            }
            opts.seed = static_cast<std::uint64_t>(value);
        } else if (arg == "--plot") {
            opts.plot = true;
        } else if (arg == "--no-plot") {
            opts.plot = false;
        } else if (arg == "--tree") {
            opts.show_tree = true;
        } else if (arg == "--no-tree") {
            opts.show_tree = false;
        } else if (arg == "--mst") {
            opts.show_mst = true;
        } else if (arg == "--no-mst") {
            opts.show_mst = false;
        } else if (arg == "--lengths") {
            opts.print_lengths = true;
        } else if (arg == "--no-lengths") {
            opts.print_lengths = false;
        } else if (arg == "--dot-size") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--dot-size requires a value");
            }
            opts.dot_size_pt = std::stod(argv[++i]);
            if (!(opts.dot_size_pt > 0.0)) {
                throw std::invalid_argument("--dot-size must be positive");
            }
        } else if (arg == "--tree-width") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--tree-width requires a value");
            }
            opts.tree_line_width_pt = std::stod(argv[++i]);
            if (!(opts.tree_line_width_pt > 0.0)) {
                throw std::invalid_argument("--tree-width must be positive");
            }
        } else if (arg == "--mst-width") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--mst-width requires a value");
            }
            opts.mst_line_width_pt = std::stod(argv[++i]);
            if (!(opts.mst_line_width_pt > 0.0)) {
                throw std::invalid_argument("--mst-width must be positive");
            }
        } else if (arg == "--output") {
            if (i + 1 >= argc) {
                throw std::invalid_argument("--output requires a value");
            }
            opts.output_path = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(EXIT_SUCCESS);
        } else {
            std::ostringstream oss;
            oss << "Unknown argument: " << arg;
            throw std::invalid_argument(oss.str());
        }
    }

    if (!alpha_set) {
        throw std::invalid_argument("--alpha must be provided and positive");
    }
    if (!(opts.alpha > 0.0)) {
        throw std::invalid_argument("--alpha must be positive");
    }
    return opts;
}

double euclidean_distance(const Point &a, const Point &b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<Point> generate_points(const Options &opts) {
    std::mt19937_64 rng;
    if (opts.seed.has_value()) {
        rng.seed(opts.seed.value());
    } else {
        std::random_device rd;
        rng.seed(rd());
    }

    std::vector<Point> points;
    points.reserve(opts.steps + 1);
    points.push_back(Point{0.0, 0.0, 0});

    for (std::size_t t = 1; t <= opts.steps; ++t) {
        std::uniform_int_distribution<std::size_t> parent_dist(0, t - 1);
        std::size_t parent = parent_dist(rng);

        double variance = std::pow(static_cast<double>(t), -opts.alpha);
        double stddev = std::sqrt(variance);
        std::normal_distribution<double> normal(0.0, stddev);

        double dx = normal(rng);
        double dy = normal(rng);

        const Point &parent_point = points[parent];
        Point p{parent_point.x + dx, parent_point.y + dy, parent};
        points.push_back(p);
    }

    return points;
}

double symmetric_extent(const std::vector<Point> &points) {
    double max_abs = 0.0;
    for (const auto &p : points) {
        max_abs = std::max({max_abs, std::fabs(p.x), std::fabs(p.y)});
    }
    return std::max(max_abs, kMinimumExtent);
}

double construction_tree_length(const std::vector<Point> &points) {
    double length = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        length += euclidean_distance(points[i], points[points[i].parent]);
    }
    return length;
}

MstResult compute_mst(const std::vector<Point> &points) {
    const std::size_t n = points.size();
    if (n <= 1) {
        return {};
    }

    std::vector<double> key(n, kInfinity);
    std::vector<std::size_t> parent(n, 0);
    std::vector<bool> in_mst(n, false);
    key[0] = 0.0;

    double total_length = 0.0;
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.reserve(n - 1);

    for (std::size_t count = 0; count < n; ++count) {
        // Pick minimum key vertex not yet included
        double min_key = kInfinity;
        std::size_t u = 0;
        for (std::size_t v = 0; v < n; ++v) {
            if (!in_mst[v] && key[v] < min_key) {
                min_key = key[v];
                u = v;
            }
        }

        in_mst[u] = true;
        if (u != 0) {
            edges.emplace_back(u, parent[u]);
            total_length += key[u];
        }

        for (std::size_t v = 0; v < n; ++v) {
            if (in_mst[v] || v == u) {
                continue;
            }
            double dist = euclidean_distance(points[u], points[v]);
            if (dist < key[v]) {
                key[v] = dist;
                parent[v] = u;
            }
        }
    }

    return {total_length, edges};
}

double point_set_diameter(const std::vector<Point> &points) {
    if (points.size() <= 1) {
        return 0.0;
    }

    double max_dist = 0.0;
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        for (std::size_t j = i + 1; j < points.size(); ++j) {
            max_dist = std::max(max_dist, euclidean_distance(points[i], points[j]));
        }
    }
    return max_dist;
}

std::string render_tikz(const std::vector<Point> &points,
                        const MstResult *mst,
                        bool show_tree,
                        bool show_mst,
                        double plot_width_cm,
                        double dot_size_pt,
                        double tree_line_width_pt,
                        double mst_line_width_pt) {
    std::ostringstream tikz;
    tikz << std::fixed << std::setprecision(5);

    double extent = symmetric_extent(points);
    double box_width = 2.0 * extent;
    double scale_cm = plot_width_cm / box_width;

    tikz << "\\begin{tikzpicture}[x=" << scale_cm << "cm,y=" << scale_cm << "cm]\n";
    tikz << "  \\path[use as bounding box] (" << -extent << "," << -extent << ") rectangle (" << extent << ","
         << extent << ");\n";

    if (show_tree) {
        for (std::size_t i = 1; i < points.size(); ++i) {
            const Point &child = points[i];
            const Point &parent = points[child.parent];
            tikz << "  \\draw[black!70, line width=" << tree_line_width_pt << "pt] (" << parent.x << ","
                 << parent.y << ") -- (" << child.x << "," << child.y << ");\n";
        }
    }

    if (show_mst && mst != nullptr) {
        for (const auto &edge : mst->edges) {
            const Point &a = points[edge.first];
            const Point &b = points[edge.second];
            tikz << "  \\draw[blue, line width=" << mst_line_width_pt << "pt] (" << a.x << "," << a.y
                 << ") -- (" << b.x << "," << b.y << ");\n";
        }
    }

    for (const auto &p : points) {
        tikz << "  \\filldraw[black] (" << p.x << "," << p.y << ") circle (" << dot_size_pt << "pt);\n";
    }

    tikz << "\\end{tikzpicture}\n";
    return tikz.str();
}

} // namespace

int main(int argc, char **argv) {
    Options opts;
    try {
        opts = parse_arguments(argc, argv);
    } catch (const std::exception &ex) {
        std::cerr << ex.what() << "\n\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    std::vector<Point> points = generate_points(opts);
    double tree_length = construction_tree_length(points);
    double diameter = point_set_diameter(points);

    MstResult mst;
    bool need_mst = opts.show_mst || opts.print_lengths;
    if (need_mst) {
        mst = compute_mst(points);
    }

    if (opts.print_lengths) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Construction tree length: " << tree_length << "\n";
        std::cout << "MST length: " << mst.length << "\n";
        std::cout << "Point set diameter: " << diameter << "\n";
    }

    if (opts.plot) {
        const MstResult *mst_ptr = (opts.show_mst || opts.print_lengths) ? &mst : nullptr;
        std::string tikz = render_tikz(points,
                                       mst_ptr,
                                       opts.show_tree,
                                       opts.show_mst,
                                       opts.plot_width_cm,
                                       opts.dot_size_pt,
                                       opts.tree_line_width_pt,
                                       opts.mst_line_width_pt);

        if (!opts.output_path.empty()) {
            std::ofstream out(opts.output_path);
            if (!out) {
                std::cerr << "Failed to open output file: " << opts.output_path << "\n";
                return EXIT_FAILURE;
            }
            out << tikz;
        } else {
            std::cout << tikz;
        }
    }

    return EXIT_SUCCESS;
}
