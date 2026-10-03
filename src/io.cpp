// Input and output files. A tiny JSON reader for the input, and writers for
// the mesh (VTK), the report and the recorded steps (JSON).

#include "io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace mesher {
namespace {

// A small JSON reader: enough for the input file, nothing more. This keeps
// the program free of outside libraries.
struct Json {
    enum class Type { null, boolean, number, string, array, object };
    Type type = Type::null;
    double number = 0.0;
    std::string text;
    std::vector<Json> items;
    std::vector<std::pair<std::string, Json>> members;

    const Json* find(const std::string& key) const {
        for (const auto& [name, value] : members) {
            if (name == key) return &value;
        }
        return nullptr;
    }
};

// Reads one JSON value from text, character by character.
class Parser {
public:
    explicit Parser(const std::string& source) : s_(source) {}

    // The whole text must be one value, with nothing after it.
    Json parse_document() {
        Json value = parse_value();
        skip_space();
        if (i_ != s_.size()) fail("unexpected trailing characters");
        return value;
    }

private:
    // Stop with an error that says where in the file the problem is.
    [[noreturn]] void fail(const std::string& message) const {
        throw std::runtime_error("JSON error at byte " + std::to_string(i_) + ": " + message);
    }

    // Move past spaces, tabs and line breaks.
    void skip_space() {
        while (i_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[i_]))) ++i_;
    }

    // Move past character c if it comes next; say whether it did.
    bool consume(char c) {
        skip_space();
        if (i_ < s_.size() && s_[i_] == c) {
            ++i_;
            return true;
        }
        return false;
    }

    // Character c must come next.
    void expect(char c) {
        if (!consume(c)) fail(std::string("expected '") + c + "'");
    }

    // One value: object, array, string, true/false, null or number.
    Json parse_value() {
        skip_space();
        if (i_ >= s_.size()) fail("unexpected end of input");
        Json value;
        const char c = s_[i_];
        if (c == '{') {
            value.type = Json::Type::object;
            ++i_;
            if (consume('}')) return value;
            do {
                skip_space();
                std::string key = parse_string();
                expect(':');
                value.members.emplace_back(std::move(key), parse_value());
            } while (consume(','));
            expect('}');
        } else if (c == '[') {
            value.type = Json::Type::array;
            ++i_;
            if (consume(']')) return value;
            do {
                value.items.push_back(parse_value());
            } while (consume(','));
            expect(']');
        } else if (c == '"') {
            value.type = Json::Type::string;
            value.text = parse_string();
        } else if (s_.compare(i_, 4, "true") == 0 || s_.compare(i_, 5, "false") == 0) {
            value.type = Json::Type::boolean;
            value.number = s_[i_] == 't' ? 1.0 : 0.0;
            i_ += s_[i_] == 't' ? 4 : 5;
        } else if (s_.compare(i_, 4, "null") == 0) {
            i_ += 4;
        } else {
            value.type = Json::Type::number;
            value.number = parse_number();
        }
        return value;
    }

    /// A JSON number: optional minus, integer part, optional fraction and exponent.
    /// strtod alone would also take "+1", ".5", "0x10", "inf" and "nan".
    double parse_number() {
        const std::size_t begin = i_;
        const auto digits = [this] {
            const std::size_t start = i_;
            while (i_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[i_]))) ++i_;
            return i_ > start;
        };
        if (i_ < s_.size() && s_[i_] == '-') ++i_;
        if (i_ < s_.size() && s_[i_] == '0') {
            ++i_;
        } else if (!digits()) {
            i_ = begin;
            fail("expected a value");
        }
        if (i_ < s_.size() && s_[i_] == '.') {
            ++i_;
            if (!digits()) fail("expected digits after the decimal point");
        }
        if (i_ < s_.size() && (s_[i_] == 'e' || s_[i_] == 'E')) {
            ++i_;
            if (i_ < s_.size() && (s_[i_] == '+' || s_[i_] == '-')) ++i_;
            if (!digits()) fail("expected digits in the exponent");
        }
        return std::strtod(s_.substr(begin, i_ - begin).c_str(), nullptr);
    }

    // A quoted string.
    std::string parse_string() {
        if (i_ >= s_.size() || s_[i_] != '"') fail("expected a string");
        ++i_;
        std::string result;
        while (i_ < s_.size() && s_[i_] != '"') {
            if (s_[i_] == '\\') {
                // Only keys are strings in this schema, so an escape just keeps the next character.
                ++i_;
                if (i_ >= s_.size()) break;
            }
            result += s_[i_++];
        }
        if (i_ >= s_.size()) fail("unterminated string");
        ++i_;
        return result;
    }

    const std::string& s_;
    std::size_t i_ = 0;
};

// The value as a number, or an error that names the input key.
double as_number(const Json& value, const std::string& what) {
    if (value.type != Json::Type::number || !std::isfinite(value.number)) {
        throw std::runtime_error(what + " must be a finite number");
    }
    return value.number;
}

// The value as a list of [x, y] points, or an error.
std::vector<Point> as_loop(const Json& value, const std::string& what) {
    if (value.type != Json::Type::array) throw std::runtime_error(what + " must be an array of [x, y] points");
    std::vector<Point> loop;
    for (const Json& item : value.items) {
        if (item.type != Json::Type::array || item.items.size() != 2) {
            throw std::runtime_error(what + " must contain [x, y] pairs");
        }
        loop.push_back({as_number(item.items[0], what), as_number(item.items[1], what)});
    }
    return loop;
}

// Open a file for writing, with enough digits that every double reads back exactly.
std::ofstream open_output(const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    return out;
}

// [x,y]
void write_point(std::ostream& out, const Point& p) { out << '[' << p[0] << ',' << p[1] << ']'; }

// [[a,b,c],...]
void write_triangles(std::ostream& out, const std::vector<Tri>& triangles) {
    out << '[';
    for (std::size_t i = 0; i < triangles.size(); ++i) {
        const Tri& t = triangles[i];
        out << (i ? "," : "") << '[' << t[0] << ',' << t[1] << ',' << t[2] << ']';
    }
    out << ']';
}

// The word for each reason, as used in steps.json and on the web page.
const char* reason_name(Cause::Reason reason) {
    switch (reason) {
    case Cause::Reason::encroached: return "encroached";
    case Cause::Reason::small_angle: return "small_angle";
    case Cause::Reason::too_long: return "too_long";
    case Cause::Reason::centre_encroaches: return "centre_encroaches";
    case Cause::Reason::centre_outside: return "centre_outside";
    case Cause::Reason::centre_on_segment: return "centre_on_segment";
    }
    throw std::logic_error("unknown refinement reason");
}

}  // namespace

// Read and check the input file. Any problem stops the run with a message that names it.
Input read_input(const std::string& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot read " + path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    const std::string source = buffer.str();
    const Json root = Parser(source).parse_document();
    if (root.type != Json::Type::object) throw std::runtime_error("input must be a JSON object");

    // Every key must be known: a misspelt "max_egde" would otherwise run silently with no size limit.
    static const std::set<std::string> known{"outer", "holes", "min_angle_deg", "max_edge"};
    std::set<std::string> seen;
    for (const auto& [name, value] : root.members) {
        if (known.count(name) == 0) {
            throw std::runtime_error("unknown input key \"" + name +
                                     "\"; the keys are outer, holes, min_angle_deg and max_edge");
        }
        if (!seen.insert(name).second) throw std::runtime_error("input key \"" + name + "\" appears twice");
    }

    Input input;
    const Json* outer = root.find("outer");
    if (outer == nullptr) throw std::runtime_error("input needs an \"outer\" loop");
    input.domain.outer = as_loop(*outer, "\"outer\"");
    if (const Json* holes = root.find("holes")) {
        if (holes->type != Json::Type::array) throw std::runtime_error("\"holes\" must be an array of loops");
        for (const Json& hole : holes->items) input.domain.holes.push_back(as_loop(hole, "a hole"));
    }
    const Json* angle = root.find("min_angle_deg");
    if (angle == nullptr) throw std::runtime_error("input needs \"min_angle_deg\"");
    input.settings.min_angle_deg = as_number(*angle, "\"min_angle_deg\"");
    if (!(input.settings.min_angle_deg > 0.0 && input.settings.min_angle_deg < 60.0)) {
        throw std::runtime_error("\"min_angle_deg\" must be between 0 and 60");
    }
    if (const Json* edge = root.find("max_edge")) {
        input.settings.max_edge = as_number(*edge, "\"max_edge\"");
        if (!(input.settings.max_edge > 0.0)) throw std::runtime_error("\"max_edge\" must be positive");
    }
    return input;
}

// The mesh as a VTK file (opens in ParaView), with each triangle's smallest angle.
void write_vtk(const std::string& path, const std::vector<Point>& points, const std::vector<Tri>& triangles) {
    std::ofstream out = open_output(path);
    out << "# vtk DataFile Version 3.0\nmesher triangle mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << points.size() << " double\n";
    for (const Point& p : points) out << p[0] << ' ' << p[1] << " 0\n";
    out << "CELLS " << triangles.size() << ' ' << 4 * triangles.size() << '\n';
    for (const Tri& t : triangles) out << "3 " << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
    out << "CELL_TYPES " << triangles.size() << '\n';
    for (std::size_t i = 0; i < triangles.size(); ++i) out << "5\n";  // VTK_TRIANGLE
    out << "CELL_DATA " << triangles.size() << "\nSCALARS min_angle double 1\nLOOKUP_TABLE default\n";
    for (const Tri& t : triangles) {
        const auto angles = triangle_angles(points, t);
        out << *std::min_element(angles.begin(), angles.end()) << '\n';
    }
}

// report.json, written by hand in a fixed order.
void write_report(const std::string& path, const Report& r) {
    std::ofstream out = open_output(path);
    const Settings& s = r.input.settings;
    const Quality& q = r.quality;
    out << "{\n";
    out << "  \"targets\": {\"min_angle_deg\": " << s.min_angle_deg << ", \"max_edge\": ";
    if (s.max_edge > 0.0) out << s.max_edge; else out << "null";
    out << ", \"met\": " << (r.stats.targets_met() ? "true" : "false") << "},\n";
    out << "  \"points\": " << r.points << ",\n";
    out << "  \"triangles\": " << r.triangles << ",\n";
    out << "  \"boundary_segments\": " << r.boundary_segments << ",\n";
    out << "  \"min_angle_deg\": " << q.min_angle << ",\n";
    out << "  \"max_angle_deg\": " << q.max_angle << ",\n";
    out << "  \"max_radius_edge_ratio\": " << q.max_radius_edge << ",\n";
    out << "  \"max_edge\": " << q.max_edge << ",\n";
    out << "  \"area\": {\"min\": " << q.min_area << ", \"max\": " << q.max_area << ", \"total\": "
        << q.total_area << ", \"domain\": " << r.domain_area << "},\n";
    out << "  \"angle_histogram\": {\"bin_width_deg\": " << kAngleHistogramBinDeg << ", \"counts\": [";
    for (std::size_t i = 0; i < q.angle_histogram.size(); ++i) out << (i ? ", " : "") << q.angle_histogram[i];
    out << "]},\n";
    const RefineStats& st = r.stats;
    out << "  \"refinement\": {\"inserted\": " << st.inserted << ", \"splits\": " << st.splits
        << ", \"refused_splits\": " << st.refused_splits << ", \"rejected_points\": " << st.rejected_points
        << ", \"hit_point_limit\": " << (st.hit_point_limit ? "true" : "false")
        << ", \"triangles_below_min_angle\": " << st.below_min_angle
        << ", \"of_which_at_sharp_input_corners\": " << st.at_sharp_corner
        << ", \"triangles_over_max_edge\": " << st.over_max_edge << "},\n";
    out << "  \"time_ms\": {\"triangulate\": " << r.triangulate_ms << ", \"refine\": " << r.refine_ms << "}\n";
    out << "}\n";
}

// steps.json: the input, steps 2 to 4 change by change, the first full mesh,
// then every step of step 5 with its reason.
void write_steps(const std::string& path, const Replay& r) {
    std::ofstream out = open_output(path);
    const auto write_loop = [&out](const std::vector<Point>& loop) {
        out << '[';
        for (std::size_t i = 0; i < loop.size(); ++i) {
            if (i) out << ',';
            write_point(out, loop[i]);
        }
        out << ']';
    };
    const auto write_delta = [&out](const Delta& d) {
        out << "\"removed\": ";
        write_triangles(out, d.removed);
        out << ", \"added\": ";
        write_triangles(out, d.added);
    };

    const Domain& domain = r.input.domain;
    out << "{\"input\": {\"outer\": ";
    write_loop(domain.outer);
    out << ", \"holes\": [";
    for (std::size_t i = 0; i < domain.holes.size(); ++i) {
        if (i) out << ',';
        write_loop(domain.holes[i]);
    }
    out << "], \"min_angle_deg\": " << r.input.settings.min_angle_deg << ", \"max_edge\": ";
    if (r.input.settings.max_edge > 0.0) out << r.input.settings.max_edge; else out << "null";
    out << "},\n";

    out << "\"points\": ";
    write_loop(r.boundary_points);
    out << ",\n\"build\": {\"enclosing\": ";
    write_loop(std::vector<Point>(r.build.enclosing.begin(), r.build.enclosing.end()));
    out << ",\n\"insertions\": [";
    for (std::size_t i = 0; i < r.build.insertions.size(); ++i) {
        out << (i ? ",\n{" : "\n{");
        write_delta(r.build.insertions[i]);
        out << '}';
    }
    out << "],\n\"recoveries\": [";
    for (std::size_t i = 0; i < r.build.recoveries.size(); ++i) {
        const auto& [segment, flips] = r.build.recoveries[i];
        out << (i ? ",\n" : "\n") << "{\"segment\": [" << segment.first << ',' << segment.second << "], ";
        write_delta(flips);
        out << '}';
    }
    out << "],\n\"exterior\": {";
    write_delta(r.build.exterior);
    out << "},\n\"legalize\": {";
    write_delta(r.build.legalize);
    out << "}},\n";

    out << "\"triangles\": ";
    write_triangles(out, r.initial_triangles);
    out << ",\n\"steps\": [";
    for (std::size_t i = 0; i < r.steps.size(); ++i) {
        const Step& step = r.steps[i];
        const Cause& cause = step.cause;
        out << (i ? ",\n" : "\n") << "{\"kind\": \"" << (step.kind == Step::Kind::insert ? "insert" : "split")
            << "\", \"index\": " << step.point << ", \"point\": ";
        write_point(out, r.final_points[static_cast<std::size_t>(step.point)]);
        out << ", \"reason\": \"" << reason_name(cause.reason) << '"';
        if (cause.triangle) {
            const Tri& t = *cause.triangle;
            out << ", \"triangle\": [" << t[0] << ',' << t[1] << ',' << t[2] << ']';
        }
        if (cause.segment) out << ", \"segment\": [" << cause.segment->first << ',' << cause.segment->second << ']';
        if (cause.encroacher) out << ", \"encroacher\": " << *cause.encroacher;
        if (cause.centre) {
            out << ", \"centre\": ";
            write_point(out, *cause.centre);
        }
        if (cause.shell) out << ", \"shell\": true";
        out << ", ";
        write_delta(step);
        out << '}';
    }
    out << "\n]}\n";
}

}  // namespace mesher
