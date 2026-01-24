#include "utils.hpp"

using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

#include "options.hpp"
extern std::shared_ptr<Options> _options;

#include <opencv2/opencv.hpp>
#include <string>
#include <tuple>
#include <vector>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <memory>
#include "options.hpp"
extern std::shared_ptr<Options> _options;

#include <random>
#include <unordered_set>
#include <opencv2/opencv.hpp>

#include "reef.hpp"
#include "options.hpp"
#include "utils.hpp"

// Small helper: map substrate to BGR (matching your existing scheme)
static inline cv::Scalar substrate_bgr(SubstrateType t) {
  switch (t) {
    case SubstrateType::CORAL_WITH_ALGAE:  return cv::Scalar(0, 180, 165);   // olive-ish (B,G,R)
    case SubstrateType::CORAL_NO_ALGAE:    return cv::Scalar(0, 0, 255);     // red
    case SubstrateType::SAND_WITH_ALGAE:   return cv::Scalar(0, 255, 0);     // green
    case SubstrateType::SAND_NO_ALGAE:     return cv::Scalar(0, 255, 255);   // yellow
    case SubstrateType::NONE:
    default:                               return cv::Scalar(0, 0, 0);       // black
  }
}

static inline const char* substrate_short(SubstrateType t) {
  switch (t) {
    case SubstrateType::CORAL_WITH_ALGAE:  return "CWA";
    case SubstrateType::CORAL_NO_ALGAE:    return "CNA";
    case SubstrateType::SAND_WITH_ALGAE:   return "SWA";
    case SubstrateType::SAND_NO_ALGAE:     return "SNA";
    case SubstrateType::NONE:
    default:                               return "NON";
  }
}

// Show distances
namespace utils {
static inline const char* substrate_name(SubstrateType t) {
  switch (t) {
    case SubstrateType::CORAL_WITH_ALGAE:  return "CWA";
    case SubstrateType::CORAL_NO_ALGAE:    return "CNA";
    case SubstrateType::SAND_WITH_ALGAE:   return "SWA";
    case SubstrateType::SAND_NO_ALGAE:     return "SNA";
    default:                               return "UNK";
  }
}

// Forward decls assumed to exist in utils.cpp already:
//   static inline cv::Scalar substrate_bgr(SubstrateType t);

static inline const std::vector<uint16_t>& pick_distance_field(const Reef& reef, SubstrateType t) {
  // Adjust these member names to match your Reef class
  switch (t) {
    case SubstrateType::CORAL_WITH_ALGAE: return reef.D_coral_w_algae;
    case SubstrateType::CORAL_NO_ALGAE:   return reef.D_coral_no_algae;
    case SubstrateType::SAND_WITH_ALGAE:  return reef.D_sand_w_algae;
    case SubstrateType::SAND_NO_ALGAE:    return reef.D_sand_no_algae;
    default:                              return reef.D_coral_w_algae;
  }
}

std::pair<int64_t, int64_t>
nearest_grid_point(double x, double y, const std::shared_ptr<GridCoords>& grid_size)
{
    // Round to nearest integer cell centre
    int64_t gx = static_cast<int64_t>(std::floor(x + 0.5));
    int64_t gy = static_cast<int64_t>(std::floor(y + 0.5));

    // Wrap indices into valid [0, size) range using modular arithmetic
    gx = (gx % grid_size->x + grid_size->x) % grid_size->x;
    gy = (gy % grid_size->y + grid_size->y) % grid_size->y;

    return {gx, gy};
}

  
void showSubstrateDistanceContours(
    const Reef& reef,
    const Options& opt,
    SubstrateType target,
    const std::string& out_png,
    float contour_step,
    uint32_t /*rng_seed*/)
{
  const int W = opt.ecosystem_dims[0];
  const int H = opt.ecosystem_dims[1];
  if (W <= 0 || H <= 0) {
    SWARN("showSubstrateDistanceContours(): invalid ecosystem dims W=", W, " H=", H, "\n");
    return;
  }
  if (contour_step <= 0.0f) contour_step = 10.0f;

  const auto& cells = reef.get_ecosystem_cells();
  if ((int64_t)cells.size() < (int64_t)W * (int64_t)H) {
    SWARN("showSubstrateDistanceContours(): ecosystem_cells too small: ",
          cells.size(), " < ", (int64_t)W * (int64_t)H, "\n");
  }

  // ------------------------------------------------------------
  // 1) Base substrate image (preserve substrate colours)
  // ------------------------------------------------------------
  cv::Mat img(H, W, CV_8UC3, cv::Scalar(0, 0, 0));
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const int64_t id = GridCoords::to_1d(x, y, 0);
      if (id < 0 || id >= (int64_t)cells.size()) continue;
      const cv::Scalar col = substrate_bgr(cells[id]);
      img.at<cv::Vec3b>(y, x) = cv::Vec3b((uchar)col[0], (uchar)col[1], (uchar)col[2]);
    }
  }

  // ------------------------------------------------------------
  // 2) Pick the correct distance field (these are uint16_t in your build)
  // ------------------------------------------------------------
  const std::vector<uint16_t>* D = nullptr;
  switch (target) {
    case SubstrateType::CORAL_WITH_ALGAE: D = &reef.D_coral_w_algae; break;
    case SubstrateType::CORAL_NO_ALGAE:   D = &reef.D_coral_no_algae; break;
    case SubstrateType::SAND_WITH_ALGAE:  D = &reef.D_sand_w_algae; break;
    case SubstrateType::SAND_NO_ALGAE:    D = &reef.D_sand_no_algae; break;
    default: break;
  }
  if (!D) {
    SWARN("showSubstrateDistanceContours(): unknown target substrate\n");
    return;
  }

  // Convert to float image for thresholding/contours
  cv::Mat dist(H, W, CV_32F, cv::Scalar(0));
  float max_d = 0.0f;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const int64_t id = GridCoords::to_1d(x, y, 0);
      if (id < 0 || id >= (int64_t)D->size()) continue;
      const float v = (float)(*D)[id];
      dist.at<float>(y, x) = v;
      if (v > max_d) max_d = v;
    }
  }
  if (max_d <= 0.0f) {
    SWARN("showSubstrateDistanceContours(): max distance is 0; nothing to contour\n");
    return;
  }

  // ------------------------------------------------------------
  // 3) Draw contours at multiple levels
  //    We do this by thresholding dist <= level and extracting the boundary.
  // ------------------------------------------------------------
  const cv::Scalar contour_col = substrate_bgr(target);

  auto put_label = [&](cv::Mat& canvas, const std::string& text, cv::Point p, double font_scale) {
    int baseline = 0;
    const int font_face = cv::FONT_HERSHEY_SIMPLEX;
    cv::Size ts = cv::getTextSize(text, font_face, font_scale, 1, &baseline);

    // Keep text on canvas
    p.x = std::max(1, std::min(p.x, W - ts.width - 2));
    p.y = std::max(ts.height + 2, std::min(p.y, H - 2));

    // Black halo then coloured text
    cv::putText(canvas, text, p, font_face, font_scale, cv::Scalar(0, 0, 0), 3, cv::LINE_AA);
    cv::putText(canvas, text, p, font_face, font_scale, contour_col, 1, cv::LINE_AA);
  };

  // Controls for readability
  const int contour_thickness = 1;
  const int label_min_area_px = 600;          // skip tiny contours
  const int labels_per_contour = 3;           // >1 label per contour, as requested
  const int min_point_separation_px = 40;     // spacing between labels on same contour
  const float level_start = contour_step;     // skip 0 (too busy)
  const float level_end = max_d;

  // Font scale tuned for large images; tweak if needed
  const double font_scale = std::max(0.35, std::min(0.9, 700.0 / (double)std::max(W, H)));

  for (float level = level_start; level <= level_end; level += contour_step) {
    // mask = dist <= level
    cv::Mat mask;
    cv::compare(dist, level, mask, cv::CMP_LE); // mask is 0/255

    // Find boundaries of all regions for this level
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(mask, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

    if (contours.empty()) continue;

    // Draw contours
    cv::drawContours(img, contours, -1, contour_col, contour_thickness, cv::LINE_AA);

    // Add labels at multiple places on each contour
    const std::string label = cv::format("%.0f", level);

    for (const auto& c : contours) {
      if ((int)c.size() < 40) continue;

      const double area = std::fabs(cv::contourArea(c));
      if (area < (double)label_min_area_px) continue;

      // Estimate arc length to pick roughly even spacing
      const double arc = cv::arcLength(c, true);
      if (arc < 200.0) continue;

      // Decide how many labels to place (at least 2, typically 3)
      int nlab = labels_per_contour;
      if (arc < 400.0) nlab = 2;

      // Choose points along the contour: indices spaced by arc fraction.
      // We also avoid putting labels too close together in pixel space.
      std::vector<cv::Point> placed;
      placed.reserve((size_t)nlab);

      for (int i = 0; i < nlab; ++i) {
        const double frac = (i + 1.0) / (nlab + 1.0); // avoid endpoints
        int idx = (int)std::llround(frac * (double)(c.size() - 1));
        idx = std::max(0, std::min(idx, (int)c.size() - 1));
        cv::Point p = c[(size_t)idx];

        // Simple nudge off the line so text doesn’t sit exactly on the contour
        // Use local normal approximation from neighbouring points
        const cv::Point p0 = c[(size_t)std::max(0, idx - 5)];
        const cv::Point p1 = c[(size_t)std::min((int)c.size() - 1, idx + 5)];
        cv::Point2f t((float)(p1.x - p0.x), (float)(p1.y - p0.y));
        const float len = std::sqrt(t.x * t.x + t.y * t.y);
        if (len > 1e-3f) {
          t.x /= len; t.y /= len;
          // normal = (-ty, tx)
          const cv::Point2f n(-t.y, t.x);
          p.x = (int)std::llround((double)p.x + 8.0 * (double)n.x);
          p.y = (int)std::llround((double)p.y + 8.0 * (double)n.y);
        }

        // Enforce minimum spacing between labels
        bool too_close = false;
        for (const auto& q : placed) {
          const int dx = p.x - q.x;
          const int dy = p.y - q.y;
          if (dx * dx + dy * dy < min_point_separation_px * min_point_separation_px) {
            too_close = true;
            break;
          }
        }
        if (too_close) continue;

        placed.push_back(p);
        put_label(img, label, p, font_scale);
      }
    }
  }

  // Small legend in the corner (optional but helpful)
  {
    const std::string title = std::string("Contours for ") + substrate_short(target) +
                              " (step=" + cv::format("%.0f", contour_step) + ")";
    put_label(img, title, cv::Point(10, 24), font_scale);
  }

  std::vector<int> params = { cv::IMWRITE_PNG_COMPRESSION, 9 };
  if (!cv::imwrite(out_png, img, params)) {
    SWARN("showSubstrateDistanceContours(): failed to write ", out_png, "\n");
  } else {
    SLOG("🖼️  Wrote distance contour debug image: ", out_png, "\n");
  }
}
  
void showSubstrateDistances(
    const Reef& reef,
    const Options& opt,
    double fraction_to_show,
    const std::string& out_png,
    int marker_radius_px,
    uint32_t rng_seed)
{
  if (fraction_to_show <= 0.0) return;
  if (marker_radius_px < 2) marker_radius_px = 2;

  const int W = opt.ecosystem_dims[0];
  const int H = opt.ecosystem_dims[1];
  if (W <= 0 || H <= 0) {
    SWARN("showSubstrateDistances(): invalid ecosystem dims W=", W, " H=", H, "\n");
    return;
  }

  const auto& cells = reef.get_ecosystem_cells();

  // Distance fields (names as in your Reef)
  const auto& D_cwa = reef.D_coral_w_algae;
  const auto& D_cna = reef.D_coral_no_algae;
  const auto& D_swa = reef.D_sand_w_algae;
  const auto& D_sna = reef.D_sand_no_algae;

  auto dist_for = [&](SubstrateType target, int64_t id) -> float {
    switch (target) {
      case SubstrateType::CORAL_WITH_ALGAE: return (id >= 0 && id < (int64_t)D_cwa.size()) ? D_cwa[(size_t)id] : -1.0f;
      case SubstrateType::CORAL_NO_ALGAE:   return (id >= 0 && id < (int64_t)D_cna.size()) ? D_cna[(size_t)id] : -1.0f;
      case SubstrateType::SAND_WITH_ALGAE:  return (id >= 0 && id < (int64_t)D_swa.size()) ? D_swa[(size_t)id] : -1.0f;
      case SubstrateType::SAND_NO_ALGAE:    return (id >= 0 && id < (int64_t)D_sna.size()) ? D_sna[(size_t)id] : -1.0f;
      default:                              return -1.0f;
    }
  };

  // Base image exactly HxW (one pixel per cell)
  cv::Mat img(H, W, CV_8UC3, cv::Scalar(0, 0, 0));
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const int64_t id = GridCoords::to_1d(x, y, 0);
      if (id < 0 || id >= (int64_t)cells.size()) continue;
      const cv::Scalar col = substrate_bgr(cells[(size_t)id]);
      img.at<cv::Vec3b>(y, x) = cv::Vec3b((uchar)col[0], (uchar)col[1], (uchar)col[2]);
    }
  }

  // Choose random (x,y) directly to preserve the GridCoords mapping assumption
  const int64_t N2 = int64_t(W) * int64_t(H);
  int64_t k = (int64_t)std::llround(fraction_to_show * (double)N2);
  if (k < 1) k = 1;
  if (k > N2) k = N2;

  std::mt19937 rng(rng_seed);
  std::uniform_int_distribution<int> pick_x(0, W - 1);
  std::uniform_int_distribution<int> pick_y(0, H - 1);
  std::uniform_int_distribution<int> pick_target(0, 3);

  auto target_from_idx = [&](int i) -> SubstrateType {
    switch (i) {
      case 0: return SubstrateType::CORAL_WITH_ALGAE;
      case 1: return SubstrateType::CORAL_NO_ALGAE;
      case 2: return SubstrateType::SAND_WITH_ALGAE;
      default:return SubstrateType::SAND_NO_ALGAE;
    }
  };

  // Avoid duplicates
  std::unordered_set<int64_t> chosen;
  chosen.reserve((size_t)k * 2);

  auto key2d = [&](int x, int y) -> int64_t {
    return int64_t(y) * int64_t(W) + int64_t(x);
  };

  while ((int64_t)chosen.size() < k) {
    const int x = pick_x(rng);
    const int y = pick_y(rng);
    chosen.insert(key2d(x, y));
  }

  // Text tuned to be visible inside the circle
  const int font_face = cv::FONT_HERSHEY_SIMPLEX;
  const int outline_thickness = 2;

  for (int64_t idx2d : chosen) {
    const int x = int(idx2d % W);
    const int y = int(idx2d / W);

    const int64_t id = GridCoords::to_1d(x, y, 0);
    if (id < 0 || id >= (int64_t)cells.size()) continue;

    const SubstrateType cell_t   = cells[(size_t)id];
    const SubstrateType target_t = target_from_idx(pick_target(rng));

    const float d = dist_for(target_t, id);
    if (d < 0.0f) continue;

    // White fill
    cv::circle(img, cv::Point(x, y), marker_radius_px,
               cv::Scalar(255, 255, 255), cv::FILLED, cv::LINE_AA);

    // Coloured outline (use black when target==cell so it doesn't vanish on white fill)
    cv::Scalar outline = substrate_bgr(target_t);
    if (target_t == cell_t) outline = cv::Scalar(0, 0, 0);

    cv::circle(img, cv::Point(x, y), marker_radius_px,
               outline, outline_thickness, cv::LINE_AA);

    // Text inside circle: show target tag + distance
    // If you only want the number, remove the prefix.
    std::string text = std::string(substrate_short(target_t)) + ":" + cv::format("%.1f", d);

    // Scale based on radius (works for radius ~6..50+)
    const double font_scale = std::max(0.35, 0.06 * marker_radius_px);
    int baseline = 0;
    const cv::Size ts = cv::getTextSize(text, font_face, font_scale, 1, &baseline);

    cv::Point org(x - ts.width / 2, y + ts.height / 2);

    // Keep origin in bounds
    org.x = std::max(0, std::min(org.x, W - ts.width - 1));
    org.y = std::max(ts.height + 1, std::min(org.y, H - 1));

    // Black halo + coloured (or black) text for legibility on white
    cv::Scalar text_col = outline;
    if (text_col == cv::Scalar(255, 255, 255)) text_col = cv::Scalar(0, 0, 0);

    cv::putText(img, text, org, font_face, font_scale, cv::Scalar(0, 0, 0), 2, cv::LINE_AA);
    cv::putText(img, text, org, font_face, font_scale, text_col, 1, cv::LINE_AA);
  }

  std::vector<int> params = { cv::IMWRITE_PNG_COMPRESSION, 9 };
  if (!cv::imwrite(out_png, img, params)) {
    SWARN("showSubstrateDistances(): failed to write ", out_png, "\n");
  } else {
    SLOG("🖼️  Wrote substrate distance debug image: ", out_png, "\n");
  }
}

} // namespace utils

double sigmoid(double x, double midpoint, double steepness) {
    return 1.0 / (1.0 + std::exp(-steepness * (x - midpoint)));
}

static cv::VideoWriter video_writer;

// Converts intuitive RGB(r,g,b) values to OpenCV's internal BGR order
inline cv::Scalar RGB(int r, int g, int b) {
    return cv::Scalar(b, g, r);
}

cv::Mat render_frame(
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar, float, float, float, int>> &fish_points,
    int scale)
{
    int scaled_width  = width  / scale;
    int scaled_height = height / scale;

    // --------------------------------------------------
    // Base frame (substrate)
    // --------------------------------------------------
    cv::Mat frame(scaled_height, scaled_width, CV_8UC3, cv::Scalar(0, 0, 0));

    auto draw_points_px = [&](const std::vector<std::tuple<int, int, cv::Scalar>> &points) {
        for (const auto &[x_raw, y_raw, color] : points) {
            int x = x_raw / scale;
            int y = y_raw / scale;
            if ((unsigned)x < (unsigned)scaled_width &&
                (unsigned)y < (unsigned)scaled_height) {
                frame.at<cv::Vec3b>(y, x) =
                    cv::Vec3b(
                        static_cast<uchar>(color[0]),
                        static_cast<uchar>(color[1]),
                        static_cast<uchar>(color[2]));
            }
        }
    };

    draw_points_px(coral_w_algae_points);
    draw_points_px(coral_no_algae_points);
    draw_points_px(sand_w_algae_points);
    draw_points_px(sand_no_algae_points);

    // --------------------------------------------------
    // Fish overlay (ONLY for translucent discs)
    // --------------------------------------------------
    cv::Mat fish_overlay = frame.clone();

    double radius_percent = 0.005;
    int base_radius = std::max(
        2,
        static_cast<int>(
            radius_percent *
            std::min(scaled_width, scaled_height)));

    // PASS 1: draw fish discs onto overlay
    for (const auto &[x_raw, y_raw, base_color, kappa, step_length, density, thickness]
         : fish_points)
    {
        int x = x_raw / scale;
        int y = y_raw / scale;

        if ((unsigned)x >= (unsigned)scaled_width ||
            (unsigned)y >= (unsigned)scaled_height)
            continue;

        cv::circle(
            fish_overlay,
            cv::Point(x, y),
            base_radius,
            base_color,
            cv::FILLED,
            cv::LINE_AA
        );
    }

    // Blend ONCE
    const double alpha = 0.8;
    cv::addWeighted(
        fish_overlay, alpha,
        frame,        1.0 - alpha,
        0.0,
        frame
    );

    // --------------------------------------------------
    // PASS 2: rings + text (fully opaque)
    // --------------------------------------------------
    for (const auto &[x_raw, y_raw, base_color, kappa, step_length, density, thickness]
         : fish_points)
    {
        int x = x_raw / scale;
        int y = y_raw / scale;

        if ((unsigned)x >= (unsigned)scaled_width ||
            (unsigned)y >= (unsigned)scaled_height)
            continue;

        int social_radius_px = _options->social_density_radius / scale;

        // Social interaction ring
        cv::circle(
            frame,
            cv::Point(x, y),
            social_radius_px,
            cv::Scalar(255, 0, 0),
            2,
            cv::LINE_AA
        );

        // Text colours
        cv::Scalar text_color(0, 0, 0);      // black
        cv::Scalar density_color(0, 255, 0); // green

        std::string line1 = cv::format("%.2f", kappa);
        std::string line2 = cv::format("%.2f", step_length);
        std::string line3 = cv::format("%.2f", density);

        int font_face = cv::FONT_HERSHEY_SIMPLEX;
        double font_scale = 0.35;
        int font_thickness = 1;

        int b1, b2, b3;
        cv::Size s1 = cv::getTextSize(line1, font_face, font_scale, font_thickness, &b1);
        cv::Size s2 = cv::getTextSize(line2, font_face, font_scale, font_thickness, &b2);
        cv::Size s3 = cv::getTextSize(line3, font_face, font_scale, font_thickness, &b3);

        int line_spacing = 2;

        cv::Point p1(x - s1.width / 2, y - line_spacing);
        cv::Point p2(x - s2.width / 2, y + s2.height + line_spacing);
        cv::Point p3(x - s3.width / 2,
                     y + s2.height + s3.height + 2 * line_spacing);

        cv::putText(frame, line1, p1, font_face, font_scale,
                    text_color, font_thickness, cv::LINE_8);
        cv::putText(frame, line2, p2, font_face, font_scale,
                    text_color, font_thickness, cv::LINE_8);
        cv::putText(frame, line3, p3, font_face, font_scale,
                    density_color, font_thickness, cv::LINE_8);
    }

    return frame;
}

// Writes a full frame to the video composed of coral, algae, sand, and fish

void write_full_frame_to_video(
    const std::string &video_path,
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_w_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_no_algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar, float, float, float, int>> &fish_points, // <-- 7 elements
    int scale)
{
    static bool initialized = false;
    int fps = 5;

    if (!initialized) {
        video_writer.open(
            video_path,
            cv::VideoWriter::fourcc('a', 'v', 'c', '1'),
            fps,
            cv::Size(width, height),
            true);

        if (!video_writer.isOpened()) {
            SWARN("❌ Failed to open video writer at ", video_path, "\n");
            return;
        }

        SLOG("🎞️  Video writer initialised at ", video_path, "\n");
        initialized = true;
    }

    //cv::Mat frame = render_frame(width, height, coral_w_algae_points, coral_no_algae_points, sand_w_algae_points, sand_no_algae_points, fish_points, 1);

    cv::Mat frame = render_frame(width, height,
                                 coral_w_algae_points,
                                 coral_no_algae_points,
                                 sand_w_algae_points,
                                 sand_no_algae_points,
                                 fish_points,
                                 scale); //
    
    SLOG("cv::Mat frame size: ", frame.cols, "x", frame.rows);
    SLOG("cv::VideoWriter expects: ", width, "x", height);
    SLOG("Frame is continuous: ", frame.isContinuous());

    if (frame.cols == 1878 && frame.rows == 1304) {
      cv::imwrite("debug_frame_1878x1304.png", frame);
      SLOG("🖼️ Wrote debug PNG: debug_frame_1878x1304.png\n");
    }
    
    video_writer.write(frame);
    SLOG("🖼️  Wrote frame to video ", video_path, " (", coral_w_algae_points.size() + coral_no_algae_points.size() + sand_w_algae_points.size() + sand_no_algae_points.size() + fish_points.size(), " points)\n");
}

// Closes the video file when done
void finalize_video_writer() {
    if (video_writer.isOpened()) {
        video_writer.release();
        SLOG("✅ Finalized video writer\n");
    }
}


// Function to read a 24-bit BMP file into a 2D vector where:
//  0 = Black pixel (0, 0, 0)
//  1 = Red-only pixel (255, 0, 0)
//  2 = Green-only pixel (0, 255, 0)
//  3 = Blue-only pixel (0, 0, 255)
// Assumes all pixel values are exactly 0 or 255.
vector<vector<uint8_t>> readBMPColorMap(const string& file_name) {
  ifstream file(file_name, std::ios::binary);
    if (!file) {
        DIE("Failed to open file ", file_name, "\n");
        return {};
    }

    // Read standard BMP header (54 bytes)
    uint8_t header[54];
    file.read(reinterpret_cast<char*>(header), 54);

    // Extract image metadata
    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int data_offset = *(int*)&header[10];
    int row_padded = (width * 3 + 3) & (~3);  // Each row is padded to a 4-byte boundary

    // Prepare the output 2D vector to hold the encoded pixel values
    vector<vector<uint8_t>> color_map(height, vector<uint8_t>(width, 0));

    // Seek to the beginning of pixel data
    file.seekg(data_offset);

    // BMP stores rows bottom-up, so iterate in reverse row order
    for (int i = height - 1; i >= 0; --i) {
        vector<uint8_t> row(row_padded);
        file.read(reinterpret_cast<char*>(row.data()), row_padded);

        for (int j = 0; j < width; ++j) {
            uint8_t blue  = row[j * 3 + 0];
            uint8_t green = row[j * 3 + 1];
            uint8_t red   = row[j * 3 + 2];

            // Map RGB values to encoded colour class
            if(red == 255 && green == 255 && blue == 0) {
	      // Sand no algae
	      color_map[i][j] = 0;
            }else if (red == 0 && green == 255 && blue == 0) {
	      // Sand with Algae
	      color_map[i][j] = 1;
            } else if (red == 255 && green == 0 && blue == 0) {
	      // Coral no algae
	      color_map[i][j] = 2;
            } else if (red == 165 && green == 180 && blue == 0) {
	      // Coral with algae
                color_map[i][j] = 3;
            } else if (red == 0 && green == 0 && blue == 0) {
	      // None
                color_map[i][j] = 4;
            } else {
                DIE("readBMPColorMap(): Unexpected pixel value at (", i, ", ", j, "): ",
                    "R=", red, " G=", green, " B=", blue, "\n");
            }
        }
    }

    return color_map;
}

// Function to validate that a BMP was correctly read into a 2D vector color map.
// Confirms that the encoded pixel values (0–3) match the original RGB pixels.
void debugColorMapData(const string& file_name, const vector<vector<uint8_t>>& color_map) {
  ifstream file(file_name, std::ios::binary);
    if (!file) {
        DIE("Failed to reopen file ", file_name, "\n");
        return;
    }

    // Read BMP header
    uint8_t header[54];
    file.read(reinterpret_cast<char*>(header), 54);

    // Extract image info
    int width = *(int*)&header[18];
    int height = *(int*)&header[22];
    int data_offset = *(int*)&header[10];
    int row_padded = (width * 3 + 3) & (~3);

    // Verify dimensions match
    if (color_map.size() != height || color_map[0].size() != width) {
        DIE("ERROR: Image dimensions don't match BMP header.\n");
        return;
    }

    // Seek to pixel data
    file.seekg(data_offset);

    // Stats
    int count_red = 0, count_green = 0, count_blue = 0, count_olive = 0, count_black = 0, count_yellow = 0;
    int mismatches = 0;
    uint64_t total_original_pixel_value = 0;
    uint64_t total_mapped_value_sum = 0;

    for (int i = height - 1; i >= 0; --i) {
        vector<uint8_t> row(row_padded);
        file.read(reinterpret_cast<char*>(row.data()), row_padded);

        for (int j = 0; j < width; ++j) {
            uint8_t blue  = row[j * 3 + 0];
            uint8_t green = row[j * 3 + 1];
            uint8_t red   = row[j * 3 + 2];

            // Accumulate total pixel brightness
            total_original_pixel_value += red + green + blue;

	    // Determine the expected encoding based on RGB colour
	    uint8_t expected_value = 0;
	    
	    if (red == 255 && green == 255 && blue == 0) {
	      // Sand (light yellow)
	      expected_value = 0;
	      count_yellow++;
	    } else if (red == 0 && green == 255 && blue == 0) {
	      // Algae on sand (vibrant green)
	      expected_value = 1;
	      count_green++;
	    } else if (red == 255 && green == 0 && blue == 0) {
	      // Coral with no algae (pure red)
	      expected_value = 2;
	      count_red++;
	    } else if (red == 165 && green == 180 && blue == 0) {
	      // Algae on coral (dark olive green)
	      expected_value = 3;
	      count_olive++;
	    } else if (red == 0 && green == 0 && blue == 0) {
	      // None / background
	      expected_value = 4;
	      count_black++;
	    } else {
	      DIE("Unexpected pixel value at (", i, ", ", j, "): ",
		  "R=", red, " G=", green, " B=", blue, "\n");
	    }
	    
            // Compare against the decoded color map
            total_mapped_value_sum += color_map[i][j];

            if (color_map[i][j] != expected_value) {
                mismatches++;
                if (mismatches <= 10) {
                    DIE("Decode Mismatch at (", i, ", ", j, "): ",
                        "Expected ", (int)expected_value,
                        " Got ", (int)color_map[i][j], "\n");
                }
            }
        }
    }

    // Print summary statistics
    SLOG("BMP Color Map Statistics:\n",
	 "  Dimensions: ", width, " x ", height, "\n",
	 "  Yellow Pixels (0): ", count_yellow, "\n",
	 "  Green Pixels  (1): ", count_green, "\n",
	 "  Red Pixels    (2): ", count_red, "\n",
	 "  Olive Pixels  (3): ", count_olive, "\n",
	 "  Black Pixels  (—): ", count_black, "\n",
	 "  Total Mismatches: ", mismatches, "\n",
	 "  Total raw pixel value sum: ", total_original_pixel_value, "\n",
	 "  Total mapped value sum:    ", total_mapped_value_sum, "\n");
    
    if (mismatches == 0) {
      SLOG("✅ Color map matches original BMP pixel data.\n");
    } else {
      DIE("❌ Found ", mismatches, " mismatches in pixel data.\n");
    }
    
    // Verify that total encoded value sum matches the zero-based encoding
    // 0 = yellow (sand), 1 = green (algae on sand),
    // 2 = red (coral), 3 = olive (algae on coral)
    uint64_t expected_sum =
      count_yellow * 0
      + count_green  * 1
      + count_red    * 2
      + count_olive  * 3;
    
    if (total_mapped_value_sum == expected_sum) {
      SLOG("✅ Encoded value sum matches expected pixel encoding.\n");
    } else {
      DIE("❌ Mismatch in total encoded value sum.\n",
	  "Expected: ", expected_sum, ", Got: ", total_mapped_value_sum, "\n");
    }
}


void writeBMPColorMap(const string& file_name, const vector<vector<uint8_t>>& substrate_array) {
    int height = substrate_array.size();
    if (height == 0) return;
    int width = substrate_array[0].size();

    int row_padded = (width * 3 + 3) & (~3);
    int data_size = row_padded * height;
    int file_size = 54 + data_size;

    uint8_t header[54] = {0};

    // BMP Header
    header[0] = 'B';
    header[1] = 'M';
    *(int*)&header[2] = file_size;
    *(int*)&header[10] = 54;
    *(int*)&header[14] = 40;       // DIB header size
    *(int*)&header[18] = width;
    *(int*)&header[22] = height;
    *(short*)&header[26] = 1;      // Color planes
    *(short*)&header[28] = 24;     // Bits per pixel
    *(int*)&header[34] = data_size;

    ofstream file(file_name, std::ios::binary);
    file.write(reinterpret_cast<char*>(header), 54);

    // Write pixel data bottom-up
    vector<uint8_t> row(row_padded);
    for (int i = height - 1; i >= 0; --i) {
        for (int j = 0; j < width; ++j) {
            uint8_t code = substrate_array[i][j];
            uint8_t r = 0, g = 0, b = 0;

            switch (code) {
                case 1: r = g = b = 0; break;
                case 2: b = 255; break;
                case 3: g = 255; break;
                case 4: r = 255; break;
                case 5: r=255; g=255; b=0; break;  // yellow (NEW)
                default:
		  throw std::runtime_error("Invalid color code at " + to_string(static_cast<unsigned int>(code)) + " (" + to_string(i) + ", " + to_string(j) + ")");
            }

            row[j * 3 + 0] = b;
            row[j * 3 + 1] = g;
            row[j * 3 + 2] = r;
        }

        // Pad the end of the row with zeroes if needed
        for (int p = width * 3; p < row_padded; ++p) {
            row[p] = 0;
        }

        file.write(reinterpret_cast<char*>(row.data()), row_padded);
    }

    file.close();
}

int pin_thread(pid_t pid, int cid) {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(cid, &cpu_set);
  if (sched_setaffinity(pid, sizeof(cpu_set), &cpu_set) == -1) {
    if (errno == 3) WARN("%s, pid: %d", strerror(errno), pid);
    return -1;
  }
  return 0;
}

void dump_single_file(const string &fname, const string &out_str) {
  auto sz = out_str.length();
  upcxx::atomic_domain<size_t> ad({upcxx::atomic_op::fetch_add, upcxx::atomic_op::load});
  upcxx::global_ptr<size_t> fpos = nullptr;
  // always give rank 0 the first chunk so it can write header info if needed
  if (!upcxx::rank_me()) fpos = upcxx::new_<size_t>(sz);
  fpos = upcxx::broadcast(fpos, 0).wait();
  size_t my_fpos = 0;
  if (upcxx::rank_me()) my_fpos = ad.fetch_add(fpos, sz, std::memory_order_relaxed).wait();
  // wait until all ranks have updated the global counter
  upcxx::barrier();
  int fileno = -1;
  size_t fsize = 0;
  if (!upcxx::rank_me()) {
    fsize = ad.load(fpos, std::memory_order_relaxed).wait();
    // rank 0 creates the file and truncates it to the correct length
    fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) WARN("Error trying to create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, fsize) == -1)
      WARN("Could not truncate ", fname, " to ", fsize, " bytes\n");
  }
  upcxx::barrier();
  ad.destroy();
  // wait until rank 0 has finished setting up the file
  if (rank_me()) fileno = open(fname.c_str(), O_WRONLY, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (fileno == -1) WARN("Error trying to open file ", fname, ": ", strerror(errno), "\n");
  auto bytes_written = pwrite(fileno, out_str.c_str(), sz, my_fpos);
  close(fileno);
  if (bytes_written != sz)
    DIE("Could not write all ", sz, " bytes; only wrote ", bytes_written, "\n");
  upcxx::barrier();
  auto tot_bytes_written = upcxx::reduce_one(bytes_written, upcxx::op_fast_add, 0).wait();
  SLOG_VERBOSE("Successfully wrote ", get_size_str(tot_bytes_written), " to ", fname, "\n");
}

std::shared_ptr<Random> _rnd_gen;

void write_test_video(const std::string& output_path) {
    int width = 512;
    int height = 512;
    int num_frames = 50;
    int fps = 4;

    // Create VideoWriter with dynamic output path
    cv::VideoWriter writer(output_path,
                       cv::VideoWriter::fourcc('a', 'v', 'c', '1'),
                       fps,
                       cv::Size(width, height),
                       true); // color
    
    if (!writer.isOpened()) {
        SWARN("Failed to open video writer for path: ", output_path, "\n");
        return;
    }

    for (int frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
        cv::Mat frame(height, width, CV_8UC3, cv::Scalar(0, 0, 0)); // black background

        int x = frame_idx;
        int y = frame_idx;
        cv::Rect dot_rect(x, y, 10, 10);

        cv::rectangle(frame, dot_rect, cv::Scalar(0, 0, 255), cv::FILLED); // red square

        writer.write(frame);
    }

    SLOG("Wrote test video to: ", output_path, "\n");
}

//grazer trajectory logger
namespace {
  // Keep one append stream per grazer id per process to avoid re-open costs
  static std::unordered_map<std::string, std::unique_ptr<std::ofstream>> g_streams;
  static std::filesystem::path g_tracks_dir;
  static bool g_dir_ready = false;

  inline void ensure_tracks_dir() {
    if (g_dir_ready) return;
    // Use the run's output directory if available
    std::filesystem::path base = _options ? std::filesystem::path(_options->output_dir)
                                          : std::filesystem::current_path();
    g_tracks_dir = base / "grazer_tracks";
    std::error_code ec;
    std::filesystem::create_directories(g_tracks_dir, ec);
    g_dir_ready = true;
  }
}


void log_grazer_step(const std::string& fish_id, int timestep,
                     int64_t x, int64_t y, int64_t z, int substrate, float kappa) {
  ensure_tracks_dir();
  auto it = g_streams.find(fish_id);
  if (it == g_streams.end()) {
    const auto fname = (g_tracks_dir / ("grazer_" + fish_id + ".csv")).string();

    // If file doesn't exist, we’ll write CSV header first
    const bool exists = std::filesystem::exists(fname);
    auto ofs = std::make_unique<std::ofstream>(fname, std::ios::app);
    if (!ofs->good()) {
      SWARN("Could not open trajectory file ", fname, " for fish ", fish_id, "\n");
      return;
    }
    if (!exists) (*ofs) << "timestep,x,y,z,substrate,kappa\n";
    it = g_streams.emplace(fish_id, std::move(ofs)).first;
  }

  auto& out = *(it->second);
  out << timestep << ',' << x << ',' << y << ',' << z << "," << substrate << "," << kappa << '\n';
}

void finalize_grazer_logs() {
  for (auto &kv : g_streams) {
    if (kv.second && kv.second->is_open()) kv.second->flush();
    // Let unique_ptr close on destruction
  }
  g_streams.clear();
}


// ======================================================================
// Algae consumption + grazer time distribution logger
//   Writes a single CSV in: <output_dir>/population_stats/algae_grazer_summary.csv
// ======================================================================

namespace {
  static std::ofstream g_algae_grazer_ofs;
  static bool g_algae_grazer_ready = false;

  inline void ensure_algae_grazer_file() {
    if (g_algae_grazer_ready) return;

    std::filesystem::path base =
        _options ? std::filesystem::path(_options->output_dir)
                 : std::filesystem::current_path();

    std::filesystem::path stats_dir = base / "population_stats";
    std::error_code ec;
    std::filesystem::create_directories(stats_dir, ec);

    std::filesystem::path csv_path = stats_dir / "population_summary.csv";

    const bool exists = std::filesystem::exists(csv_path);
    g_algae_grazer_ofs.open(csv_path, std::ios::app);

    if (!g_algae_grazer_ofs.good()) {
      SWARN("Could not open summary file ",
            csv_path.string(), "\n");
      return;
    }

    // Write header once
    if (!exists) {
      g_algae_grazer_ofs
          << "total_timesteps,"
          << "grazer_count,"
          << "coral_w_algae_timesteps,coral_w_algae_percent,"
          << "coral_no_algae_timesteps,coral_no_algae_percent,"
          << "sand_w_algae_timesteps,sand_w_algae_percent,"
          << "sand_no_algae_timesteps,sand_no_algae_percent,"
          << "initial_algae,final_algae,difference,reduction_percent\n";
    }

    g_algae_grazer_ready = true;
  }
} // anonymous namespace

void log_algae_and_grazer_stats(int64_t total_timesteps,
                                int64_t grazer_count,
                                int64_t coral_w_algae_steps,
                                int64_t coral_no_algae_steps,
                                int64_t sand_w_algae_steps,
                                int64_t sand_no_algae_steps,
                                double initial_algae,
                                double final_algae) {
  ensure_algae_grazer_file();
  if (!g_algae_grazer_ofs.good()) return;

  auto pct = [](int64_t a, int64_t t) -> double {
    return (t > 0) ? (100.0 * double(a) / double(t)) : 0.0;
  };

  const double diff = initial_algae - final_algae;
  const double reduction_pct =
      (initial_algae > 0.0) ? (diff / initial_algae * 100.0) : 0.0;

  double coral_w_steps_per_grazer  = double(coral_w_algae_steps)  / grazer_count;
  double coral_no_steps_per_grazer = double(coral_no_algae_steps) / grazer_count;
  double sand_w_steps_per_grazer   = double(sand_w_algae_steps)   / grazer_count;
  double sand_no_steps_per_grazer  = double(sand_no_algae_steps)  / grazer_count;

  g_algae_grazer_ofs << total_timesteps << ',' 
                     << grazer_count << ','
                     << coral_w_steps_per_grazer << ','
                     << std::fixed << std::setprecision(2)
                     << pct(coral_w_algae_steps, (total_timesteps* grazer_count)) << ','
                     << coral_no_steps_per_grazer << ','
                     << pct(coral_no_algae_steps, (total_timesteps* grazer_count)) << ','
                     << sand_w_steps_per_grazer << ','
                     << pct(sand_w_algae_steps, (total_timesteps* grazer_count)) << ','
                     << sand_no_steps_per_grazer << ','
                     << pct(sand_no_algae_steps, (total_timesteps* grazer_count)) << ','
                     << initial_algae << ','
                     << final_algae << ','
                     << diff << ','
                     << reduction_pct << '\n';

  g_algae_grazer_ofs.flush();
}



