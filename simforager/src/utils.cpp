#include "utils.hpp"

using namespace upcxx_utils;

using std::min;
using std::string;
using std::string_view;
using std::to_string;

#include <opencv2/opencv.hpp>
#include <string>
#include <tuple>
#include <vector>
#include <filesystem>
#include "utils.hpp"  // Make sure this includes the function declarations

static cv::VideoWriter video_writer;

cv::Mat render_frame(
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &fish_points,
    int scale = 1)
{
    int scaled_width = width / scale;
    int scaled_height = height / scale;
    cv::Mat frame(scaled_height, scaled_width, CV_8UC3, cv::Scalar(0, 0, 0));
    
    // Determine dynamic fish radius based on frame size
    double radius_percent = 0.005;  // 0.5%
    int base_radius = std::max(1, static_cast<int>(radius_percent * std::min(scaled_width, scaled_height)));
    
    // Lambda for drawing points with specified radius
    auto draw_points = [&](const std::vector<std::tuple<int, int, cv::Scalar>> &points, int radius = 0) {
			 for (const auto &[x_raw, y_raw, color] : points) {
			   int x = x_raw / scale;
			   int y = y_raw / scale;
			   if (x >= 0 && x < scaled_width && y >= 0 && y < scaled_height) {
			     cv::circle(frame, cv::Point(x, y), radius, color, -1);  // filled circle
			   }
			 }
		       };
    
    draw_points(coral_points, 0);     // substrate: 1 pixel point
    draw_points(algae_points, 0);
    draw_points(sand_points, 0);
    draw_points(fish_points, base_radius);  // fish: scaled radius
    return frame;
}

/*
// Safely render a frame with clamped 2D points
cv::Mat render_frame(
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &fish_points,
    int scale  // allows rendering a lower-resolution image
) {
    int scaled_width = width / scale;
    int scaled_height = height / scale;
    cv::Mat frame(scaled_height, scaled_width, CV_8UC3, cv::Scalar(0, 0, 0));

    auto draw_points = [&](const std::vector<std::tuple<int, int, cv::Scalar>> &points, int radius) {
        for (const auto &[x_raw, y_raw, color] : points) {
            int x = x_raw / scale;
            int y = y_raw / scale;
            if (x >= 0 && x < scaled_width && y >= 0 && y < scaled_height) {
                if (radius > 0) {
                    cv::circle(frame, cv::Point(x, y), radius, color, cv::FILLED);  // BGR colour
                } else {
                    frame.at<cv::Vec3b>(y, x) = cv::Vec3b(color[2], color[1], color[0]);  // fallback
                }
            }
        }
    };

    draw_points(coral_points, 0); // Single pixel
    draw_points(algae_points, 0);
    draw_points(sand_points, 0);
    draw_points(fish_points, 2);  // Fish drawn larger with radius 2

    return frame;
}

*/

// Writes a full frame to the video composed of coral, algae, sand, and fish
void write_full_frame_to_video(
    const std::string &video_path,
    int width,
    int height,
    const std::vector<std::tuple<int, int, cv::Scalar>> &coral_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &algae_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &sand_points,
    const std::vector<std::tuple<int, int, cv::Scalar>> &fish_points)
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

    cv::Mat frame = render_frame(width, height, coral_points, algae_points, sand_points, fish_points, 1);

    SLOG("cv::Mat frame size: ", frame.cols, "x", frame.rows);
    SLOG("cv::VideoWriter expects: ", width, "x", height);
    SLOG("Frame is continuous: ", frame.isContinuous());

    if (frame.cols == 1878 && frame.rows == 1304) {
      cv::imwrite("debug_frame_1878x1304.png", frame);
      SLOG("🖼️ Wrote debug PNG: debug_frame_1878x1304.png\n");
    }
    
    video_writer.write(frame);
    SLOG("🖼️  Wrote frame to video ", video_path, " (", coral_points.size() + algae_points.size() + sand_points.size() + fish_points.size(), " points)\n");
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
    ifstream file(file_name, ios::binary);
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
            if (red == 255 && green == 0 && blue == 0) {
                color_map[i][j] = 4;
            } else if (green == 255 && red == 0 && blue == 0) {
                color_map[i][j] = 3;
            } else if (blue == 255 && red == 0 && green == 0) {
                color_map[i][j] = 2;
            } else if (red == 0 && green == 0 && blue == 0) {
                color_map[i][j] = 1;
            } else {
                DIE("Unexpected pixel value at (", i, ", ", j, "): ",
                    "R=", red, " G=", green, " B=", blue, "\n");
            }
        }
    }

    return color_map;
}

// Function to validate that a BMP was correctly read into a 2D vector color map.
// Confirms that the encoded pixel values (0–3) match the original RGB pixels.
void debugColorMapData(const string& file_name, const vector<vector<uint8_t>>& color_map) {
    ifstream file(file_name, ios::binary);
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
    int count_red = 0, count_green = 0, count_blue = 0, count_black = 0;
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

            // Determine the expected encoding
            uint8_t expected_value = 0;
            if (red == 255 && green == 0 && blue == 0) {
                expected_value = 4;
                count_red++;
            } else if (green == 255 && red == 0 && blue == 0) {
                expected_value = 3;
                count_green++;
            } else if (blue == 255 && red == 0 && green == 0) {
                expected_value = 2;
                count_blue++;
            } else if (red == 0 && green == 0 && blue == 0) {
                expected_value = 1;
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
                    DIE("Mismatch at (", i, ", ", j, "): ",
                        "Expected ", (int)expected_value,
                        " Got ", (int)color_map[i][j], "\n");
                }
            }
        }
    }

    // Print summary statistics
    SLOG("BMP Color Map Statistics:\n",
         "  Dimensions: ", width, " x ", height, "\n",
         "  Red Pixels (1):   ", count_red, "\n",
         "  Green Pixels (2): ", count_green, "\n",
         "  Blue Pixels (3):  ", count_blue, "\n",
         "  Black Pixels (0): ", count_black, "\n",
         "  Total Mismatches: ", mismatches, "\n",
         "  Total raw pixel value sum:     ", total_original_pixel_value, "\n",
         "  Total mapped value sum:        ", total_mapped_value_sum, "\n");

    if (mismatches == 0) {
        SLOG("✅ Color map matches original BMP pixel data.\n");
    } else {
        DIE("❌ Found ", mismatches, " mismatches in pixel data.\n");
    }

    // Verify that value sum matches what we expect from encoding
    if (total_mapped_value_sum == count_red * 4 + count_green * 3 + count_blue * 2 + count_black * 1) {
        SLOG("✅ Encoded value sum matches expected pixel encoding.\n");
    } else {
        DIE("❌ Mismatch in total encoded value sum.\n");
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

    ofstream file(file_name, ios::binary);
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
                default:
		  throw runtime_error("Invalid color code at " + to_string(static_cast<unsigned int>(code)) + " (" + to_string(i) + ", " + to_string(j) + ")");
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
