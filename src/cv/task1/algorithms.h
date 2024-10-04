#ifndef CGCV_ALGORITHMS_H
#define CGCV_ALGORITHMS_H

#include <numeric>
#include <opencv2/opencv.hpp>

typedef std::array<double, 256> Histogram;

class algorithms
{
   public:
    static void compute_grayscale(const cv::Mat &input_image, cv::Mat &grayscale_image);

    static void compute_histogram(const cv::Mat &grayscale_image, Histogram &normalized_histogram);

    static void compute_threshold(const Histogram &normalized_histogram, int &threshold);

    static void compute_binary(const cv::Mat &grayscale_image, const int threshold, cv::Mat &binary_image);

    static void draw_contours(const cv::Mat &binary_image, cv::Mat &filled_image,
                              std::vector<std::vector<cv::Point>> &contours);

    static void apply_morph_operation(const cv::Mat &filled_image, const int kernel_size, const cv::MorphTypes mode,
                                      cv::Mat &morphed_image);

    static void calc_area(const cv::Mat &normalized_image, const float conversion_factor_m2, int &area_px,
                          float &area_km2);

    static void calc_perimeter(const std::vector<std::vector<cv::Point>> &contours, const float &conversion_factor_m,
                               int &perimeter_px, float &perimeter_km);

    static void calc_center_of_mass(const cv::Mat &normalized_image, const int area_px, cv::Point &center_of_mass);

    static void crop_center_of_mass(const cv::Mat &input_image, const cv::Mat &grayscale_image,
                                    const cv::Point &center_of_mass, cv::Mat &cropped_center);

    // Bonus

    static void bonus_find_contours(cv::Mat &binary_image, std::vector<std::vector<cv::Point>> &contours);
};

#endif  // CGCV_ALGORITHMS_H
