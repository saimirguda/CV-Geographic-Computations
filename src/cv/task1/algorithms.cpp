#include "algorithms.h"

//===============================================================================
// compute_grayscale()
//-------------------------------------------------------------------------------
// TODO:  - iterate over all pixels in the image.
//        - multiply [R,G,B] with [0.2989, 0.5870, 0.1140] and accumulate.
// hints: - use cv::Vec3b type to access the color values of a 3-channel image.
//        - be aware that OpenCV treats matrix accesses in row-major order;
//          iterate through rows, then columns!
//
// parameters:
//  - input_image: [CV_8UC3] input image for the grayscale calculation
//  - grayscale_image: [CV_8UC1] output grayscale image
// return: void
//===============================================================================
void algorithms::compute_grayscale(const cv::Mat &input_image, cv::Mat &grayscale_image)
{
    grayscale_image.create(input_image.rows, input_image.cols, CV_8UC1);

    cv::Vec3f transform_weights(0.1140f, 0.5870f, 0.2989f);

    for (int i = 0; i < input_image.rows; i++) {
        for (int j = 0; j < input_image.cols; j++) {
            const auto& rgb = input_image.at<cv::Vec3b>(i, j);

            cv::Vec3f rgb_float(rgb[0], rgb[1], rgb[2]);

            float grayscale_value = rgb_float.dot(transform_weights);

            grayscale_image.at<uint8_t>(i, j) = static_cast<uint8_t>(grayscale_value); // Ayy it worked
        }
    }
}

//===============================================================================
// compute_histogram()
//-------------------------------------------------------------------------------
// TODO:  - calculate the histogram of pixel intensities in the grayscale image.
//        - normalize the histogram values.
// hints: - Histogram is a 2D array which is initialized to zero, where each
//          index, i.e. bin, corresponds to a different pixel intensity.
//
// parameters:
//  - grayscale_image: [CV_8UC1] input grayscale image
//  - normalized_histogram: [Histogram] output normalized histogram data
// return: void
//===============================================================================
void algorithms::compute_histogram(const cv::Mat &grayscale_image, Histogram &normalized_histogram)
{
    Histogram histogram = {}; // Initialize with 0

    for (int i = 0; i < grayscale_image.rows; i++) {
        for (int j = 0; j < grayscale_image.cols; j++) {
            uint8_t val = grayscale_image.at<uint8_t>(i, j);
            histogram[val]++;
        }
    }
    double total = grayscale_image.rows * grayscale_image.cols;

    for (int i = 0; i < histogram.size(); i++) {
        normalized_histogram[i] = histogram[i] / total; // normalized
    }
}

//===============================================================================
// compute_threshold()
//-------------------------------------------------------------------------------
// TODO:  - implement OTSU's method (according to the assignment) sheet to
//          calculate an optimal threshold value for binary image segmentation.
//        - use the normalized histogram to determine the inter-class variance
//          for each possible threshold value t ∈ [1,255].
//        - find the optimal threshold value by maximizing the inter-class
//          variance and save it in 'threshold'.
// hints: - the inter-class variance is calculated using the probablities and the
//          mean intensities of each class (background, foreground).
//
// parameters:
//  - normalized_histogram: [Histogram] input normalized histogram data
//  - threshold: [int] output optimal threshold value
// return: void
//===============================================================================
void algorithms::compute_threshold(const Histogram &normalized_histogram, int &threshold)
{
    const int L = 256;

    std::vector<double> sums(L, 0);
    std::vector<double> means(L, 0);

    sums[0] = normalized_histogram[0];
    for (int i = 1; i < L; ++i) {
        sums[i] = sums[i - 1] + normalized_histogram[i];
        means[i] = means[i - 1] + i * normalized_histogram[i];
    }

    threshold = 0;
    double totalMean = means[L - 1];
    double maxVariance = 0;

    for (int t = 1; t < L; t++) {
        double weight_b = sums[t - 1]; // background
        double weight_F = 1 - weight_b; // foreground


        if (weight_b == 0 || weight_F == 0) continue;

        double mB = means[t - 1] / weight_b;
        double mF = (totalMean - means[t - 1]) / weight_F;

        double variance = weight_b * weight_F * pow((mB - mF), 2);

        if (variance > maxVariance) {
            maxVariance = variance;
            threshold = t;
        }
    }
}

//===============================================================================
// compute_binary()
//-------------------------------------------------------------------------------
// TODO:  - apply the optimal threshold value to convert the grayscale image into
//          an inverted binary image.
//        - for pixel intensities above the threshold, the corresponding pixel in
//          the binary image should be set to black (0); otherwise white (255).
// hints: - compare each pixel value against the threshold to determine its new
//          color.
//
// parameters:
//  - grayscale_image: [CV_8UC1] input grayscale image
//  - threshold: [int] calculated optimal threshold value
//  - binary_image: [CV_8UC1] output binary image
// return: void
//===============================================================================
void algorithms::compute_binary(const cv::Mat &grayscale_image, const int threshold, cv::Mat &binary_image)
{
    binary_image = cv::Mat(grayscale_image.size(), CV_8UC1, 255); // makes it white

    for (int i = 0; i < grayscale_image.rows; i++) {
        for (int j = 0; j < grayscale_image.cols; j++) {
            uint8_t pixelValue = grayscale_image.at<uint8_t>(i, j);
            if (pixelValue > threshold) {
                binary_image.at<uint8_t>(i, j) = 0; // makes it black
            }
        }
    }}

//===============================================================================
// draw_contours()
//-------------------------------------------------------------------------------
// TODO:  - find contours in the binary image using cv::findContours(); use mode
//          cv::RETR_EXTERNAL to obtain outer contours only.
//        - draw all identified contours using cv::drawContours(); use the
//          cv::FILLED parameter to fill enclosed areas.
//        - save the filled image and contours to 'filled_image' and 'contours',
//          respectively.
// hints: - choose the fill color such that the enclosed area is white (255).
//
// parameters:
//  - binary_image: [CV_8UC1] input binary image
//  - filled_image: [CV_8UC1] output image with filled area
//  - contours: [std::vector<std::vector<cv::Point>>] output 2D vector to store
//              all detected contours
// return: void
//===============================================================================
void algorithms::draw_contours(const cv::Mat &binary_image, cv::Mat &filled_image,
                               std::vector<std::vector<cv::Point>> &contours)
{
    filled_image = cv::Mat::zeros(binary_image.size(), CV_8UC1);

    std::vector<cv::Vec4i> hierarchy; // unused
    cv::findContours(binary_image, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    cv::drawContours(filled_image, contours, -1, cv::Scalar(255), cv::FILLED);
}

//===============================================================================
// apply_morph_operation()
//-------------------------------------------------------------------------------
// TODO:  - apply a morphological operation (either erosion or dilation) to the
//          filled image using the specified kernel size and mode.
//        - implement the logic behind erosion and dilation according to the
//          assignment sheet.
//        - store the result in 'morphed_image'.
// hints: - the kernel size represents the area around each pixel in the image
//          considered during the morpholocial operation; the kernel size
//          is given as an odd number to ensure proper centering.
//
// parameters:
//  - filled_image: [CV_8UC1] noisy output image with filled area
//  - kernel_size: [int] size of the (square-shaped) kernel
//  - mode: [cv::MorphTypes] type of morphological operation; this will either
//          be MORPH_ERODE or MORPH_DILATE
//  - morphed_image: [CV_8UC1] output image after applied morphological operation;
//                   either dilated or eroded image, depending on mode
// return: void
//===============================================================================
void algorithms::apply_morph_operation(const cv::Mat &filled_image, const int kernel_size, const cv::MorphTypes mode,
                                       cv::Mat &morphed_image)
{
    morphed_image = cv::Mat::zeros(filled_image.size(), filled_image.type());

    int center = kernel_size / 2;

    for (int i = center; i < filled_image.rows - center; i++) {
        for (int j = center; j < filled_image.cols - center; j++) {

            bool conditionMet = (mode == cv::MORPH_ERODE);

            for (int k_i = -center; k_i <= center; k_i++) {
                for (int k_j = -center; k_j <= center; k_j++) {
                    int pixelVal = filled_image.at<uchar>(i + k_i, j + k_j);

                    if (mode == cv::MORPH_ERODE && pixelVal == 0) {
                        conditionMet = false;
                        break;
                    }

                    if (mode == cv::MORPH_DILATE && pixelVal == 255) {
                        conditionMet = true;
                        break;
                    }
                }
            }

            if (conditionMet) {
                morphed_image.at<uchar>(i, j) = 255;
            } else {
                morphed_image.at<uchar>(i, j) =  0;
            }
        }
    }

}

//===============================================================================
// calc_area()
//-------------------------------------------------------------------------------
// TODO:  - calculate the area of the filled region as the sum of pixels in the
//          normalized image.
//        - convert the resulting area from pixel units to real-world units using
//          the provided conversion factor.
//        - save the calculated areas to 'area_px' and 'area_km2'.
// hints: - the final area must be stored in square kilometers; convert your
//          result accordingly.
//
// parameters:
//  - normalized_image: [CV_8UC1] input normalized image with filled area
//  - conversion_factor_m2: [float] conversion factor from pixels to square meters;
//                          this considers the scale of the image
//  - area_px: [float] output area in pixels
//  - area_km2: [float] output area in square kilometers
// return: void
//===============================================================================
void algorithms::calc_area(const cv::Mat &normalized_image, const float conversion_factor_m2, int &area_px,
                           float &area_km2)
{
    area_px = 0;
    for (int i = 0; i < normalized_image.rows; i++) {
        for (int j = 0; j < normalized_image.cols; j++) {
            if (normalized_image.at<uint8_t>(i, j) == 1) {
                area_px += 1;
            }
        }
    }
    float area_m2 = area_px * conversion_factor_m2;

    area_km2 = area_m2 / 1000000;
}

//===============================================================================
// calc_perimeter()
//-------------------------------------------------------------------------------
// TODO:  - calculate the perimeter of the largest detected contour in the input
//          vector.
//        - For each pixel in the contour, count the outer edges by checking its
//          4-neighborhood; the total perimeter is the sum of outer edges.
//        - convert the resulting perimeter from pixel units to real-world units
//          using the provided conversion factor.
//        - save the calculated perimeters to 'perimeter_px' and 'perimeter_km'.
// hints: - find the largest contour in the given 2D vector of contours.
//        - utilize cv::pointPolygonTest() to check if neighboring pixels, i.e.,
//          edges of the pixel, lie outside the contour.
//
// parameters:
//  - contours: [std::vector<std::vector<cv::Point>>] input 2D vector with all
//              detected contours
//  - conversion_factor_m: [float] conversion factor from pixels to meters;
//                         this considers the scale of the image
//  - perimeter_px: [float] output perimeter in pixels
//  - perimeter_km: [float] output perimeter in kilometers
// return: void
//===============================================================================
void algorithms::calc_perimeter(const std::vector<std::vector<cv::Point>> &contours, const float &conversion_factor_m,
                                int &perimeter_px, float &perimeter_km)
{
    auto largestContourIt = std::max_element(contours.begin(), contours.end(),[](const std::vector<cv::Point>& a, const std::vector<cv::Point>& b) {
                                                 return cv::contourArea(a) < cv::contourArea(b);
                                             });

    if (largestContourIt == contours.end()) {
        perimeter_px = 0;
        perimeter_km = 0;
        return;
    }

    const std::vector<cv::Point> &largestContour = *largestContourIt;
    perimeter_px = 0;

    for (const auto &point : largestContour) {
        std::vector<cv::Point> neighbors = {
                cv::Point(point.x - 1, point.y), // majtas
                cv::Point(point.x + 1, point.y), // djathtas
                cv::Point(point.x, point.y - 1), // lart
                cv::Point(point.x, point.y + 1)  // posht
        };

        for (const auto &neighbor : neighbors) {
            if (cv::pointPolygonTest(largestContour, neighbor, false) < 0) {
                perimeter_px += 1;
            }
        }
    }
    float perimeter_m = perimeter_px * conversion_factor_m;
    perimeter_km = perimeter_m / 1000;
}

//===============================================================================
// calc_center_of_mass()
//-------------------------------------------------------------------------------
// TODO:  - calculate the center of mass for the filled region using pixel
//          intensities in the normalized input image.
//        - the center of mass is computed as the weighted average of the
//          spatial coordinates.
//        - save the center of mass in 'center_of_mass'.
// hints: - utilize the provided area in pixels to compute the center of mass.
//
// parameters:
//  - normalized_image: [CV_8UC1] input normalized image with filled area
//  - area_px: [float] area in pixels
//  - center_of_mass: [cv::Point] output center of mass
// return: void
//===============================================================================
void algorithms::calc_center_of_mass(const cv::Mat &normalized_image, const int area_px, cv::Point &center_of_mass)
{
    if (area_px <= 0) {
        center_of_mass = cv::Point(-1, -1);
        return;
    }
    cv::Moments m = cv::moments(normalized_image, true);
    center_of_mass.x = static_cast<int>(m.m10 / m.m00);
    center_of_mass.y = static_cast<int>(m.m01 / m.m00);
}

//===============================================================================
// crop_center_of_mass()
//-------------------------------------------------------------------------------
// TODO:  - convert the 1-channel *grayscale* image to a 3-channel image.
//        - use a circular mask, centered at the calculated center of mass, for
//          cropping; the size is given as 'crop_radius'.
//        - copy or replace the masked circular area in the grayscale image with
//          the *colored* input image.
//        - highlight the masked circular region by drawing a red circle around
//          the center.
//        - crop the final image to the correct size, specified by 'crop_size'.
//        - store the final image in 'cropped_center'.
// hints: - utilize cv::cvtColor() to convert the grayscale image; remember that
//          OpenCV uses BGR format.
//
// parameters:
//  - input_image: [CV_8UC3] input color image
//  - grayscale_image: [CV_8UC1] input grayscale image
//  - center_of_mass: [cv::Point] center of mass
//  - cropped_center: [CV_8UC3] cropped output image; highlighted center of mass
// return: void
//===============================================================================
void algorithms::crop_center_of_mass(const cv::Mat &input_image, const cv::Mat &grayscale_image,
                                     const cv::Point &center_of_mass, cv::Mat &cropped_center)
{
    const int crop_radius = 50;
    const cv::Size crop_size(200, 200);
    const cv::Scalar highlight_color(0, 0, 255);

    cv::Mat bgr_image;
    cv::cvtColor(grayscale_image, bgr_image, cv::COLOR_GRAY2BGR);

    cv::Mat mask = cv::Mat::zeros(input_image.size(), CV_8UC1);
    cv::circle(mask, center_of_mass, crop_radius, 255, cv::FILLED);
    input_image.copyTo(bgr_image, mask);

    cv::circle(bgr_image, center_of_mass, crop_radius, highlight_color, 1);
    cv::Point crop_origin(std::max(center_of_mass.x - crop_size.width / 2, 0), std::max(center_of_mass.y - crop_size.height / 2, 0));
    crop_origin.x = std::min(crop_origin.x, bgr_image.cols - crop_size.width);
    crop_origin.y = std::min(crop_origin.y, bgr_image.rows - crop_size.height);

    cv::Rect crop_rect(crop_origin, crop_size);
    cropped_center = bgr_image(crop_rect).clone();


}

void traceContour(const cv::Mat& binary_image, cv::Mat& visited, int i, int j,
                  std::vector<cv::Point>& contour) {

}



//===============================================================================
// bonus_find_contours()
//-------------------------------------------------------------------------------
// TODO:  - optional BONUS task
//        - implement the simplified contour tracing algorithm by Suzuki and Abe
//          to find contours in the binary image.
//        - save all detected contours as vectors of 2D points in 'contours'.
// hints: - follow the assignment sheet closely.
//
// parameters:
//  - binary_image: [CV_32SC1] input normalized binary image
//  - contours: [std::vector<std::vector<cv::Point>>] output 2D vector to store
//              all detected contours
// return: void
//===============================================================================
void algorithms::bonus_find_contours(cv::Mat &binary_image, std::vector<std::vector<cv::Point>> &contours)
{

}
