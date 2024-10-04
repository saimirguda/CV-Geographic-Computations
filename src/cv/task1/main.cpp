#include <dirent.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "algorithms.h"
#include "opencv2/opencv.hpp"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

// #define FULL_VERSION 1
// #define FINAL_RUN 1
// #define GENERATE_REF 1

#define RST "\x1B[0m"
#define KRED "\x1B[31m"
#define KGRN "\x1B[32m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST

#define BOLD(x) "\x1B[1m" x RST

#if GENERATE_REF || FINAL_RUN
struct reference
{
    cv::Mat *mat;
    cv::Point *point;
    int *value_px;
    float *value_real;

    reference(cv::Mat *m) : mat(m), point(nullptr), value_px(nullptr), value_real(nullptr) {}

    reference(cv::Point *pt) : mat(nullptr), point(pt), value_px(nullptr), value_real(nullptr) {}

    reference(int *val) : mat(nullptr), point(nullptr), value_px(val), value_real(nullptr) {}

    reference(float *val) : mat(nullptr), point(nullptr), value_px(nullptr), value_real(val) {}
};
#endif

//===============================================================================
// Configuration
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
struct Config
{
    // morph operations
    int kernel_size_morph = 0;
    int kernel_size_perimeter = 0;

    // scaling factors
    float px_to_m = 0.f;
    float scale_adjustment_divisor_area = 0.f;
    float scale_adjustment_divisor_perimeter = 0.f;
};

//===============================================================================
// make_directory()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void make_directory(const char *path)
{
#if defined(_WIN32)
    _mkdir(path);
#else
    mkdir(path, 0777);
#endif
}

//===============================================================================
// is_path_existing()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
bool is_path_existing(const char *path)
{
    struct stat buffer
    {
    };
    return (stat(path, &buffer)) == 0;
}

#if GENERATE_REF
//===============================================================================
// generate_ref()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void generate_ref(std::string path, reference ref)
{
    cv::FileStorage fs(path, cv::FileStorage::WRITE);
    if (ref.mat != nullptr)
    {
        fs << "image" << *ref.mat;
    }
    else if (ref.point != nullptr)
    {
        fs << "image" << *ref.point;
    }
    else if (ref.value_px != nullptr)
    {
        fs << "image" << *ref.value_px;
    }
    else if (ref.value_real != nullptr)
    {
        fs << "image" << *ref.value_real;
    }
}
#endif

#if FINAL_RUN
//===============================================================================
// get_ref_image()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void get_ref_image(std::string ref_directory, std::string name, reference ref)
{
    struct dirent *entry;
    DIR *dir = opendir((ref_directory).c_str());
    while ((entry = readdir(dir)) != NULL)
    {
        std::string entry_name = entry->d_name;
        if (entry_name.find(name) != std::string::npos)
        {
            std::string full_path = ref_directory + name;
            cv::FileStorage fs(full_path, cv::FileStorage::READ);
            if (ref.mat != nullptr)
            {
                fs["image"] >> *ref.mat;
            }
            else if (ref.point != nullptr)
            {
                fs["image"] >> *ref.point;
            }
            else if (ref.value_px != nullptr)
            {
                fs["image"] >> *ref.value_px;
            }
            else if (ref.value_real != nullptr)
            {
                fs["image"] >> *ref.value_real;
            }
            break;
        }
    }
    closedir(dir);
}
#endif

//===============================================================================
// save_image()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void save_image(const std::string &out_directory, const std::string &name, const size_t number, const cv::Mat &image)
{
    std::stringstream number_stringstream;
    number_stringstream << std::setfill('0') << std::setw(2) << number;
    std::string path = out_directory + number_stringstream.str() + "_" + name + ".png";
    cv::imwrite(path, image);
    std::cout << "saving image: " << path << std::endl;
}

//===============================================================================
// save_value()
//-------------------------------------------------------------------------------
// std::ios_base::openmode
//  + std::ios_base::trunc  -> for first call, clear file before write
//  + std::ios_base::app    -> for second call, if you want to append to the file (new line)
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void save_value(const std::string &out_directory, const std::string &name, const size_t number,
                const std::string &field_name, const float value, const int precision,
                const std::ios_base::openmode mode = std::ios_base::trunc)
{
    std::stringstream number_stringstream;
    number_stringstream << std::setfill('0') << std::setw(2) << number;
    std::string path = out_directory + number_stringstream.str() + "_" + name + ".txt";
    std::ofstream file(path, mode);
    if (file.is_open())
    {
        file << field_name << ": " << std::fixed << std::setprecision(precision) << value << std::endl;
        file.close();
        std::cout << "saving value: " << path << std::endl;
    }
    else
    {
        std::cout << BOLD(FRED("[ERROR]")) << " Could not save value to (" << path << ")" << std::endl;
    }
}

//===============================================================================
// run()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void run(const cv::Mat &input_image, const std::string &out_directory, const std::string &ref_directory, Config config)
{
    size_t output_counter = 0;

    //=============================================================================
    // Grayscale image
    //=============================================================================
    std::cout << "Step  1 - calculating grayscale image... " << std::endl;
    cv::Mat grayscale = cv::Mat::zeros(input_image.size(), CV_8UC1);
    algorithms::compute_grayscale(input_image, grayscale);
    save_image(out_directory, "grayscale", ++output_counter, grayscale);
#if GENERATE_REF
    generate_ref(ref_directory + "01_grayscale.json", reference{&grayscale});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "01_grayscale.json", reference{&grayscale});
#endif

    //=============================================================================
    // Binary image
    //=============================================================================
    std::cout << "Step  2a - calculating histogram... " << std::endl;
    Histogram histogram = {0.f};
    algorithms::compute_histogram(grayscale, histogram);

    std::cout << "Step  2b - calculating optimal threshold... " << std::endl;
    // init at middle of 8-bit grayscale histogram (0 to 255)
    int threshold = 128;
    algorithms::compute_threshold(histogram, threshold);
    std::cout << "+--- Using threshold value " << threshold << std::endl;

    std::cout << "Step  2c - calculating binary image... " << std::endl;
    cv::Mat binary = cv::Mat::zeros(input_image.size(), CV_8UC1);
    algorithms::compute_binary(grayscale, threshold, binary);
    save_image(out_directory, "binary", ++output_counter, binary);
#if GENERATE_REF
    generate_ref(ref_directory + "02_binary.json", reference{&binary});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "02_binary.json", reference{&binary});
#endif

    //=============================================================================
    // Find & draw contours
    //=============================================================================
    std::cout << "Step  3 - calculate and fill contour of border... " << std::endl;
    std::vector<std::vector<cv::Point>> contours;
    cv::Mat filled = cv::Mat::zeros(input_image.size(), CV_8UC1);
    algorithms::draw_contours(binary, filled, contours);
    // save intermediate
    cv::Mat img_contours = cv::Mat::zeros(input_image.size(), CV_8UC1);
    cv::drawContours(img_contours, contours, -1, cv::Scalar(255));
    save_image(out_directory, "contours", ++output_counter, img_contours);
    save_image(out_directory, "contours_filled", output_counter, filled);
#if GENERATE_REF
    generate_ref(ref_directory + "03_contours.json", reference{&img_contours});
    generate_ref(ref_directory + "03_contours_filled.json", reference{&filled});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "03_contours.json", reference{&img_contours});
    get_ref_image(ref_directory, "03_contours_filled.json", reference{&filled});
#endif

    //=============================================================================
    // Morphological operation(s)
    //=============================================================================
    std::cout << "Step  4 - applying morphological operations (ERODE)... " << std::endl;
    cv::Mat eroded_image = cv::Mat::zeros(input_image.size(), CV_8UC1);
    algorithms::apply_morph_operation(filled, config.kernel_size_morph, cv::MORPH_ERODE, eroded_image);
    save_image(out_directory, "erosion_applied", ++output_counter, eroded_image);

    std::cout << "Step  5 - applying morphological operations (DILATE)... " << std::endl;
    cv::Mat opened_image = cv::Mat::zeros(input_image.size(), CV_8UC1);
    algorithms::apply_morph_operation(eroded_image, config.kernel_size_morph, cv::MORPH_DILATE, opened_image);
    save_image(out_directory, "opening_applied", ++output_counter, opened_image);
#if GENERATE_REF
    generate_ref(ref_directory + "04_erosion_applied.json", reference{&eroded_image});
    generate_ref(ref_directory + "05_opening_applied.json", reference{&opened_image});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "04_erosion_applied.json", reference{&eroded_image});
    get_ref_image(ref_directory, "05_opening_applied.json", reference{&opened_image});
#endif

    // normalize image
    cv::Mat normalized;
    cv::normalize(opened_image, normalized, 0, 1.0, cv::NORM_MINMAX);

    //=============================================================================
    // Calculate area
    //=============================================================================
    std::cout << "Step  6 - calculating area... " << std::endl;
    int area_px = 0;
    float area_km2 = 0.0f;
    float conversion_factor_m2 = config.px_to_m * config.px_to_m / config.scale_adjustment_divisor_area;
    algorithms::calc_area(normalized, conversion_factor_m2, area_px, area_km2);
    std::cout << "+--- Area in pixels: " << area_px << std::endl;
    std::cout << "+--- Area in km^2: " << area_km2 << std::endl;
    save_value(out_directory, "area", ++output_counter, "Area (pixels)", area_px, 0);
    save_value(out_directory, "area", output_counter, "Area (km^2)", area_km2, 3, std::ios_base::app);
#if GENERATE_REF
    generate_ref(ref_directory + "06_area_px.json", reference{&area_px});
    generate_ref(ref_directory + "06_area_km2.json", reference{&area_km2});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "06_area_px.json", reference{&area_px});
    get_ref_image(ref_directory, "06_area_km2.json", reference{&area_km2});
#endif

    //=============================================================================
    // Calculate perimeter
    //=============================================================================
    std::cout << "Step  7 - calculating perimeter... " << std::endl;
    cv::Mat perimeter_image = cv::Mat::zeros(input_image.size(), CV_8UC1);
    algorithms::apply_morph_operation(binary, config.kernel_size_perimeter, cv::MORPH_ERODE, perimeter_image);
    std::vector<std::vector<cv::Point>> contours_perimeter;
    cv::findContours(perimeter_image, contours_perimeter, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_NONE);
    int perimeter_px = 0;
    float perimeter_km = 0.0f;
    float conversion_factor_m = config.px_to_m / config.scale_adjustment_divisor_perimeter;
    algorithms::calc_perimeter(contours_perimeter, conversion_factor_m, perimeter_px, perimeter_km);
    std::cout << "+--- Perimeter in pixels: " << perimeter_px << std::endl;
    std::cout << "+--- Perimeter in km: " << perimeter_km << std::endl;
    save_value(out_directory, "perimeter", ++output_counter, "Perimeter (pixels)", perimeter_px, 0);
    save_value(out_directory, "perimeter", output_counter, "Perimeter (km)", perimeter_km, 3, std::ios_base::app);
#if GENERATE_REF
    generate_ref(ref_directory + "07_perimeter_px.json", reference{&perimeter_px});
    generate_ref(ref_directory + "07_perimeter_km.json", reference{&perimeter_km});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "07_perimeter_px.json", reference{&perimeter_px});
    get_ref_image(ref_directory, "07_perimeter_km.json", reference{&perimeter_km});
#endif

    //=============================================================================
    // Calculate center of mass
    //=============================================================================
    std::cout << "Step  8 - calculating center of mass... " << std::endl;
    // init at image center
    cv::Point center_of_mass = {normalized.cols / 2, normalized.rows / 2};
    algorithms::calc_center_of_mass(normalized, area_px, center_of_mass);
    std::cout << "+--- Center of mass: " << center_of_mass << std::endl;
    save_value(out_directory, "center_of_mass", ++output_counter, "Center x", center_of_mass.x, 0);
    save_value(out_directory, "center_of_mass", output_counter, "Center y", center_of_mass.y, 0, std::ios_base::app);
#if GENERATE_REF
    generate_ref(ref_directory + "08_center_of_mass.json", reference{&center_of_mass});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "08_center_of_mass.json", reference{&center_of_mass});
#endif

    //=============================================================================
    // Crop center of mass
    //=============================================================================
    std::cout << "Step 9 - cropping center of mass... " << std::endl;
    cv::Mat cropped_center = cv::Mat::zeros(input_image.size(), CV_8UC3);
    algorithms::crop_center_of_mass(input_image, grayscale, center_of_mass, cropped_center);
    save_image(out_directory, "cropped_center", ++output_counter, cropped_center);
#if GENERATE_REF
    generate_ref(ref_directory + "09_cropped_center.json", reference{&cropped_center});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "09_cropped_center.json", reference{&cropped_center});
#endif

    // ==============================================================================
    // Bonus: Find contours
    // ==============================================================================
    std::cout << "Bonus  - computing contours... " << std::endl;
    std::vector<std::vector<cv::Point>> contours_bonus;
    cv::Mat binary_bonus = cv::Mat::zeros(input_image.size(), CV_32SC1);
    cv::normalize(binary, binary_bonus, 0, 1, cv::NORM_MINMAX, CV_32SC1);
    algorithms::bonus_find_contours(binary_bonus, contours_bonus);
    // get the largest contour and draw it
    cv::Mat img_contours_bonus = cv::Mat::zeros(input_image.size(), CV_8UC1);
    cv::drawContours(img_contours_bonus, contours_bonus, -1, cv::Scalar(255));
    save_image(out_directory + "bonus/", "contours_bonus", ++output_counter, img_contours_bonus);
#if GENERATE_REF
    generate_ref(ref_directory + "10_bonus_contours.json", reference{&img_contours_bonus});
#endif
#if FINAL_RUN
    get_ref_image(ref_directory, "10_bonus_contours.json", reference{&img_contours_bonus});
#endif
}

//===============================================================================
// execute_testcase()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
void execute_testcase(const rapidjson::Value &config_data)
{
    //=============================================================================
    // Parse input data
    //=============================================================================
    std::string name = config_data["name"].GetString();
    std::string image_path = config_data["image_path"].GetString();

    Config config;

    // morph operations
    config.kernel_size_morph = (int)config_data["kernel_size_morph"].GetUint();
    config.kernel_size_perimeter = (int)config_data["kernel_size_perimeter"].GetUint();

    // scaling factors
    config.px_to_m = (float)config_data["px_to_m"].GetDouble();
    config.scale_adjustment_divisor_area = (float)config_data["scale_adjustment_divisor_area"].GetDouble();
    config.scale_adjustment_divisor_perimeter = (float)config_data["scale_adjustment_divisor_perimeter"].GetDouble();

    //=============================================================================
    // Load input images
    //=============================================================================
    std::cout << BOLD(FGRN("[INFO]")) << " Input image: " << image_path << std::endl;

    cv::Mat img = cv::imread(image_path);

    if (!img.data)
    {
        std::cout << BOLD(FRED("[ERROR]")) << " Could not load image (" << image_path << ")" << std::endl;
        throw std::runtime_error("Could not load file");
    }

    //=============================================================================
    // Create output directory
    //=============================================================================
    std::string output_directory = "output/" + name + "/";

    std::cout << BOLD(FGRN("[INFO]")) << " Output path: " << output_directory << std::endl;

    make_directory("output/");
    make_directory(output_directory.c_str());
    // create bonus directory
    make_directory((output_directory + "/bonus/").c_str());

    std::string ref_path = "data/intm/";
    std::string ref_directory = ref_path + name + "/";

#if FINAL_RUN
    if (!is_path_existing(ref_directory.c_str()))
    {
        std::cout << BOLD(FRED("[ERROR]")) << " ref directory does not exist!" << std::endl;
        std::cout << BOLD(FGRN("[INFO]")) << " execute with GENERATE_REF 1 first" << std::endl;
        throw std::runtime_error("Could not load ref files");
    }
    else
    {
        std::cout << "opening ref directory" << std::endl;
    }
#endif

#if GENERATE_REF
    make_directory(ref_path.c_str());
    make_directory(ref_directory.c_str());
#endif

    //=============================================================================
    // Starting default task
    //=============================================================================
    std::cout << "Starting MAIN Task..." << std::endl;
    run(img, output_directory, ref_directory, config);
}

//===============================================================================
// main()
//-------------------------------------------------------------------------------
// TODO:
//  - Nothing!
//  - Do not change anything here
//===============================================================================
int main(int argc, char *argv[])
{
    std::cout << "CV/task1 framework version 1.0" << std::endl;  // DO NOT REMOVE THIS LINE!!!
    std::cout << "===================================" << std::endl;
    std::cout << "               CV Task 1           " << std::endl;
    std::cout << "===================================" << std::endl;

    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <config-file>" << std::endl;
        return 1;
    }

    std::string path = std::string(argv[1]);
    std::ifstream fs(path);
    if (!fs)
    {
        std::cout << "Error: Failed to open file " << path << std::endl;
        return 2;
    }
    std::stringstream buffer;
    buffer << fs.rdbuf();

    rapidjson::Document doc;
    rapidjson::ParseResult check;
    check = doc.Parse<0>(buffer.str().c_str());

    if (check)
    {
        if (doc.HasMember("testcases"))
        {
            rapidjson::Value &testcases = doc["testcases"];
            for (rapidjson::SizeType i = 0; i < testcases.Size(); i++)
            {
                rapidjson::Value &testcase = testcases[i];
                try
                {
                    execute_testcase(testcase);
                }
                catch (const std::exception &e)
                {
                    std::cout << e.what() << std::endl;
                    std::cout << BOLD(FRED("[ERROR]")) << " Program exited with errors!" << std::endl;
                    return -1;
                }
            }
        }
        std::cout << "Program exited normally!" << std::endl;
    }
    else
    {
        std::cout << "Error: Failed to parse file " << argv[1] << ":" << check.Offset() << std::endl;
        return 3;
    }
    return 0;
}
