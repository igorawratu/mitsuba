#include "matrixreconstruction.h"

#include <mitsuba/core/plugin.h>
#include <chrono>
#include <algorithm>
#include <set>
#include <fstream>
#include <string>
#include "definitions.h"
#include <utility>
#include <eigen3/Eigen/Dense>

MTS_NAMESPACE_BEGIN

MatrixReconstructionRenderer::MatrixReconstructionRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, std::pair<std::uint32_t, std::uint32_t> bucket_size,
    std::uint32_t light_samples, float min_dist, float step_size_factor, float tolerance, float tau, 
    std::uint32_t max_iterations) : clusterer_(std::move(clusterer)), bucket_size_(bucket_size), light_samples_(light_samples), min_dist_(min_dist),
    step_size_factor_(step_size_factor), tolerance_(tolerance), tau_(tau), max_iterations_(max_iterations), cancel_(false){
}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other) : clusterer_(std::move(other.clusterer_)), 
    bucket_size_(other.bucket_size_), light_samples_(other.light_samples_), min_dist_(other.min_dist_), 
    step_size_factor_(other.step_size_factor_), tolerance_(other.tolerance_), tau_(other.tau_), max_iterations_(other.max_iterations_), cancel_(other.cancel_){
}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (MatrixReconstructionRenderer&& other){
    if(this != &other){
        clusterer_ = std::move(other.clusterer_);
        bucket_size_ = other.bucket_size_; 
        light_samples_ = other.light_samples_;
        min_dist_ = other.min_dist_;
        step_size_factor_ = other.step_size_factor_;
        tolerance_ = other.tolerance_; 
        tau_ = other.tau_;
        max_iterations_ = other.max_iterations_;
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::vector<std::pair<std::uint32_t, std::uint32_t>> generateComputableIndices(
    const std::pair<std::uint32_t, std::uint32_t>& bucket_size, const Vector2i& image_size, std::uint32_t num_lights, 
    std::uint32_t light_samples){

    std::uint32_t x_buckets = image_size.x / bucket_size.first + image_size.x % (bucket_size.first > 0 ? 1 : 0);
    std::uint32_t y_buckets = image_size.y / bucket_size.second + image_size.y % (bucket_size.second > 0 ? 1 : 0);
    light_samples = std::min(light_samples, num_lights);

    std::vector<std::pair<std::uint32_t, std::uint32_t>> indices_to_compute;
    indices_to_compute.reserve(x_buckets * y_buckets * light_samples);

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    for(std::uint32_t i = 0; i < y_buckets; ++i){
        for(std::uint32_t j = 0; j < x_buckets; ++j){
            std::uint32_t x_start = j * bucket_size.first;
            std::uint32_t y_start = i * bucket_size.second;

            std::uniform_int_distribution<std::uint32_t> gen_y(0, std::min(bucket_size.second, image_size.y - y_start) - 1);
            std::uniform_int_distribution<std::uint32_t> gen_x(0, std::min(bucket_size.first, image_size.x - x_start) - 1);
            std::uniform_int_distribution<std::uint32_t> gen_light(0, num_lights - 1);

            std::uint32_t max_tries = num_lights * 10;
            std::uint32_t current_tries = 0;
            std::set<std::pair<std::uint32_t, std::uint32_t>> light_indices;
            
            while(current_tries++ < max_tries && light_indices.size() < light_samples){
                light_indices.insert(std::make_pair(
                    (gen_x(rng) + x_start) + (gen_y(rng) + y_start) * image_size.x, gen_light(rng)));
            }

            for(auto iter = light_indices.begin(); iter != light_indices.end(); ++iter){
                indices_to_compute.emplace_back(*iter);
            }
        }
    }

    return indices_to_compute;
}

void calculateSparseSamples(Scene* scene, const std::vector<VPL>& vpls, Eigen::MatrixXf& matrix,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& indices, const Vector2i& size, float min_dist){
    
    ref<Sensor> sensor = scene->getSensor();
	ref<Film> film = sensor->getFilm();

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
    Point2 aperture_sample(0.5f, 0.5f);
    Float time_sample(0.5f);

    for(size_t i = 0; i < indices.size(); ++i){
        std::uint32_t x = indices[i].first % size.x;
        std::uint32_t y = indices[i].first / size.x;
        const VPL& vpl = vpls[indices[i].second];

        Ray ray;

        Point2 sample_position(x + 0.5f, y + 0.5f);

        sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            continue;
        }

        Normal n = its.geoFrame.n;

        Point ray_origin = its.p;
        Ray shadow_ray(ray_origin, normalize(vpl.its.p - ray_origin), ray.time);

        Float t;
        ConstShapePtr shape;
        Normal norm;
        Point2 uv;
        if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
            if(abs((ray_origin - vpl.its.p).length() - t) > 0.0001f ){
                continue;
            }
        }

        BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
        bsdf_sample_record.wi = its.toLocal(normalize(vpl.its.p - its.p));
        bsdf_sample_record.wo = its.toLocal(n);

        Spectrum albedo = its.getBSDF()->eval(bsdf_sample_record);

        float d = std::max((its.p - vpl.its.p).length(), min_dist);
        float attenuation = 1.f / (d * d);

        float n_dot_ldir = std::max(0.f, dot(normalize(n), normalize(vpl.its.p - its.p)));
        float ln_dot_ldir = std::max(0.f, dot(normalize(vpl.its.shFrame.n), normalize(its.p - vpl.its.p)));

        Spectrum lightContribution = (vpl.P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
        Float r, g, b;
        lightContribution.toSRGB(r, g, b);
        matrix(indices[i].first * 3, indices[i].second) = r;
        matrix(indices[i].first * 3 + 1, indices[i].second) = g;
        matrix(indices[i].first * 3 + 2, indices[i].second) = b;
    }
}

void copyMatrixToBuffer(std::uint8_t* output_image, Eigen::MatrixXf& light_matrix, Vector2i image_size){
    for(std::uint32_t i = 0; i < light_matrix.rows() / 3; ++i){
        float r = 0, g = 0, b = 0;
        for(int j = 0; j < light_matrix.cols(); ++j){
            r += light_matrix(i*3, j);
            g += light_matrix(i*3 + 1, j);
            b += light_matrix(i*3 + 2, j);
        }

        output_image[i*3] = std::min(1.f, r) * 255 + 0.5f;
        output_image[i*3 + 1] = std::min(1.f, g) * 255 + 0.5f;
        output_image[i*3 + 2] = std::min(1.f, b) * 255 + 0.5f;
    }
}

void printToFile(const std::vector<float>& vals, std::string filename, std::ios_base::openmode mode, bool new_line, 
    std::streamsize precision = 0){

    std::ofstream output_file;
    output_file.open(filename, mode);

    if(precision > 0){
        output_file.precision(precision);
    }

    for(std::uint32_t i = 0; i < vals.size(); ++i){
        output_file << vals[i] << " ";
        if(new_line) 
            output_file << std::endl;
    }

    if(!new_line) output_file << std::endl;

    output_file.close();
}

void svt(Eigen::MatrixXf& reconstructed_matrix, const Eigen::MatrixXf& lighting_matrix, float step_size, 
    float tolerance, float tau, std::uint32_t max_iterations, 
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& sampled_indices){

    std::uint32_t k0 = tau / (step_size * lighting_matrix.norm()) + 1.5f; //extra .5 for rounding in case of float error
    Eigen::MatrixXf y = step_size * (float)k0 * lighting_matrix;

    std::vector<float> error;
    std::vector<float> trace;

    for(std::uint32_t i = 0; i < max_iterations; ++i){
        std::cout << "Computing svd for iteration " << i << std::endl;
        auto svd = y.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);

        std::cout << "Reconstructing matrix for iteration " << i << std::endl;
        auto singular_values = svd.singularValues();
        Eigen::MatrixXf diagonal_singular = Eigen::MatrixXf::Zero(svd.matrixU().rows(), svd.matrixV().rows());
        
        for(std::uint32_t j = 0; j < singular_values.rows(); ++j){
            diagonal_singular(j, j) = std::max(0.f, singular_values(j, 0) - tau);
        }

        reconstructed_matrix = svd.matrixU() * diagonal_singular * svd.matrixV().transpose();

        auto reconsvd = reconstructed_matrix.jacobiSvd();
        auto sv = reconsvd.singularValues();
        
        float tn = 0.f;
        for(std::uint32_t j = 0; j < sv.rows(); ++j){
            tn += sv(j, 0);
        }
        trace.push_back(tn);

        float numer_total_dist = 0.f;

        for(std::uint32_t j = 0; j < sampled_indices.size(); ++j){
            for(std::uint32_t k = 0; k < 3; ++k){
                std::uint32_t row = sampled_indices[j].first * 3 + k;
                std::uint32_t col = sampled_indices[j].second;

                float d = reconstructed_matrix(row, col) - lighting_matrix(row, col);
                numer_total_dist += d * d;
            }
        }

        float ratio = sqrt(numer_total_dist) / lighting_matrix.norm();

        error.push_back(ratio);

        std::cout << ratio << std::endl;
        if(ratio < tolerance){
            break;
        }

        for(std::uint32_t j = 0; j < sampled_indices.size(); ++j){
            for(std::uint32_t k = 0; k < 3; ++k){
                std::uint32_t row = sampled_indices[j].first * 3 + k;
                std::uint32_t col = sampled_indices[j].second;

                float step = lighting_matrix(row, col) - reconstructed_matrix(row, col);
                y(row, col) += step_size * step;
            }
        }
    }

    printToFile(trace, "trace.txt", std::ios::out, true);
    printToFile(error, "error.txt", std::ios::out, true);
}

bool MatrixReconstructionRenderer::render(Scene* scene){
    Intersection its;
    auto vpls = clusterer_->getClusteringForPoint(its);

    if(scene == nullptr || vpls.size() == 0){
        return true;
    }

    {
        std::lock_guard<std::mutex> lock(cancel_lock_);
        cancel_ = false;
    }

    ref<Sensor> sensor = scene->getSensor();
	ref<Film> film = sensor->getFilm();

    auto size = film->getSize();
    if(size.x == 0 || size.y == 0){
        return true;
    }

    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, size);
    std::uint8_t* output_image = output_bitmap->getUInt8Data();

    Eigen::MatrixXf lighting_matrix = Eigen::MatrixXf::Zero(size.x * size.y * 3, vpls.size());

    auto indices_to_compute = generateComputableIndices(bucket_size_, size, vpls.size(), light_samples_);
    calculateSparseSamples(scene, vpls, lighting_matrix, indices_to_compute, size, min_dist_);

    auto svd = lighting_matrix.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto singular_values = svd.singularValues();
    Eigen::MatrixXf diagonal_singular = Eigen::MatrixXf::Zero(svd.matrixU().rows(), svd.matrixV().rows());

    std::vector<float> sv_to_write;

    for(std::uint32_t j = 0; j < singular_values.rows(); ++j){
        diagonal_singular(j, j) = singular_values(j, 0);
        sv_to_write.push_back(singular_values(j, 0));
    }

    printToFile(sv_to_write, "originalsv.txt", std::ios::out, true);

    lighting_matrix = svd.matrixU() * diagonal_singular * svd.matrixV().transpose();
    
    Eigen::MatrixXf reconstructed_matrix;

    float step_size = 1.5f;//step_size_factor * (float)(lighting_matrix.rows() * lighting_matrix.cols()) / 
        //(float)(indices_to_compute.size() * 3); 

    svt(reconstructed_matrix, lighting_matrix, step_size, tolerance_, tau_, max_iterations_, indices_to_compute);

    svd = reconstructed_matrix.jacobiSvd();
    auto rsv = svd.singularValues();
    sv_to_write.clear();

    for(std::uint32_t j = 0; j < rsv.rows(); ++j){
        sv_to_write.push_back(rsv(j, 0));
    }

    printToFile(sv_to_write, "reconstructedsv.txt", std::ios::out, true);

    copyMatrixToBuffer(output_image, reconstructed_matrix, size);
    
    printToFile(std::vector<float>(), "initial_matrix.txt", std::ios::out, true);
    printToFile(std::vector<float>(), "reconstructed_matrix.txt", std::ios::out, true);
    printToFile(std::vector<float>(), "full_matrix.txt", std::ios::out, true);

    Eigen::MatrixXf full_matrix = Eigen::MatrixXf::Zero(size.x * size.y * 3, vpls.size());
    auto full_indices = generateComputableIndices(std::make_pair(1, 1), size, vpls.size(), vpls.size());
    calculateSparseSamples(scene, vpls, full_matrix, full_indices, size, min_dist_);

    svd = full_matrix.jacobiSvd();
    auto fsv = svd.singularValues();
    sv_to_write.clear();

    for(std::uint32_t j = 0; j < fsv.rows(); ++j){
        sv_to_write.push_back(fsv(j, 0));
    }

    printToFile(sv_to_write, "fullsv.txt", std::ios::out, true);

    std::vector<float> init_row(lighting_matrix.cols());
    std::vector<float> recon_row(lighting_matrix.cols());
    std::vector<float> full_row(lighting_matrix.cols());

    for(std::uint32_t i = 0; i < lighting_matrix.rows(); ++i){
        for(std::uint32_t j = 0; j < lighting_matrix.cols(); ++j){
            init_row[j] = lighting_matrix(i, j);
            recon_row[j] = reconstructed_matrix(i, j);
            full_row[j] = full_matrix(i, j);
        }
        printToFile(init_row, "initial_matrix.txt", std::ios::out | std::ios::app, false, 2);
        printToFile(recon_row, "reconstructed_matrix.txt", std::ios::out | std::ios::app, false, 2);
        printToFile(full_row, "full_matrix.txt", std::ios::out | std::ios::app, false, 2);
    }

    film->setBitmap(output_bitmap);

    return cancel_;
}

MTS_NAMESPACE_END