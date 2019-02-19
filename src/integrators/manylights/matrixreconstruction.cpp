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
#include "common.h"

MTS_NAMESPACE_BEGIN

MatrixReconstructionRenderer::MatrixReconstructionRenderer(std::unique_ptr<ManyLightsClusterer> clusterer,
    float sample_percentage, float min_dist, float step_size_factor, float tolerance, float tau, 
    std::uint32_t max_iterations, std::uint32_t slice_size) : clusterer_(std::move(clusterer)), 
    sample_percentage_(sample_percentage), min_dist_(min_dist), step_size_factor_(step_size_factor), 
    tolerance_(tolerance), tau_(tau), max_iterations_(max_iterations), slice_size_(slice_size), cancel_(false){
}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other) : 
    clusterer_(std::move(other.clusterer_)), sample_percentage_(other.sample_percentage_), min_dist_(other.min_dist_), 
    step_size_factor_(other.step_size_factor_), tolerance_(other.tolerance_), tau_(other.tau_), max_iterations_(other.max_iterations_), 
    slice_size_(other.slice_size_), cancel_(other.cancel_){
}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (MatrixReconstructionRenderer&& other){
    if(this != &other){
        clusterer_ = std::move(other.clusterer_);
        sample_percentage_ = other.sample_percentage_;
        min_dist_ = other.min_dist_;
        step_size_factor_ = other.step_size_factor_;
        tolerance_ = other.tolerance_; 
        tau_ = other.tau_;
        max_iterations_ = other.max_iterations_;
        slice_size_ = other.slice_size_;
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::unique_ptr<KDTNode<ReconstructionSample>> constructKDTree(Scene* scene, std::uint32_t size_threshold, 
    std::vector<ReconstructionSample>& samples, float min_dist){

    auto kdt_root = std::unique_ptr<KDTNode<ReconstructionSample>>(new KDTNode<ReconstructionSample>(&samples));

    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    samples.resize(film->getSize().y * film->getSize().x);

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for (std::int32_t y = 0; y < film->getSize().y; ++y) {
        #pragma omp parallel for
        for (std::int32_t x = 0; x < film->getSize().x; ++x) {
            Ray ray;

            Point2 sample_position(x + 0.5f, y + 0.5f);

            Point2 aperture_sample(0.5f, 0.5f);
            Float time_sample(0.5f);

            sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

            Intersection its;

            bool intersected = scene->rayIntersect(ray, its);
            if(!intersected){
                continue;
            }

            ReconstructionSample curr_sample(x, y, intersected, its);;

            if(its.isEmitter()){
                curr_sample.emitter_color = its.Le(-ray.d);
            }

            samples[y * film->getSize().x + x] = std::move(curr_sample);
        }
    }

    kdt_root->sample_indices.reserve(samples.size());
    for(std::uint32_t i = 0; i < samples.size(); ++i){
        if(samples[i].intersected_scene){
            kdt_root->sample_indices.push_back(i);
        }
    }

    splitKDTree(kdt_root.get(), size_threshold, min_dist);

    return kdt_root;
}

std::vector<std::uint32_t> calculateSparseSamples(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, Eigen::MatrixXf& matrix, std::uint32_t num_samples, float min_dist, 
    std::mt19937& rng){
    assert(matrix.rows() > 0 && matrix.cols() > 0);

    std::uint32_t total_samples = matrix.rows() * matrix.cols() / 3;
    num_samples = std::min(num_samples, total_samples);
    std::vector<std::uint32_t> indices(total_samples);
    std::iota(indices.begin(), indices.end(), 0);

    std::vector<std::uint32_t> indices_to_compute(num_samples);

    for(std::uint32_t i = 0; i < indices_to_compute.size(); ++i){
        std::uniform_int_distribution<std::uint32_t> gen(0, indices.size() - 1);

        std::uint32_t idx = gen(rng);
        indices_to_compute[i] = indices[idx];
        indices[idx] = indices.back();
        indices.pop_back();
    }

    ref<Sensor> sensor = scene->getSensor();
	ref<Film> film = sensor->getFilm();

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < indices_to_compute.size(); ++i){
        std::uint32_t light_index = indices_to_compute[i] % vpls.size();
        std::uint32_t sample_index = indices_to_compute[i] / vpls.size();
        const VPL& vpl = vpls[light_index];
        ReconstructionSample& sample_to_compute = slice->sample(sample_index);

        Spectrum lightContribution = sample(scene, sampler, sample_to_compute.its, vpl, min_dist, true);

        Float r, g, b;
        lightContribution.toSRGB(r, g, b);
        matrix(sample_index, light_index * 3) = r;
        matrix(sample_index, light_index * 3 + 1) = g;
        matrix(sample_index, light_index * 3 + 2) = b;
    }

    return indices_to_compute;
}

void copyMatrixToBuffer(std::uint8_t* output_image, Eigen::MatrixXf& light_matrix, KDTNode<ReconstructionSample>* slice, 
    Vector2i image_size){
    for(std::uint32_t i = 0; i < light_matrix.rows(); ++i){
        float r = 0, g = 0, b = 0;
        for(int j = 0; j < light_matrix.cols() / 3; ++j){
            r += light_matrix(i, j * 3);
            g += light_matrix(i, j * 3 + 1);
            b += light_matrix(i, j * 3 + 2);
        }

        ReconstructionSample& sample = slice->sample(i);

        std::uint32_t buffer_pos = sample.image_x + sample.image_y * image_size.y;

        output_image[buffer_pos*3] = std::min(1.f, r) * 255 + 0.5f;
        output_image[buffer_pos*3 + 1] = std::min(1.f, g) * 255 + 0.5f;
        output_image[buffer_pos*3 + 2] = std::min(1.f, b) * 255 + 0.5f;
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
    float tolerance, float tau, std::uint32_t max_iterations, const std::vector<std::uint32_t>& sampled_indices,
    std::mutex& mutex){

    std::uint32_t k0 = tau / (step_size * lighting_matrix.norm()) + 1.5f; //extra .5 for rounding in case of float error
    Eigen::MatrixXf y = step_size * (float)k0 * lighting_matrix;

    for(std::uint32_t i = 0; i < max_iterations; ++i){
        reconstructed_matrix = softThreshRank(y, tau, std::min(y.cols(), y.rows()) / 10, 
            std::min(y.cols(), y.rows()) / 20, mutex);

        float numer_total_dist = 0.f;

        for(std::uint32_t j = 0; j < sampled_indices.size(); ++j){
            for(std::uint32_t k = 0; k < 3; ++k){
                std::uint32_t row = sampled_indices[j] / (lighting_matrix.cols() / 3);
                std::uint32_t col = sampled_indices[j] % (lighting_matrix.cols() / 3) * 3 + k;

                float d = reconstructed_matrix(row, col) - lighting_matrix(row, col);
                numer_total_dist += d * d;
            }
        }

        float ratio = sqrt(numer_total_dist) / lighting_matrix.norm();

        if(ratio < tolerance){
            break;
        }

        for(std::uint32_t j = 0; j < sampled_indices.size(); ++j){
            for(std::uint32_t k = 0; k < 3; ++k){
                std::uint32_t row = sampled_indices[j] / (lighting_matrix.cols() / 3);
                std::uint32_t col = sampled_indices[j] % (lighting_matrix.cols() / 3) * 3 + k;

                float step = lighting_matrix(row, col) - reconstructed_matrix(row, col);
                y(row, col) += step_size * step;
            }
        }
    }
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

    std::cout << "constructing kd tree" << std::endl;
    auto kdt_root = constructKDTree(scene, slice_size_, samples_, min_dist_);

    std::cout << "getting slices" << std::endl;
    std::vector<KDTNode<ReconstructionSample>*> slices;
    getSlices(kdt_root.get(), slices);

    float step_size = 1.5f;//step_size_factor * (float)(lighting_matrix.rows() * lighting_matrix.cols()) / 
        //(float)(indices_to_compute.size() * 3); 

    std::cout << "reconstructing slices" << std::endl;
    std::mutex mutex;
    #pragma omp parallel for
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        Eigen::MatrixXf lighting_matrix = Eigen::MatrixXf::Zero(slices[i]->sample_indices.size(), vpls.size() * 3);
        std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * i);

        std::uint32_t num_samples = slices[i]->sample_indices.size() * vpls.size() * sample_percentage_;
        auto indices = calculateSparseSamples(scene, slices[i], vpls, lighting_matrix, num_samples, min_dist_, rng);

        Eigen::MatrixXf reconstructed_matrix;
        svt(reconstructed_matrix, lighting_matrix, step_size, tolerance_, tau_, max_iterations_, indices, mutex);

        copyMatrixToBuffer(output_image, reconstructed_matrix, slices[i], size);
    }

    film->setBitmap(output_bitmap);

    return cancel_;
}

MTS_NAMESPACE_END