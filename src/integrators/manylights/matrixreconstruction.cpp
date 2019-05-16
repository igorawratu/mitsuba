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
#include <thread>
#include <random>

#include "common.h"

MTS_NAMESPACE_BEGIN

MatrixReconstructionRenderer::MatrixReconstructionRenderer(std::unique_ptr<ManyLightsClusterer> clusterer,
    float sample_percentage, float min_dist, float step_size_factor, float tolerance, float tau, 
    std::uint32_t max_iterations, std::uint32_t slice_size, bool visibility_only, bool adaptive_col, 
    bool adaptive_importance_sampling, bool adaptive_force_resample, bool adaptive_recover_transpose,
    bool truncated, bool show_slices, bool vsl) : 
        clusterer_(std::move(clusterer)), 
        sample_percentage_(sample_percentage), 
        min_dist_(min_dist), 
        step_size_factor_(step_size_factor), 
        tolerance_(tolerance), 
        tau_(tau), 
        max_iterations_(max_iterations), 
        slice_size_(slice_size), 
        visibility_only_(visibility_only), 
        adaptive_col_sampling_(adaptive_col), 
        adaptive_importance_sampling_(adaptive_importance_sampling),
        adaptive_force_resample_(adaptive_force_resample),
        adaptive_recover_transpose_(adaptive_recover_transpose),
        truncated_(truncated),
        show_slices_(show_slices),
        vsl_(vsl),
        cancel_(false){
}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other) : 
    clusterer_(std::move(other.clusterer_)),
    sample_percentage_(other.sample_percentage_), 
    min_dist_(other.min_dist_), 
    step_size_factor_(other.step_size_factor_), 
    tolerance_(other.tolerance_), 
    tau_(other.tau_), 
    max_iterations_(other.max_iterations_), 
    slice_size_(other.slice_size_), 
    visibility_only_(other.visibility_only_), 
    adaptive_col_sampling_(other.adaptive_col_sampling_),
    adaptive_importance_sampling_(other.adaptive_importance_sampling_),
    adaptive_force_resample_(other.adaptive_force_resample_),
    truncated_(other.truncated_),
    show_slices_(other.show_slices_),
    vsl_(other.vsl_),
    cancel_(other.cancel_){
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
        visibility_only_ = other.visibility_only_;
        adaptive_col_sampling_ = other.adaptive_col_sampling_;
        adaptive_importance_sampling_ = other.adaptive_importance_sampling_;
        adaptive_force_resample_ = other.adaptive_force_resample_;
        adaptive_recover_transpose_ = other.adaptive_recover_transpose_;
        truncated_ = other.truncated_;
        show_slices_ = other.show_slices_;
        vsl_ = other.vsl_;
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::unique_ptr<KDTNode<ReconstructionSample>> constructKDTree(Scene* scene, std::uint32_t size_threshold, 
    std::vector<ReconstructionSample>& samples, float min_dist, bool calc_unoccluded_samples,
    const std::vector<VPL>& vpls, bool vsl){

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

            ReconstructionSample curr_sample;
            curr_sample.image_x = x;
            curr_sample.image_y = y;
            curr_sample.ray = ray;

            if(calc_unoccluded_samples){
                curr_sample.unoccluded_samples.resize(vpls.size());

                for(std::uint32_t i = 0; i < vpls.size(); ++i){
                    curr_sample.unoccluded_samples[i] = sample(scene, sampler, curr_sample.its, ray, vpls[i], min_dist, false, 
                        10, i == 0, curr_sample.intersected_scene, true, vsl);

                    if(!curr_sample.intersected_scene || curr_sample.its.isEmitter()){
                        break;
                    }
                }
            }
            else{
                //call to sample primarily to get intersection details
                sample(scene, sampler, curr_sample.its, ray, vpls[0], min_dist, false, 10, true, 
                    curr_sample.intersected_scene, true, true);
            }

            if(curr_sample.intersected_scene && curr_sample.its.isEmitter()){
                curr_sample.emitter_color = curr_sample.its.Le(-ray.d);
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
    const std::vector<VPL>& vpls, Eigen::MatrixXd& rmat, Eigen::MatrixXd& gmat, Eigen::MatrixXd& bmat,
    std::uint32_t num_samples, float min_dist, std::mt19937& rng, bool vsl){
    assert(rmat.rows() * rmat.cols() > 0 && gmat.rows() * gmat.cols() > 0 && bmat.rows() * bmat.cols() > 0);

    std::uint32_t total_samples = slice->sample_indices.size() * vpls.size();
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

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < indices_to_compute.size(); ++i){
        std::uint32_t light_index = indices_to_compute[i] % vpls.size();
        std::uint32_t sample_index = indices_to_compute[i] / vpls.size();
        const VPL& vpl = vpls[light_index];
        ReconstructionSample& sample_to_compute = slice->sample(sample_index);


        Spectrum lightContribution = sample(scene, sampler, sample_to_compute.its, sample_to_compute.ray, vpl, 
            min_dist, true, 10, false, sample_to_compute.intersected_scene, true, vsl);

        Float r, g, b;
        lightContribution.toLinearRGB(r, g, b);
        rmat(sample_index, light_index) = r;
        gmat(sample_index, light_index) = g;
        bmat(sample_index, light_index) = b;
    }

    return indices_to_compute;
}

void copyMatrixToBuffer(std::uint8_t* output_image, const Eigen::MatrixXd& rmat, const Eigen::MatrixXd& gmat,
    const Eigen::MatrixXd& bmat, KDTNode<ReconstructionSample>* slice, Vector2i image_size){
    for(std::uint32_t i = 0; i < rmat.rows(); ++i){
        float r = 0, g = 0, b = 0;
        if(slice->sample(i).its.isEmitter()){
            slice->sample(i).emitter_color.toSRGB(r, g, b);
        }
        else{
            for(int j = 0; j < rmat.cols(); ++j){
                r += rmat(i, j);
                g += gmat(i, j);
                b += bmat(i, j);
            }

            Spectrum converter;
            converter.fromLinearRGB(r, g, b);
            converter.toSRGB(r, g, b);
        }
        
        ReconstructionSample& sample = slice->sample(i);

        std::uint32_t buffer_pos = sample.image_x + sample.image_y * image_size.x;

        output_image[buffer_pos*3] = std::max(0.f, std::min(1.f, r)) * 255 + 0.5f;
        output_image[buffer_pos*3 + 1] = std::max(0.f, std::min(1.f, g)) * 255 + 0.5f;
        output_image[buffer_pos*3 + 2] = std::max(0.f, std::min(1.f, b)) * 255 + 0.5f;
    }
}

void copyMatrixToBuffer(std::uint8_t* output_image, const Eigen::MatrixXd& mat, KDTNode<ReconstructionSample>* slice, 
    Vector2i image_size, bool visibility_only, bool recover_transpose){

    std::uint32_t total_samples = slice->sample_indices.size();
    for(std::uint32_t i = 0; i < total_samples; ++i){
        float r = 0, g = 0, b = 0;
        if(slice->sample(i).its.isEmitter()){
            slice->sample(i).emitter_color.toSRGB(r, g, b);
        }
        else{
            if(visibility_only){
                Spectrum light_contributions(0.f);
                for(std::uint32_t j = 0; j < slice->sample(i).unoccluded_samples.size(); ++j){
                    float coeff = /*mat(i, j);*/mat(i, j) > 0.f ? 1.f : 0.f;
                    light_contributions += slice->sample(i).unoccluded_samples[j] * coeff;
                }

                light_contributions.toSRGB(r, g, b);
            }
            else{
                std::uint32_t lights = recover_transpose ? mat.rows() / 3 : mat.cols();
                for(std::uint32_t j = 0; j < lights; ++j){
                    if(recover_transpose){
                        r += mat(j * 3, i);
                        g += mat(j * 3 + 1, i);
                        b += mat(j * 3 + 2, i);
                    }
                    else{
                        r += mat(i * 3, j);
                        g += mat(i * 3 + 1, j);
                        b += mat(i * 3 + 2, j);
                    }
                }

                Spectrum converter;
                converter.fromLinearRGB(r, g, b);
                converter.toSRGB(r, g, b);
            }
        }
        
        ReconstructionSample& sample = slice->sample(i);

        std::uint32_t buffer_pos = sample.image_x + sample.image_y * image_size.x;

        output_image[buffer_pos*3] = std::max(0.f, std::min(1.f, r)) * 255 + 0.5f;
        output_image[buffer_pos*3 + 1] = std::max(0.f, std::min(1.f, g)) * 255 + 0.5f;
        output_image[buffer_pos*3 + 2] = std::max(0.f, std::min(1.f, b)) * 255 + 0.5f;
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

std::vector<std::uint32_t> importanceSample(KDTNode<ReconstructionSample>* slice, std::uint32_t sample_index, 
    std::uint32_t num_samples, std::mt19937& rng, bool lights){
    assert((lights && slice->sample(sample_index).unoccluded_samples.size() > 0) ||
        (!lights && slice->sample_indices.size()) > 0);
    
    std::uint32_t total_vals = lights ? slice->sample(sample_index).unoccluded_samples.size() : 
        slice->sample_indices.size();
    std::vector<std::uint32_t> available(total_vals);
    std::iota(available.begin(), available.end(), 0);

    num_samples = std::min((size_t)num_samples, available.size());
    if(num_samples == available.size()){
        return available;
    }

    std::vector<std::uint32_t> sampled(num_samples);
    for(std::uint32_t i = 0; i < num_samples; ++i){
        double total_contrib = 0.;
        for(std::uint32_t j = 0; j < available.size(); ++j){
            total_contrib += lights ? slice->sample(sample_index).unoccluded_samples[available[j]].getLuminance() :
                slice->sample(available[j]).unoccluded_samples[sample_index].getLuminance();
        }

        std::uniform_real_distribution<double> gen(0., total_contrib);
        double selection = gen(rng);

        std::uint32_t idx = 0;
        for(; idx < available.size(); ++idx){
            selection -= lights ? slice->sample(sample_index).unoccluded_samples[available[idx]].getLuminance() :
                slice->sample(available[idx]).unoccluded_samples[sample_index].getLuminance();
            if(selection <= 0.){
                break;
            }
        }

        sampled[i] = available[idx];
        available[idx] = available.back();
        available.pop_back();
    }

    std::sort(sampled.begin(), sampled.end());

    return sampled;
}

std::vector<std::uint32_t> sampleRow(Scene* scene, KDTNode<ReconstructionSample>* slice, const std::vector<VPL>& vpls, 
    std::uint32_t col, float min_dist, std::uint32_t num_samples, std::mt19937& rng, Eigen::MatrixXd& mat, 
    const std::vector<std::uint32_t>& sample_set, bool resample, bool visibility_only, bool recover_transpose,
    bool adaptive_sampling, bool vsl){
    
    std::uint32_t expected_row_length = visibility_only ? num_samples : num_samples * 3;
    std::uint32_t max_samples = recover_transpose ? vpls.size() : slice->sample_indices.size();

    assert((size_t)mat.rows() == expected_row_length && mat.cols() == 1 && 
        num_samples <= max_samples);

    std::vector<std::uint32_t> sampled_indices;
    if(resample){
        if(adaptive_sampling){
            sampled_indices = importanceSample(slice, col, num_samples, rng, recover_transpose);
        }
        else{
            if(num_samples == max_samples){
                sampled_indices.resize(num_samples);
                std::iota(sampled_indices.begin(), sampled_indices.end(), 0);
            }
            else{
                std::uint32_t samples_per_bucket = float(max_samples) / 10 + 0.5f;
                if(max_samples < 50 || (samples_per_bucket * 10) < num_samples){
                    std::vector<std::uint32_t> indices(max_samples);
                    std::iota(indices.begin(), indices.end(), 0);
                    std::random_shuffle(indices.begin(), indices.end());
                    sampled_indices.insert(sampled_indices.begin(), indices.begin(), indices.begin() + num_samples);
                    std::sort(sampled_indices.begin(), sampled_indices.end());
                }
                else{
                    std::array<std::vector<std::uint32_t>, 10> buckets;
                    for(std::uint8_t i = 0; i < 10; ++i){
                        std::uint32_t num_bucket_samples = i < 9 ? samples_per_bucket : 
                            max_samples - 9 * samples_per_bucket;
                        buckets[i] = std::vector<std::uint32_t>(num_bucket_samples);
                        std::iota(buckets[i].begin(), buckets[i].end(), samples_per_bucket * i);
                        std::random_shuffle(buckets[i].begin(), buckets[i].end());
                    }

                    sampled_indices.resize(num_samples);
                    std::uint32_t pos = 0;

                    while(pos < num_samples){
                        std::uint32_t bucket = pos % 10;
                        std::uint32_t level = pos / 10;

                        sampled_indices[pos++] = buckets[bucket][level];
                        //std::cout << buckets[bucket][level] << " " << max_samples << std::endl;
                    }
                }
            }  
        }
    }
    else{
        sampled_indices = sample_set;
    }

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < sampled_indices.size(); ++i){
        const VPL& vpl = recover_transpose ? vpls[sampled_indices[i]] : vpls[col];
        ReconstructionSample& scene_sample = recover_transpose ? slice->sample(col) : slice->sample(sampled_indices[i]);

        if(visibility_only){
            Point ray_origin = scene_sample.its.p;
            Ray shadow_ray(ray_origin, normalize(vpl.its.p - ray_origin), 0.f);

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            mat(i, 0) = 1.;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if(abs((ray_origin - vpl.its.p).length() - t) > 0.0001f ){
                    mat(i, 0) = -1.f;
                }
            }
        }
        else{
            Spectrum lightContribution = sample(scene, sampler, scene_sample.its, scene_sample.ray, vpl, 
                min_dist, true, 10, false, scene_sample.intersected_scene, true, vsl);

            Float r, g, b;
            lightContribution.toLinearRGB(r, g, b);
            mat(i * 3, 0) = r;
            mat(i * 3 + 1, 0) = g;
            mat(i * 3 + 2, 0) = b;
        }
    }

    return sampled_indices;
}

std::uint32_t adaptiveMatrixReconstruction(Eigen::MatrixXd& mat, Scene* scene,
    KDTNode<ReconstructionSample>* slice, const std::vector<VPL>& vpls, float min_dist, float sample_perc,
    std::mt19937& rng, bool visibility_only, bool recover_transpose, bool adaptive_sampling, bool force_resample,
    bool vsl){
    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);

    std::uint32_t total_row_samples = recover_transpose ? vpls.size() : slice->sample_indices.size();
    std::uint32_t expected_row_size = visibility_only ? total_row_samples : total_row_samples * 3;
    std::uint32_t num_cols = recover_transpose ? slice->sample_indices.size() : vpls.size();
    
    if((size_t)mat.cols() != num_cols || (size_t)mat.rows() != expected_row_size){
        mat = Eigen::MatrixXd(expected_row_size, num_cols);
    }

    std::uint32_t num_samples = total_row_samples * sample_perc + 0.5f;
    if(num_samples == 0){
        num_samples = total_row_samples;
    }
    std::uint32_t expected_omega_rows = visibility_only ? num_samples : num_samples * 3;

    std::vector<std::uint32_t> order(num_cols);
    std::iota(order.begin(), order.end(), 0);
    std::random_shuffle(order.begin(), order.end());

    Eigen::MatrixXd reconstructed(expected_row_size, 1);
    std::vector<std::uint32_t> sampled;

    Eigen::MatrixXd q;
    Eigen::MatrixXd sample_omega(expected_omega_rows, 1);
    Eigen::MatrixXd q_omega;

    for(std::uint32_t i = 0; i < num_cols; ++i){
        if(q.cols() > 0){
            sample_omega.setZero();
            //previous full sample was used, which means that there is now a need to regenerate the sample indices and 
            //the omega matrices
            if(num_samples != sampled.size() || num_samples == slice->sample_indices.size() || force_resample){
                sampled = sampleRow(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    true, visibility_only, recover_transpose, adaptive_sampling, vsl);
                q_omega.resize(expected_omega_rows, q.cols());

                q_omega.resize(sampled.size(), q.cols());

                for(std::uint32_t j = 0; j < sampled.size(); ++j){
                    if(visibility_only){
                        q_omega.row(j) = q.row(sampled[j]);
                    }
                    else{
                        q_omega.row(j * 3) = q.row(sampled[j] * 3);
                        q_omega.row(j * 3 + 1) = q.row(sampled[j] * 3 + 1);
                        q_omega.row(j * 3 + 2) = q.row(sampled[j] * 3 + 2);
                    }
                }
            }
                //no new direction was added so no need to regenerate sample indices
            else{
                sampled = sampleRow(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, 
                    sampled, false, visibility_only, recover_transpose, adaptive_sampling, vsl);
            }

            auto svd = q_omega.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
            auto sv = svd.singularValues();
            Eigen::MatrixXd singular_val_inv = Eigen::MatrixXd::Zero(sv.size(), sv.size());
            for(std::uint32_t i = 0; i < sv.size(); ++i){
                singular_val_inv(i, i) = sv(i) < 1e-10f ? 0.f : 1.f / sv(i);
            }

            Eigen::MatrixXd q_omega_pseudoInverse = svd.matrixV() * singular_val_inv * svd.matrixU().transpose();
            reconstructed = q * q_omega_pseudoInverse * sample_omega;
            
            float d = 0.f;
            for(std::uint32_t j = 0; j < sampled.size(); ++j){
                float dist = fabs(reconstructed(sampled[j], 0) - sample_omega(j, 0));
                /*if(dist > std::numeric_limits<float>::epsilon())*/{
                    d += dist;
                }
            }

            if(d > 1e-10){
                sampled = sampleRow(scene, slice, vpls, order[i], min_dist, total_row_samples, rng, 
                    reconstructed, sampled, true, visibility_only, recover_transpose, adaptive_sampling, vsl);

                q.conservativeResize(q.rows(), q.cols() + 1);
                q.col(q.cols() - 1) = reconstructed;
            }
        }
        else{
            sampled = sampleRow(scene, slice, vpls, order[i], min_dist, slice->sample_indices.size(), rng, reconstructed, 
                sampled, true, visibility_only, recover_transpose, adaptive_sampling, vsl);
            if(reconstructed.norm() > 0.f){
                q = reconstructed;
            }
        }

        mat.col(order[i]) = reconstructed.col(0);
    }

    return slice->sample_indices.size() * q.cols() + (mat.cols() - q.cols()) * num_samples;
}


void svt(Eigen::MatrixXd& reconstructed_matrix, const Eigen::MatrixXd& lighting_matrix, float step_size, 
    float tolerance, float tau, std::uint32_t max_iterations, const std::vector<std::uint32_t>& sampled_indices,
    bool truncated){

    std::uint32_t k0 = tau / (step_size * lighting_matrix.norm()) + 1.5f; //extra .5 for rounding in case of float error
    Eigen::MatrixXd y = step_size * (float)k0 * lighting_matrix;
    for(std::uint32_t i = 0; i < max_iterations; ++i){
        std::uint32_t max_possible_rank = std::min(y.cols(), y.rows());
        std::uint32_t initial_sv = std::max(1u, max_possible_rank / 10u);
        std::uint32_t increment = std::max(1u, max_possible_rank / 20u);

        if(truncated){
            reconstructed_matrix = softThreshRank(y, tau, initial_sv, increment);
        }
        else{
            reconstructed_matrix = softThreshRankNoTrunc(y, tau);
        }

        float numer_total_dist = 0.f;

        for(std::uint32_t j = 0; j < sampled_indices.size(); ++j){
            std::uint32_t row = sampled_indices[j] / lighting_matrix.cols();
            std::uint32_t col = sampled_indices[j] % lighting_matrix.cols();

            float d = reconstructed_matrix(row, col) - lighting_matrix(row, col);
            numer_total_dist += d * d;
        }

        //norm of lighting matrix is the same as the norm of its sampled entries since the other values are 0
        float ratio = sqrt(numer_total_dist) / lighting_matrix.norm();
        if(ratio < tolerance){
            break;
        }

        for(std::uint32_t j = 0; j < sampled_indices.size(); ++j){
            std::uint32_t row = sampled_indices[j] / lighting_matrix.cols();
            std::uint32_t col = sampled_indices[j] % lighting_matrix.cols();

            float step = lighting_matrix(row, col) - reconstructed_matrix(row, col);
            y(row, col) += step_size * step;
        }
    }
}

std::mutex mutex;

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
    memset(output_image, 0, output_bitmap->getBytesPerPixel() * size.x * size.y);

    std::cout << "constructing kd tree" << std::endl;
    auto kdt_root = constructKDTree(scene, slice_size_, samples_, min_dist_, true, vpls, vsl_);

    std::vector<KDTNode<ReconstructionSample>*> slices;
    getSlices(kdt_root.get(), slices);

    std::cout << "reconstructing slices" << std::endl;

    std::uint32_t amount_sampled = 0;
    std::uint32_t total_samples = 0;

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        int y = rand() % 255;
        int u = rand() % 255;
        int v = ((float)i / slices.size()) * 255.f;

        int sb = 1.164 * (y - 16) + 2.018 * (u - 128);
        int sg = 1.164 * (y - 16) - 0.813 * (v - 128) - 0.391 * (u - 128);
        int sr = 1.164 * (y - 16) + 1.596 * (v - 128);

        if(show_slices_){
            for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
                std::uint32_t offset = (slices[i]->sample(j).image_x + slices[i]->sample(j).image_y * 
                output_bitmap->getSize().x) * output_bitmap->getBytesPerPixel();

                output_image[offset] = sr;
                output_image[offset + 1] = sg;
                output_image[offset + 2] = sb;
            }

            continue;
        }
        else{
            std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * i);
        
            if(adaptive_col_sampling_){
                Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size() * 3, vpls.size());
                std::uint32_t samples = adaptiveMatrixReconstruction(mat, scene, slices[i], vpls, min_dist_, 
                    sample_percentage_, rng, visibility_only_, adaptive_recover_transpose_, 
                    adaptive_importance_sampling_, adaptive_force_resample_, vsl_);
                copyMatrixToBuffer(output_image, mat, slices[i], size, visibility_only_, adaptive_recover_transpose_);

                {
                    std::lock_guard<std::mutex> lock(mutex);
                    amount_sampled += samples;
                    total_samples += slices[i]->sample_indices.size() * vpls.size();
                }
            }
            else{
                Eigen::MatrixXd rmat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size(), vpls.size());
                Eigen::MatrixXd gmat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size(), vpls.size());
                Eigen::MatrixXd bmat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size(), vpls.size());

                std::uint32_t num_samples = slices[i]->sample_indices.size() * vpls.size() * sample_percentage_;
                auto indices = calculateSparseSamples(scene, slices[i], vpls, rmat, gmat, bmat, num_samples, min_dist_, rng, vsl_);

                float step_size = 1.9f;//(1.2f * lighting_matrix.rows() * lighting_matrix.cols()) / (indices.size() * 3.f); 

                Eigen::MatrixXd reconstructed_r, reconstructed_b, reconstructed_g;
                svt(reconstructed_r, rmat, step_size, tolerance_, tau_, max_iterations_, indices, truncated_);
                svt(reconstructed_g, gmat, step_size, tolerance_, tau_, max_iterations_, indices, truncated_);
                svt(reconstructed_b, bmat, step_size, tolerance_, tau_, max_iterations_, indices, truncated_);
                copyMatrixToBuffer(output_image, reconstructed_r, reconstructed_g, reconstructed_b, slices[i], size);
            }
        }

    }

    if(adaptive_col_sampling_){
        float sample_perc = (float)amount_sampled / total_samples;
        std::cout << "Sample percentage: " << sample_perc << std::endl;
    }

    film->setBitmap(output_bitmap);

    return cancel_;
}

MTS_NAMESPACE_END