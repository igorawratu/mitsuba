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
#include "RedSVD-h.hpp"

#include "common.h"

MTS_NAMESPACE_BEGIN

MatrixReconstructionRenderer::MatrixReconstructionRenderer(std::unique_ptr<ManyLightsClusterer> clusterer,
    float sample_percentage, float min_dist, float step_size_factor, float tolerance, float tau, 
    std::uint32_t max_iterations, std::uint32_t slice_size, bool visibility_only, bool adaptive_col, 
    bool adaptive_importance_sampling, bool adaptive_force_resample, bool adaptive_recover_transpose,
    bool truncated, bool show_slices) : 
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
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::unique_ptr<KDTNode<ReconstructionSample>> constructKDTree(Scene* scene, std::uint32_t size_threshold, 
    std::vector<ReconstructionSample>& samples, float min_dist, const std::vector<VPL>& vpls){

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

            //calculates the position the ray intersects with
            Ray ray;

            Point2 sample_position(x + 0.5f, y + 0.5f);

            Point2 aperture_sample(0.5f, 0.5f);
            Float time_sample(0.5f);

            sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

            Intersection its;

            bool intersected = scene->rayIntersect(ray, its);
            ReconstructionSample curr_sample(x, y, intersected, its, ray);

            if(intersected){
                //unoccluded colours calculated for heuristic purposes. a proper calculation is performed for the moment but it might
                //be better if something less accurate is used
                curr_sample.unoccluded_samples.resize(vpls.size());
                float total_lum = 0.f;

                for(std::uint32_t i = 0; i < vpls.size(); ++i){
                    curr_sample.unoccluded_samples[i] = sample(scene, sampler, its, ray, vpls[i], min_dist, false);
                    total_lum += curr_sample.unoccluded_samples[i].getLuminance();
                }

                if(its.isEmitter()){
                    curr_sample.emitter_color = its.Le(-ray.d);
                }
            }
            

            samples[y * film->getSize().x + x] = std::move(curr_sample);
        }
    }

    //only samples that intersected the scene are recovered. The rest are sampled from the environment map if there is one
    kdt_root->sample_indices.reserve(samples.size());
    for(std::uint32_t i = 0; i < samples.size(); ++i){
        if(samples[i].intersected_scene){
            kdt_root->sample_indices.push_back(i);
        }
    }

    splitKDTree(kdt_root.get(), size_threshold, 0.f);

    return kdt_root;
}

//Used in the proximal gradient descent version of recovery. Matrix is sparsely populated uniformly with observations
//RGB assumed, which are dealt with in separate matrices
std::vector<std::uint32_t> calculateSparseSamples(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, Eigen::MatrixXd& rmat, Eigen::MatrixXd& gmat, Eigen::MatrixXd& bmat,
    std::uint32_t num_samples, float min_dist){
    assert(rmat.rows() * rmat.cols() > 0 && gmat.rows() * gmat.cols() > 0 && bmat.rows() * bmat.cols() > 0);

    std::uint32_t total_samples = slice->sample_indices.size() * vpls.size();
    num_samples = std::min(num_samples, total_samples);
    std::vector<std::uint32_t> indices(total_samples);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    indices.erase(indices.begin() + num_samples, indices.end());

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < indices.size(); ++i){
        std::uint32_t light_index = indices[i] % vpls.size();
        std::uint32_t sample_index = indices[i] / vpls.size();
        const VPL& vpl = vpls[light_index];
        ReconstructionSample& sample_to_compute = slice->sample(sample_index);

        Spectrum lightContribution = sample(scene, sampler, sample_to_compute.its, sample_to_compute.ray, vpl, min_dist, true);

        Float r, g, b;
        lightContribution.toLinearRGB(r, g, b);
        rmat(sample_index, light_index) = r;
        gmat(sample_index, light_index) = g;
        bmat(sample_index, light_index) = b;
    }

    return indices;
}

void copyMatrixToBuffer(std::uint8_t* output_image, const Eigen::MatrixXd& rmat, const Eigen::MatrixXd& gmat,
    const Eigen::MatrixXd& bmat, KDTNode<ReconstructionSample>* slice, Vector2i image_size){
    for(std::uint32_t i = 0; i < rmat.rows(); ++i){
        float r = 0, g = 0, b = 0;

        for(int j = 0; j < rmat.cols(); ++j){
            r += rmat(i, j);
            g += gmat(i, j);
            b += bmat(i, j);
        }

        Spectrum converter;
        converter.fromLinearRGB(r, g, b);
        if(slice->sample(i).its.isEmitter()){
            converter += slice->sample(i).emitter_color;
        }

        //have to convert to SRGB because we are setting the buffer directly
        converter.toSRGB(r, g, b);
        
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
        if(visibility_only){
            Spectrum light_contributions(0.f);
            for(std::uint32_t j = 0; j < slice->sample(i).unoccluded_samples.size(); ++j){
                float coeff = ((recover_transpose ? mat(j, i) : mat(i, j)) + 1.f) / 2.f;
                light_contributions += slice->sample(i).unoccluded_samples[j] * coeff;
            }

            if(slice->sample(i).its.isEmitter()){
                light_contributions += slice->sample(i).emitter_color;
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

            if(slice->sample(i).its.isEmitter()){
                converter += slice->sample(i).emitter_color;
            }

            converter.toSRGB(r, g, b);
        }
        
        ReconstructionSample& sample = slice->sample(i);

        std::uint32_t buffer_pos = sample.image_x + sample.image_y * image_size.x;

        output_image[buffer_pos*3] = std::max(0.f, std::min(1.f, r)) * 255 + 0.5f;
        output_image[buffer_pos*3 + 1] = std::max(0.f, std::min(1.f, g)) * 255 + 0.5f;
        output_image[buffer_pos*3 + 2] = std::max(0.f, std::min(1.f, b)) * 255 + 0.5f;
    }
}

//helper for some stats
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

//generates a set of indices based on the probabilities passed through the vector
std::vector<std::uint32_t> importanceSample(std::uint32_t num_samples, std::mt19937& rng, const std::vector<float>& probabilities){
    std::vector<std::uint32_t> available(probabilities.size());
    std::iota(available.begin(), available.end(), 0);

    num_samples = std::min((size_t)num_samples, available.size());
    if(num_samples == available.size()){
        return available;
    }

    std::vector<std::uint32_t> sampled(num_samples);
    for(std::uint32_t i = 0; i < num_samples; ++i){
        double total_contrib = 0.;
        for(std::uint32_t j = 0; j < available.size(); ++j){
            total_contrib += probabilities[available[j]];
        }

        std::uniform_real_distribution<double> gen(0., total_contrib);
        double selection = gen(rng);

        std::uint32_t idx = 0;
        for(; idx < available.size(); ++idx){
            selection -= probabilities[available[idx]];
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

//sparsely samples a single column for adaptive matrix recovery
std::vector<std::uint32_t> sampleCol(Scene* scene, KDTNode<ReconstructionSample>* slice, const std::vector<VPL>& vpls, 
    std::uint32_t col, float min_dist, std::uint32_t num_samples, std::mt19937& rng, Eigen::MatrixXd& mat, 
    const std::vector<std::uint32_t>& sample_set, bool resample, bool visibility_only, bool recover_transpose,
    bool importance_sample, const std::vector<float>& probabilities){
    
    std::uint32_t num_rows = visibility_only ? num_samples : num_samples * 3;
    std::uint32_t max_samples = recover_transpose ? vpls.size() : slice->sample_indices.size();

    assert((size_t)mat.rows() == num_rows && mat.cols() == 1 && num_samples <= max_samples);

    std::vector<std::uint32_t> sampled_indices;

    if(resample){
        if(importance_sample){
            sampled_indices = importanceSample(num_samples, rng, probabilities);
        }
        else{
            if(num_samples == max_samples){
                sampled_indices.resize(num_samples);
                std::iota(sampled_indices.begin(), sampled_indices.end(), 0);
            }
            else{
                sampled_indices.resize(max_samples);
                std::iota(sampled_indices.begin(), sampled_indices.end(), 0);
                std::random_shuffle(sampled_indices.begin(), sampled_indices.end());
                sampled_indices.erase(sampled_indices.begin() + num_samples, sampled_indices.end());
                std::sort(sampled_indices.begin(), sampled_indices.end());
            }  
        }
    }
    else{
        //The previous set of indices used which is passed through as a parameter
        //If there is no need to regenerate the indices, then this is just used
        sampled_indices = sample_set;
    }

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < sampled_indices.size(); ++i){
        std::uint32_t vpl_index = recover_transpose ? sampled_indices[i] : col;
        std::uint32_t sample_index = recover_transpose ? col : sampled_indices[i];
        const VPL& vpl = vpls[vpl_index];
        ReconstructionSample& scene_sample = slice->sample(sample_index);

        //should centralize the shadow-cast to common actually
        if(visibility_only){
            Point ray_origin = scene_sample.its.p;
            Ray shadow_ray(ray_origin, normalize(vpl.its.p - ray_origin), 0.f);

            mat(i, 0) = 1.;

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if((ray_origin - vpl.its.p).length() - t > std::numeric_limits<float>::epsilon() * min_dist * 10.f){
                    mat(i, 0) = -1.;
                }
            }
        }
        else{
            Spectrum lightContribution = sample(scene, sampler, scene_sample.its, scene_sample.ray, vpl, min_dist, true);

            Float r, g, b;
            lightContribution.toLinearRGB(r, g, b);
            mat(i * 3, 0) = r;
            mat(i * 3 + 1, 0) = g;
            mat(i * 3 + 2, 0) = b;
        }
    }

    return sampled_indices;
}

std::uint32_t adaptiveMatrixReconstruction(Eigen::MatrixXd& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc,
    std::mt19937& rng, bool visibility_only, bool recover_transpose, bool importance_sample, bool force_resample,
    std::uint32_t& basis_rank){

    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);
    std::random_shuffle(slice->sample_indices.begin(), slice->sample_indices.end());

    std::uint32_t total_rows = recover_transpose ? vpls.size() : slice->sample_indices.size();
    std::uint32_t num_rows = visibility_only ? total_rows : total_rows * 3;
    std::uint32_t num_cols = recover_transpose ? slice->sample_indices.size() : vpls.size();
    
    //re-allocate matrix if it is of the incorrect size
    if((size_t)mat.cols() != num_cols || (size_t)mat.rows() != num_rows){
        mat = Eigen::MatrixXd(num_rows, num_cols);
    }
    mat.setZero();

    std::uint32_t num_samples = total_rows * sample_perc + 0.5f;
    //just in case, this shouldn't ever really happen
    if(num_samples == 0){
        num_samples = total_rows;
    }

    //calculate the contributions of each column to the rows. In the sense of lights as rows, this would be how much each light
    //contributes to the slice, whereas when the pixels are the rows, it would be the brightest pixel
    float max_contrib = -std::numeric_limits<float>::max();
    float total_contrib = 0.f;

    std::vector<float> col_max_contrib(num_cols, -std::numeric_limits<float>::max());
    std::vector<float> col_total_contrib(num_cols, 0.f);
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        for(std::uint32_t j = 0; j < slice->sample(i).unoccluded_samples.size(); ++j){
            float lum = slice->sample(i).unoccluded_samples[j].getLuminance();
            max_contrib = std::max(max_contrib, lum);
            total_contrib += lum;

            if(recover_transpose){
                col_max_contrib[i] = std::max(col_max_contrib[i], lum);
                col_total_contrib[i] += lum;
            }
            else{
                col_max_contrib[j] = std::max(col_max_contrib[j], lum);
                col_total_contrib[j] += lum;
            }
        }
    }

    //we recover from the brightest to least bright because there will be more full samples initially, allow for better coverage of
    //higher energy sections
    std::vector<std::uint32_t> order(num_cols);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), 
        [&col_total_contrib](const std::uint32_t& lhs, const std::uint32_t& rhs){
            return col_total_contrib[lhs] > col_total_contrib[rhs];
        });


    Eigen::MatrixXd reconstructed(num_rows, 1);
    std::vector<std::uint32_t> sampled;
    Eigen::MatrixXd q;
    Eigen::MatrixXd sample_omega;
    Eigen::MatrixXd q_omega;
    Eigen::MatrixXd q_omega_pseudoinverse;

    std::uint32_t total_samples = 0;
    std::uniform_real_distribution<float> gen(0.f, 1.f);
    std::vector<float> probabilities(num_rows, 1.f);
    bool resampled_required = false;

    //The actual adaptive matrix recovery algorithm
    for(std::uint32_t i = 0; i < order.size(); ++i){
        bool full_col_sampled = false;
        std::uint32_t samples_for_col = 0;

        //if the basis is not empty, we can try reproject, otherwise a full sample is required to populate the basis
        if(q.cols() > 0){
            std::uint32_t expected_omega_rows = visibility_only ? num_samples : num_samples * 3;
            sample_omega.resize(expected_omega_rows, 1);
            sample_omega.setZero();

            //we may want to regenerate the sample indices for a variety of reasons, in which case the indices are generated and
            //the pseudoinverse is recalculated
            if(num_samples != sampled.size() || num_samples == total_rows || force_resample || resampled_required){
                resampled_required = false;
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    true, visibility_only, recover_transpose, importance_sample, probabilities);

                q_omega.resize(expected_omega_rows, q.cols());

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

                auto svd = q_omega.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto sv = svd.singularValues();
                Eigen::MatrixXd singular_val_inv = Eigen::MatrixXd::Zero(sv.size(), sv.size());

                for(std::uint32_t j = 0; j < sv.size(); ++j){
                    singular_val_inv(j, j) = sv(j) < 1e-10 ? 0. : 1. / sv(j);
                }

                q_omega_pseudoinverse = svd.matrixV() * singular_val_inv * svd.matrixU().transpose();
            }
            //no new direction was added so no need to regenerate sample indices
            else{
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    false, visibility_only, recover_transpose, importance_sample, probabilities);
            }
            samples_for_col = num_samples;

            //no need to reconstruct if full sample
            reconstructed = sample_omega.rows() == reconstructed.rows() ? sample_omega : q * q_omega_pseudoinverse * sample_omega;

            double d = 0;
            for(std::uint32_t j = 0; j < sampled.size(); ++j){
                d += fabs(reconstructed(sampled[j], 0) - sample_omega(j, 0));
            }

            //this is just some heuristic on when to regenerate indices, need to investigate this better
            if(visibility_only){
                for(std::uint32_t j = 0; j < reconstructed.rows(); ++j){
                    if(fabs(fabs(reconstructed(j, 0)) - 1.) > std::numeric_limits<float>::epsilon()){
                        resampled_required = true;
                        break;
                    }
                }
            }

            //sampled values can't be reconstructed accurately so fully sample
            if(d > 1e-5){
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, 
                    reconstructed, sampled, true, visibility_only, recover_transpose, importance_sample, probabilities);
                samples_for_col = total_rows;

                q.conservativeResize(q.rows(), q.cols() + 1);
                q.col(q.cols() - 1) = reconstructed;

                full_col_sampled = true;
            }
        }
        else{
            sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, reconstructed, 
                sampled, true, visibility_only, recover_transpose, importance_sample, probabilities);
            samples_for_col = total_rows;

            //only add to basis if vector isn't a zero vector, as it is otherwise meaningless
            //might want to change this to be above some epsilon instead
            if(reconstructed.norm() > 0.f){
                q = reconstructed;
            }

            full_col_sampled = true;
        }

        //probability update for importance sampling. This needs to be revisited since this doesn't work well
        if(full_col_sampled){
            for(std::uint32_t j = 0; j < reconstructed.rows(); ++j){
                int next_idx = std::min(int(j) + 1, int(reconstructed.rows()) - 1);
                int prev_idx = std::max(int(j) - 1, 0);
                float g1 = std::abs(reconstructed(next_idx, 0) - reconstructed(j, 0));
                float g2 = std::abs(reconstructed(j, 0) - reconstructed(prev_idx, 0));
                float grad = (g1 + g2) / 2.f;

                float mul = float(order.size() - i) / order.size();
                mul = 1.f;//pow(mul, 2.f);

                probabilities[j] += std::min(grad, 1.f) * mul;
            }
        }
        
        
        mat.col(order[i]) = reconstructed.col(0);
        total_samples += samples_for_col;
    }

    basis_rank = q.cols();

    return total_samples;
}

//singular value thresholding algorithm for the non-adaptive version. Can use the RedSVD truncated svd, or just normal svd
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
    auto kdt_root = constructKDTree(scene, slice_size_, samples_, min_dist_, vpls);
    std::cout << "getting slices" << std::endl;
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
                std::uint32_t mat_rows = visibility_only_ ? slices[i]->sample_indices.size() : slices[i]->sample_indices.size() * 3;
                Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(mat_rows, vpls.size());
                std::uint32_t max_rank = std::min(mat.rows(), mat.cols());
                std::uint32_t basis_rank;
                std::uint32_t samples = adaptiveMatrixReconstruction(mat, scene, slices[i], vpls, min_dist_, 
                    sample_percentage_, rng, visibility_only_, adaptive_recover_transpose_, 
                    adaptive_importance_sampling_, adaptive_force_resample_, basis_rank);
                
                //float v = float(samples) / (mat.rows() * mat.cols());
                float v = float(basis_rank) / max_rank;
                float r, g, b;
                std::tie(r, g, b) = floatToRGB(v);
                std::uint8_t ro = r * 255.f;
                std::uint8_t go = g * 255.f;
                std::uint8_t bo = b * 255.f;

                for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
                    std::uint32_t offset = (slices[i]->sample(j).image_x + slices[i]->sample(j).image_y * 
                    output_bitmap->getSize().x) * output_bitmap->getBytesPerPixel();

                    output_image[offset] = ro;
                    output_image[offset + 1] = go;
                    output_image[offset + 2] = bo;
                }

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
                auto indices = calculateSparseSamples(scene, slices[i], vpls, rmat, gmat, bmat, num_samples, min_dist_);

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