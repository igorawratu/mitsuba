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
#include <tuple>
#include <queue>

#include "common.h"

MTS_NAMESPACE_BEGIN

std::vector<float> calculateClusterContributions(const std::vector<VPL>& vpls, KDTNode<ReconstructionSample>* slice, Scene* scene,
    float min_dist, std::uint32_t num_rows){
    num_rows = std::min(num_rows, std::uint32_t(slice->sample_indices.size()));
    std::vector<std::uint32_t> row_samples(slice->sample_indices.size());
    std::iota(row_samples.begin(), row_samples.end(), 0);
    std::random_shuffle(row_samples.begin(), row_samples.end());
    row_samples.resize(num_rows);

    std::vector<float> contributions(vpls.size(), 0.f);

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        for(std::uint32_t j = 0; j < row_samples.size(); ++j){
            auto& curr_sample = slice->sample(row_samples[j]);
            Spectrum c = sample(scene, sampler, curr_sample.its, curr_sample.ray, vpls[i], min_dist, true, 5, false, 
            curr_sample.intersected_scene, false, false);
            contributions[i] += c.getLuminance();
        }
    }

    return contributions;
}

std::vector<std::vector<std::uint32_t>> clusterVPLsBySampling(const std::vector<VPL>& vpls, std::uint32_t num_clusters,
    float min_dist){
    num_clusters = std::min(std::uint32_t(vpls.size()), num_clusters);
    std::vector<std::vector<std::uint32_t>> clusters;

    if(num_clusters == 0){
        return clusters;
    }

    std::vector<std::uint32_t> remaining_vpls(vpls.size());
    std::iota(remaining_vpls.begin(), remaining_vpls.end(), 0);

    //minus one because a cluster is added randomly at the start
    std::vector<std::pair<float, std::uint32_t>> closest_cluster_to_vpl(vpls.size() - 1, 
        std::make_pair(std::numeric_limits<float>::max(), 0));
    std::vector<float> vpl_dist_to_all_clusters(vpls.size() - 1, 0.f);

    //not anything too important so just using the lcg
    std::uint32_t idx = rand() % remaining_vpls.size();
    clusters.emplace_back(1, remaining_vpls[idx]);

    remaining_vpls[idx] = remaining_vpls.back();
    remaining_vpls.pop_back();

    while(clusters.size() < num_clusters){
        std::uint32_t cluster_to_add = 0;
        float max_dist = -std::numeric_limits<float>::max();

        //only need to update the distance info with the last
        for(std::uint32_t i = 0; i < remaining_vpls.size(); ++i){
            float dist = (vpls[remaining_vpls[i]].its.p - vpls[clusters.back()[0]].its.p).lengthSquared();
            dist += dot(vpls[remaining_vpls[i]].its.shFrame.n, vpls[clusters.back()[0]].its.shFrame.n) * min_dist;

            if(dist < closest_cluster_to_vpl[i].first){
                closest_cluster_to_vpl[i].first = dist;
                closest_cluster_to_vpl[i].second = clusters.size() - 1;
            }

            vpl_dist_to_all_clusters[i] += dist;

            if(vpl_dist_to_all_clusters[i] > max_dist){
                cluster_to_add = i;
                max_dist = vpl_dist_to_all_clusters[i];
            }
        }

        clusters.emplace_back(1, remaining_vpls[cluster_to_add]);
        
        remaining_vpls[cluster_to_add] = remaining_vpls.back();
        vpl_dist_to_all_clusters[cluster_to_add] = vpl_dist_to_all_clusters.back();
        closest_cluster_to_vpl[cluster_to_add] = closest_cluster_to_vpl.back();

        remaining_vpls.pop_back();
        vpl_dist_to_all_clusters.pop_back();
        closest_cluster_to_vpl.pop_back();
    }

    for(std::uint32_t i = 0; i < remaining_vpls.size(); ++i){
        clusters[closest_cluster_to_vpl[i].second].push_back(remaining_vpls[i]);
    }

    return clusters;
}

std::vector<std::vector<std::uint32_t>> clusterVPLsBySplitting(const std::vector<std::vector<std::uint32_t>>& sampling_clusters, 
    const std::vector<VPL>& vpls, const std::vector<float> contributions, std::uint32_t num_clusters, std::mt19937& rng){
    std::vector<std::vector<std::uint32_t>> splitting_clusters = sampling_clusters;

    typedef std::pair<std::uint32_t, float> ClusterContributionPair;

    //it's not possible to split size 1 clusters, so we put them at the back of the queue regardless of their contribution
    auto cluster_comparator = [&splitting_clusters](const ClusterContributionPair& l, const ClusterContributionPair& r){
        if(splitting_clusters[l.first].size() == 1){
            return true;
        }

        if(splitting_clusters[r.first].size() == 1){
            return false;
        }

        return l.second < r.second;
    };

    std::priority_queue<ClusterContributionPair, std::vector<ClusterContributionPair>, 
        decltype(cluster_comparator)> split_queue(cluster_comparator);

    for(std::uint32_t i = 0; i < splitting_clusters.size(); ++i){
        float cluster_contrib = 0.f;
        for(std::uint32_t j = 0; j < splitting_clusters[i].size(); ++j){
            cluster_contrib += splitting_clusters[i][j];
        }

        split_queue.push(std::make_pair(i, cluster_contrib));
    }

    auto vpl_proj_comp = [](const std::pair<std::uint32_t, float>& l, const std::pair<std::uint32_t, float>& r){
        return l.second < r.second;
    };

    std::uniform_real_distribution<float> gen_pick(0.f, 1.f);

    std::uint32_t total_clusters = splitting_clusters.size() + num_clusters;
    while(splitting_clusters.size() < total_clusters){
        auto largest = split_queue.top();
        split_queue.pop();

        //generate random line and then project VPLs onto it. Then split it according to half the amount of energy along this line
        Eigen::VectorXf rand_line(6);
        do{
            for(std::uint32_t i = 0; i < 6; ++i){
                rand_line(i) = gen_pick(rng);
            }
        }while(rand_line.norm() == 0.f);

        std::vector<std::pair<std::uint32_t, float>> vpl_projection_distances;

        std::vector<std::uint32_t>& curr_cluster = splitting_clusters[largest.first];
        for(std::uint32_t i = 0; i < curr_cluster.size(); ++i){
            Eigen::VectorXf vpl_point(6);
            for(std::uint32_t j = 0; j < 3; ++j){
                vpl_point(j) = vpls[curr_cluster[i]].its.p[j];
            }

            vpl_point(3) = vpls[curr_cluster[i]].its.shFrame.n.x;
            vpl_point(4) = vpls[curr_cluster[i]].its.shFrame.n.y;
            vpl_point(5) = vpls[curr_cluster[i]].its.shFrame.n.z;

            float d_on_line = vpl_point.dot(rand_line) / rand_line.dot(rand_line);

            vpl_projection_distances.push_back(std::make_pair(curr_cluster[i], d_on_line));
        }

        std::sort(vpl_projection_distances.begin(), vpl_projection_distances.end(), vpl_proj_comp);
        std::uint32_t split_index = 1;
        
        if(largest.second < std::numeric_limits<float>::epsilon()){
            split_index = vpl_projection_distances.size() / 2;
        }
        else{
            float half = largest.second / 2.f;
            for(std::uint32_t i = 0; i < vpl_projection_distances.size(); ++i){
                half -= contributions[vpl_projection_distances[i].second];
                if(half <= 0.f){
                    split_index = i;
                    break;
                }
            }
        }

        float new_cluster1_contrib = 0.f;
        float new_cluster2_contrib = 0.f;
        std::vector<std::uint32_t> new_cluster1;
        std::vector<std::uint32_t> new_cluster2;

        for(std::uint32_t i = 0; i < vpl_projection_distances.size(); ++i){
            if(i < split_index){
                new_cluster1.push_back(vpl_projection_distances[i].first);
                new_cluster1_contrib += contributions[new_cluster1.back()];
            }
            else{
                new_cluster2.push_back(vpl_projection_distances[i].first);
                new_cluster2_contrib += contributions[new_cluster2.back()];
            }
        }

        splitting_clusters[largest.first] = new_cluster1;
        splitting_clusters.push_back(new_cluster2);

        split_queue.push(std::make_pair(largest.first, new_cluster1_contrib));
        split_queue.push(std::make_pair(splitting_clusters.size() - 1, new_cluster2_contrib));
    }

    //return splitting_clusters;

    return sampling_clusters;
}

std::vector<VPL> sampleRepresentatives(const std::vector<VPL>& vpls, const std::vector<std::vector<std::uint32_t>>& clusters, 
    std::mt19937& rng){
    
    std::vector<VPL> representatives;

    //might need to make this more deterministic
    for(std::uint32_t i = 0; i < clusters.size(); ++i){
        std::uniform_int_distribution<std::uint32_t> gen(0, clusters[i].size() - 1);
        std::uint32_t rep_idx = gen(rng);

        float tot_cluster_power = 0.f;
        for(std::uint32_t j = 0; j < clusters[i].size(); ++j){
            tot_cluster_power += vpls[clusters[i][j]].P.getLuminance();
        }

        representatives.push_back(vpls[clusters[i][rep_idx]]);
        representatives.back().P = representatives.back().P / representatives.back().P.getLuminance() * tot_cluster_power;
    }

    return representatives;
}

std::unique_ptr<KDTNode<ReconstructionSample>> constructKDTree(Scene* scene, std::uint32_t size_threshold, 
    std::vector<ReconstructionSample>& samples){

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

            ReconstructionSample curr_sample;
            curr_sample.image_x = x;
            curr_sample.image_y = y;

            std::uint32_t num_bounces = 0;

            while(true){
                curr_sample.intersected_scene = scene->rayIntersect(ray, curr_sample.its);
                if(!curr_sample.intersected_scene){
                    break;
                }

                if(curr_sample.its.getBSDF()->getType() & BSDF::ESmooth || curr_sample.its.isEmitter()){
                    break;
                }

                if(++num_bounces > 5){
                    break;
                }

                BSDFSamplingRecord bsdf_sample_record(curr_sample.its, sampler);
                bsdf_sample_record.typeMask = BSDF::EReflection;
                curr_sample.its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());

                curr_sample.ray = Ray(curr_sample.its.p, bsdf_sample_record.its.toWorld(bsdf_sample_record.wo), ray.time);
            }

            if(!curr_sample.intersected_scene){
                if(scene->hasEnvironmentEmitter()){
                    curr_sample.emitter_color = scene->evalEnvironment(RayDifferential(curr_sample.ray));
                }
                else curr_sample.emitter_color = Spectrum(0.f);
            }
            else if(curr_sample.its.isEmitter()){
                curr_sample.emitter_color = curr_sample.its.Le(-curr_sample.ray.d);
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

void calculateUnoccludedSamples(KDTNode<ReconstructionSample>* slice, Scene* scene, float min_dist, 
    const std::vector<VPL>& vpls, Sampler* sampler){
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        auto& curr_sample = slice->sample(i);
        curr_sample.unoccluded_samples.resize(vpls.size());
        for(std::uint32_t j = 0; j < vpls.size(); ++j){
            Spectrum col = sample(scene, sampler, curr_sample.its, curr_sample.ray, vpls[i], min_dist, false, 
                5, false, curr_sample.intersected_scene, false, false);

            curr_sample.unoccluded_samples[j] = col;
        }
    }
}

//Used in the proximal gradient descent version of recovery. Matrix is sparsely populated uniformly with observations
//RGB assumed, which are dealt with in separate matrices
std::vector<std::uint32_t> calculateSparseSamples(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, Eigen::MatrixXd& rmat, Eigen::MatrixXd& gmat, Eigen::MatrixXd& bmat,
    std::uint32_t num_samples, float min_dist, bool vsl){
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

        Spectrum lightContribution = sample(scene, sampler, sample_to_compute.its, sample_to_compute.ray, vpl, min_dist, true,
            5, false, sample_to_compute.intersected_scene, false, vsl);

        Float r, g, b;
        lightContribution.toLinearRGB(r, g, b);
        rmat(sample_index, light_index) = r;
        gmat(sample_index, light_index) = g;
        bmat(sample_index, light_index) = b;
    }

    return indices;
}

void updateSliceWithMatData(const Eigen::MatrixXd& rmat, const Eigen::MatrixXd& gmat,
    const Eigen::MatrixXd& bmat, KDTNode<ReconstructionSample>* slice){
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        for(std::uint32_t j = 0; j < slice->sample(i).unoccluded_samples.size(); ++j){
            slice->sample(i).unoccluded_samples[j].fromLinearRGB(rmat(i, j), gmat(i, j), bmat(i, j));
        }
    }
}

void updateSliceWithMatData(const Eigen::MatrixXd& mat, KDTNode<ReconstructionSample>* slice, 
    bool visibility_only, bool recover_transpose, bool vsl, Scene* scene, const std::vector<VPL>& vpls,
    float min_dist){
    
    Sampler *sampler = nullptr;

    if(vsl){
        Properties props("independent");
        sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
        sampler->configure();
        sampler->generate(Point2i(0));
    }
    
    std::uint32_t total_samples = slice->sample_indices.size();
    for(std::uint32_t i = 0; i < total_samples; ++i){
        if(visibility_only){
            Spectrum light_contributions(0.f);
            for(std::uint32_t j = 0; j < slice->sample(i).unoccluded_samples.size(); ++j){
                float coeff = ((recover_transpose ? mat(j, i) : mat(i, j)) + 1.f) / 2.f;
                if(vsl){
                    if(coeff > std::numeric_limits<float>::epsilon()){
                        slice->sample(i).unoccluded_samples[j] = sample(scene, sampler, slice->sample(i).its, slice->sample(i).ray,
                            vpls[j], min_dist, false, 5, false, slice->sample(i).intersected_scene, false, true) * coeff;
                    }
                    else{
                        slice->sample(i).unoccluded_samples[j] = Spectrum(0.f);
                    }
                }
                else{
                    slice->sample(i).unoccluded_samples[j] *= coeff;
                }
            }
        }
        else{
            for(std::uint32_t j = 0; j < slice->sample(i).unoccluded_samples.size(); ++j){
                float r = recover_transpose ? mat(j * 3, i) : mat(i * 3, j);
                float g = recover_transpose ? mat(j * 3 + 1, i) : mat(i * 3 + 1, j);
                float b = recover_transpose ? mat(j * 3 + 2, i) : mat(i * 3 + 2, j);

                slice->sample(i).unoccluded_samples[j].fromLinearRGB(r, g, b);
            }
        }
    }
}

void copySamplesToBuffer(std::uint8_t* output_image, const std::vector<ReconstructionSample>& samples, Vector2i image_size){
    #pragma omp parallel for
    for(std::uint32_t i = 0; i < samples.size(); ++i){
        Spectrum col(0.f);
        if(samples[i].intersected_scene){
            for(std::uint32_t j = 0; j < samples[i].unoccluded_samples.size(); ++j){
                col += samples[i].unoccluded_samples[j];
            }

            if(samples[i].its.isEmitter()){
                col += samples[i].emitter_color;
            }
        }
        else{
            //environment map colour is also stored in emitter_color
            col = samples[i].emitter_color;
        }

        std::uint32_t buffer_pos = samples[i].image_x + samples[i].image_y * image_size.x;

        float r, g, b;
        col.toSRGB(r, g, b);

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
    bool importance_sample, const std::vector<float>& probabilities, bool vsl){
    
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
            Spectrum lightContribution = sample(scene, sampler, scene_sample.its, scene_sample.ray, vpl, 
                min_dist, true, 5, false, scene_sample.intersected_scene, false, vsl);

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
    std::uint32_t& basis_rank, bool vsl){

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
                    true, visibility_only, recover_transpose, importance_sample, probabilities, vsl);

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
                    false, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
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
                    reconstructed, sampled, true, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
                samples_for_col = total_rows;

                q.conservativeResize(q.rows(), q.cols() + 1);
                q.col(q.cols() - 1) = reconstructed;

                full_col_sampled = true;
            }
        }
        else{
            sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, reconstructed, 
                sampled, true, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
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

MatrixReconstructionRenderer::MatrixReconstructionRenderer(const std::vector<VPL>& vpls,
    float sample_percentage, float min_dist, float step_size_factor, float tolerance, float tau, 
    std::uint32_t max_iterations, std::uint32_t slice_size, bool visibility_only, bool adaptive_col, 
    bool adaptive_importance_sampling, bool adaptive_force_resample, bool adaptive_recover_transpose,
    bool truncated, bool show_slices, bool vsl) : 
        vpls_(vpls), 
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

        std::uint32_t num_clusters = std::min(std::uint32_t(vpls_.size()), std::uint32_t(1000 * 2.f / 3.f));
        clusters_by_sampling_ = clusterVPLsBySampling(vpls_, num_clusters, min_dist_);
}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other) : 
    vpls_(other.vpls_),
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
    cancel_(other.cancel_),
    clusters_by_sampling_(other.clusters_by_sampling_){

    
}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (MatrixReconstructionRenderer&& other){
    if(this != &other){
        vpls_ = other.vpls_;
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
        clusters_by_sampling_ = other.clusters_by_sampling_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::mutex mutex;

bool MatrixReconstructionRenderer::render(Scene* scene){
    if(scene == nullptr){
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
    auto kdt_root = constructKDTree(scene, slice_size_, samples_);
    std::cout << "getting slices" << std::endl;

    std::vector<KDTNode<ReconstructionSample>*> slices;
    getSlices(kdt_root.get(), slices);

    std::cout << "reconstructing slices" << std::endl;

    std::uint32_t amount_sampled = 0;
    std::uint32_t total_samples = 0;

    if(show_slices_){
        #pragma omp parallel for
        for(std::uint32_t i = 0; i < slices.size(); ++i){
            int y = rand() % 255;
            int u = rand() % 255;
            int v = ((float)i / slices.size()) * 255.f;

            int sb = 1.164 * (y - 16) + 2.018 * (u - 128);
            int sg = 1.164 * (y - 16) - 0.813 * (v - 128) - 0.391 * (u - 128);
            int sr = 1.164 * (y - 16) + 1.596 * (v - 128);

            for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
                std::uint32_t offset = (slices[i]->sample(j).image_x + slices[i]->sample(j).image_y * 
                output_bitmap->getSize().x) * output_bitmap->getBytesPerPixel();

                output_image[offset] = sr;
                output_image[offset + 1] = sg;
                output_image[offset + 2] = sb;
            }
        }
    }
    else{
        std::uint32_t split_num_clusters = std::min(vpls_.size(), 1000 - clusters_by_sampling_.size());

        std::vector<std::vector<float>> vpl_contributions(slices.size());
        std::vector<std::vector<std::vector<std::uint32_t>>> clusters_by_splitting(slices.size());
        std::vector<std::vector<VPL>> vpls(slices.size());

        Properties props("independent");
        Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
        sampler->configure();
        sampler->generate(Point2i(0));

        std::mt19937 cluster_rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());

        #pragma omp parallel for
        for(std::uint32_t i = 0; i < slices.size(); ++i){
            vpl_contributions[i] = calculateClusterContributions(vpls_, slices[i], scene, min_dist_, 30);
            clusters_by_splitting[i] = clusterVPLsBySplitting(clusters_by_sampling_, vpls_, vpl_contributions[i],
                    split_num_clusters, cluster_rng);
            vpls[i] = sampleRepresentatives(vpls_, clusters_by_splitting[i], cluster_rng);
            calculateUnoccludedSamples(slices[i], scene, min_dist_, vpls[i], sampler);
        }

        #pragma omp parallel for
        for(std::uint32_t i = 0; i < slices.size(); ++i){
            std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * i);

            if(adaptive_col_sampling_){
                std::uint32_t mat_rows = visibility_only_ ? slices[i]->sample_indices.size() : slices[i]->sample_indices.size() * 3;
                Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(mat_rows, vpls[i].size());
                std::uint32_t max_rank = std::min(mat.rows(), mat.cols());
                std::uint32_t basis_rank;
                std::uint32_t samples = adaptiveMatrixReconstruction(mat, scene, slices[i], vpls[i], min_dist_, 
                    sample_percentage_, rng, visibility_only_, adaptive_recover_transpose_, 
                    adaptive_importance_sampling_, adaptive_force_resample_, basis_rank, vsl_);
                
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

                updateSliceWithMatData(mat, slices[i], visibility_only_, adaptive_recover_transpose_, vsl_, scene, 
                    vpls[i], min_dist_);

                {
                    std::lock_guard<std::mutex> lock(mutex);
                    amount_sampled += samples;
                    total_samples += slices[i]->sample_indices.size() * vpls[i].size();
                }
            }
            else{
                Eigen::MatrixXd rmat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size(), vpls[i].size());
                Eigen::MatrixXd gmat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size(), vpls[i].size());
                Eigen::MatrixXd bmat = Eigen::MatrixXd::Zero(slices[i]->sample_indices.size(), vpls[i].size());

                std::uint32_t num_samples = slices[i]->sample_indices.size() * vpls[i].size() * sample_percentage_;
                
                auto indices = calculateSparseSamples(scene, slices[i], vpls[i], rmat, gmat, bmat, num_samples, min_dist_, vsl_);

                float step_size = 1.9f;//(1.2f * lighting_matrix.rows() * lighting_matrix.cols()) / (indices.size() * 3.f); 

                Eigen::MatrixXd reconstructed_r, reconstructed_b, reconstructed_g;
                svt(reconstructed_r, rmat, step_size, tolerance_, tau_, max_iterations_, indices, truncated_);
                svt(reconstructed_g, gmat, step_size, tolerance_, tau_, max_iterations_, indices, truncated_);
                svt(reconstructed_b, bmat, step_size, tolerance_, tau_, max_iterations_, indices, truncated_);
                updateSliceWithMatData(reconstructed_r, reconstructed_g, reconstructed_b, slices[i]);
            }
        }

        copySamplesToBuffer(output_image, samples_, size);

        if(adaptive_col_sampling_){
            float sample_perc = (float)amount_sampled / total_samples;
            std::cout << "Sample percentage: " << sample_perc << std::endl;
        }
    }
    
    film->setBitmap(output_bitmap);

    return cancel_;
}

MTS_NAMESPACE_END