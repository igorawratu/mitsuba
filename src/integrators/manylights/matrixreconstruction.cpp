#include "matrixreconstruction.h"

#include <mitsuba/core/plugin.h>
#include <chrono>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <fstream>
#include <string>
#include "definitions.h"
#include <utility>
#include <eigen3/Eigen/Dense>
#include <thread>
#include <random>
#include <tuple>
#include <queue>
#include <unordered_map>
#include <thread>
#include "lighttree.h"
#include <mitsuba/core/statistics.h>
#include <set>
#include <unistd.h>

#include "common.h"
#include "hwshader.h"
#include "blockingqueue.hpp"

MTS_NAMESPACE_BEGIN

//vpl radii need to be updated when clustered otherwise there will be bright spots again as the radii will be too small
void updateVPLRadii(std::vector<VPL>& vpls, float min_dist){
	std::uint32_t num_sl = 0;
	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(vpls[i].type != EDirectionalEmitterVPL){
			num_sl++;
		}
	} 

    //radius calculation, 11 to account 10 + 1 for adding a node's self in nearest neighbours
    std::uint32_t num_neighbours = std::min(num_sl, 11u);

    flann::Matrix<float> dataset(new float[num_sl * 3], num_sl, 3);
	std::uint32_t curr_light = 0;
	std::vector<std::uint32_t> pointlight_indices;
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(vpls[i].type != EDirectionalEmitterVPL){
			float* curr = (float*)dataset[curr_light++];
			curr[0] = vpls[i].its.p.x;
			curr[1] = vpls[i].its.p.y;
			curr[2] = vpls[i].its.p.z;
			pointlight_indices.emplace_back(i);
		}
		else{
			vpls[i].radius = 0.f;
		}
    }

    flann::Index<flann::L2<float>> index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();

	std::vector<std::vector<int>> indices;
	std::vector<std::vector<float>> distances;

    index.knnSearch(dataset, indices, distances, num_neighbours, flann::SearchParams(128));

	for(std::uint32_t i = 0; i < distances.size(); ++i){
		float max = std::numeric_limits<float>::min();
		float min = std::numeric_limits<float>::max();
		for(std::uint32_t j = 0; j < distances[i].size(); ++j){
			if(distances[i][j] > std::numeric_limits<float>::epsilon()){
				max = std::max(distances[i][j], max);
				min = std::min(distances[i][j], min);
			}
		}

		vpls[pointlight_indices[i]].radius = sqrt(max) * 2.f;
	}
}

//initial contribution matrix computation for lightslice
Eigen::MatrixXf calculateClusterContributions(const std::vector<VPL>& vpls, 
    const std::vector<KDTNode<ReconstructionSample>*>& slices, Scene* scene,
    float min_dist, std::uint32_t samples_per_slice, bool vsl){

    assert(slices.size() > 0 && vpls.size() > 0);
    Eigen::MatrixXf contributions = Eigen::MatrixXf::Zero(slices.size() * 3, vpls.size());
    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        for(std::uint32_t k = 0; k < samples_per_slice; ++k){
            auto& curr_sample = slices[i]->sample(rand() % slices[i]->sample_indices.size()); 
            for(std::uint32_t j = 0; j < vpls.size(); ++j){
                std::uint32_t num_samples;
                Spectrum c = sample(scene, sampler, curr_sample.its, curr_sample.ray, vpls[j], min_dist, true, 5, false, 
                    curr_sample.intersected_scene, false, vsl, num_samples);

                float r, g, b;
                c.toLinearRGB(r, g, b);

                contributions(i * 3, j) += r;
                contributions(i * 3 + 1, j) += g;
                contributions(i * 3 + 2, j) += b;
            }
        }
    }

    return contributions;
}

//cluster by sampling stage of lightslice
std::vector<std::vector<std::uint32_t>> clusterVPLsBySampling(const Eigen::MatrixXf& contributions,
    std::uint32_t num_clusters, const std::vector<VPL>& vpls, float min_dist){
    num_clusters = std::min(std::uint32_t(contributions.cols()), num_clusters);
    std::vector<std::vector<std::uint32_t>> clusters;

    if(num_clusters == 0){
        return clusters;
    }

    std::vector<float> col_norms(contributions.cols());
    std::vector<std::uint32_t> zero_norm_cols;
    std::vector<std::uint32_t> nonzero_norm_cols;

    for(std::uint32_t i = 0; i < contributions.cols(); ++i){
        col_norms[i] = contributions.col(i).norm();
        //change to eps for stability?
        if(col_norms[i] > 0.f){
            nonzero_norm_cols.push_back(i);
        }
        else{
            zero_norm_cols.push_back(i);
        }
    }

    Eigen::MatrixXf nz_norms(1, nonzero_norm_cols.size());
    Eigen::MatrixXf nz_cols(contributions.rows(), nonzero_norm_cols.size());
    std::vector<float> nz_norms_dist(nonzero_norm_cols.size());

    for(std::uint32_t i = 0; i < nonzero_norm_cols.size(); ++i){
        std::uint32_t& idx = nonzero_norm_cols[i];
        nz_norms_dist[i] = nz_norms(0, i) = col_norms[idx];
        nz_cols.col(i) = contributions.col(idx);
    }

    std::discrete_distribution<std::uint32_t> nzdist(nz_norms_dist.begin(), nz_norms_dist.end());
    auto nzprobs = nzdist.probabilities();
    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    std::unordered_map<std::uint32_t, float> cluster_centers;

    int centers_added = num_clusters;
    while(centers_added-- > 0){
        std::uint32_t selected = nzdist(rng);
        float pdf =  nzprobs[selected];
        if(pdf > 0.f){
            if(cluster_centers.find(selected) == cluster_centers.end()){
                cluster_centers[selected] = 1.f / pdf;
            }
            else{
                cluster_centers[selected] += 1.f / pdf;
            }
        }
    }

    Eigen::MatrixXf cluster_norms(1, cluster_centers.size());
    Eigen::MatrixXf cluster_cols(contributions.rows(), cluster_centers.size());

    std::uint32_t idx = 0;
    for(auto iter = cluster_centers.begin(); iter != cluster_centers.end(); ++iter, ++idx){
        cluster_norms(0, idx) = contributions.col(iter->first).norm() * iter->second;
        cluster_cols.col(idx) = contributions.col(iter->first) * iter->second;
    }

    Eigen::MatrixXf dist_mat = cluster_norms.transpose() * nz_norms - cluster_cols.transpose() * nz_cols;

    clusters.resize(cluster_cols.cols());
    for(std::uint32_t i = 0; i < dist_mat.cols(); ++i){
        float min = std::numeric_limits<float>::max();
        std::uint32_t min_idx = 0;
        for(std::uint32_t j = 0; j < dist_mat.rows(); ++j){
            if(min > dist_mat(j, i) && vpls[nonzero_norm_cols[i]].type == vpls[j].type){
                min = dist_mat(j, i);
                min_idx = j;
            }
        }

        clusters[min_idx].push_back(nonzero_norm_cols[i]);
    }

    if(zero_norm_cols.size() > 0){
        clusters.push_back(zero_norm_cols);
    }

    std::vector<std::vector<std::uint32_t>> non_empty_clusters;
    for(std::uint32_t i = 0; i < clusters.size(); ++i){
        if(clusters[i].size() > 0){
            non_empty_clusters.push_back(clusters[i]);
        }
    }

    return non_empty_clusters;
}

//lightslice cost computation
double computeCost(const Eigen::MatrixXf& contrib_mat, const std::vector<int>& nn, const std::vector<std::uint32_t>& cluster, 
    const std::vector<float>& norm_cache){
    if(cluster.size() < 2){
        return 0.;
    }

    double vsum = 0.f;
    double nsum = 0.f;
    for(std::uint32_t i = 0; i < cluster.size(); ++i){
        nsum += norm_cache[cluster[i]];
    }
    
    for(std::uint32_t i = 0; i < nn.size() * 3; ++i){
        double sum = 0.;

        for(std::uint32_t j = 0; j < cluster.size(); ++j){
            std::uint32_t nn_idx = i / 3;
            std::uint32_t nn_offset = i % 3;
            sum += contrib_mat(nn[nn_idx] * 3 + nn_offset, cluster[j]);
        }

        vsum += sum * sum;
    }
    
    return std::max(0., (nsum * nsum) - vsum);
}

//cluster by splitting stage of lightslice
std::vector<std::vector<std::uint32_t>> clusterVPLsBySplitting(const std::vector<std::vector<std::uint32_t>>& sampling_clusters, 
    const Eigen::MatrixXf& cluster_contributions, const std::vector<int>& nn, std::uint32_t num_clusters,
     std::mt19937& rng){

    std::vector<std::vector<std::uint32_t>> splitting_clusters = sampling_clusters;

    typedef std::pair<std::uint32_t, double> ClusterContributionPair;

    //it's not possible to split size 1 clusters, so we put them at the back of the queue regardless of their contribution
    auto cluster_comparator = [&splitting_clusters](const ClusterContributionPair& l, const ClusterContributionPair& r){
        if(splitting_clusters[l.first].size() <= 1){
            return true;
        }

        if(splitting_clusters[r.first].size() <= 1){
            return false;
        }

        return l.second < r.second;
    };

    std::priority_queue<ClusterContributionPair, std::vector<ClusterContributionPair>, 
        decltype(cluster_comparator)> split_queue(cluster_comparator);

    std::vector<float> norm_cache(cluster_contributions.cols());
    for(std::uint32_t i = 0; i < cluster_contributions.cols(); ++i){
        norm_cache[i] = 0.f;
        for(std::uint32_t j = 0; j < nn.size(); ++j){
            norm_cache[i] += cluster_contributions(nn[j] * 3, i) * cluster_contributions(nn[j] * 3, i);
            norm_cache[i] += cluster_contributions(nn[j] * 3 + 1, i) * cluster_contributions(nn[j] * 3 + 1, i);
            norm_cache[i] += cluster_contributions(nn[j] * 3 + 2, i) * cluster_contributions(nn[j] * 3 + 2, i);
        }
        norm_cache[i] = sqrt(norm_cache[i]);
        
    }

    for(std::uint32_t i = 0; i < splitting_clusters.size(); ++i){
        float cluster_contrib = computeCost(cluster_contributions, nn, splitting_clusters[i], norm_cache);
        split_queue.push(std::make_pair(i, cluster_contrib));
    }

    std::uint32_t total_clusters = splitting_clusters.size() + num_clusters;
    while(splitting_clusters.size() < total_clusters){
        auto largest = split_queue.top();
        
        split_queue.pop();

        std::uint32_t cluster_idx = largest.first;
        std::vector<std::uint32_t> cluster_cols = splitting_clusters[cluster_idx];
        std::vector<float> nz_col_norms;
        std::vector<std::uint32_t> nz_col_indices;

        for(std::uint32_t i = 0; i < cluster_cols.size(); ++i){
            if(norm_cache[cluster_cols[i]] > 0.){
                nz_col_norms.push_back(norm_cache[cluster_cols[i]]);
                nz_col_indices.push_back(cluster_cols[i]);
            }
        }

        std::set<std::uint32_t> lidx_set;

        if(nz_col_indices.size() > 1){
            std::discrete_distribution<std::uint32_t> nzdist(nz_col_norms.begin(), nz_col_norms.end());
            std::uint32_t tries = 0;
            do{
                if(tries++ > 100){
                    std::uint32_t index = 0;
                    while(lidx_set.size() < 2 && index < cluster_cols.size()){
                        lidx_set.insert(cluster_cols[index++]);
                    }

                    break;
                }

                lidx_set.insert(nz_col_indices[nzdist(rng)]);
            }while(lidx_set.size() < 2);
        }
        else{
            lidx_set.insert(cluster_cols[0]);
            lidx_set.insert(cluster_cols.back());
        }
        
        std::vector<std::uint32_t> lidx(lidx_set.begin(), lidx_set.end());

        Eigen::VectorXf line(nn.size() * 3);
        
        for(std::uint32_t i = 0; i < nn.size(); ++i){
            line(i * 3) = cluster_contributions(nn[i] * 3, lidx[1]) - cluster_contributions(nn[i] * 3, lidx[0]);
            line(i * 3 + 1) = cluster_contributions(nn[i] * 3 + 1, lidx[1]) - cluster_contributions(nn[i] * 3 + 1, lidx[0]);
            line(i * 3 + 2) = cluster_contributions(nn[i] * 3 + 2, lidx[1]) - cluster_contributions(nn[i] * 3 + 2, lidx[0]);
        }

        std::vector<std::pair<std::uint32_t, float>> vpl_projection_distances;

        float min_d = std::numeric_limits<float>::max();
        float max_d = -std::numeric_limits<float>::max();

        for(std::uint32_t i = 0; i < cluster_cols.size(); ++i){
            Eigen::VectorXf vpl_point(nn.size() * 3);
            for(std::uint32_t j = 0; j < nn.size(); ++j){
                vpl_point(j * 3) = cluster_contributions(nn[j] * 3, cluster_cols[i]) - cluster_contributions(nn[j] * 3, lidx[0]);
                vpl_point(j * 3 + 1) = cluster_contributions(nn[j] * 3 + 1, cluster_cols[i]) - cluster_contributions(nn[j] * 3 + 1, lidx[0]);
                vpl_point(j * 3 + 2) = cluster_contributions(nn[j] * 3 + 2, cluster_cols[i]) - cluster_contributions(nn[j] * 3 + 2, lidx[0]);
            }

            float d_on_line = vpl_point.dot(line);// / line.dot(line);
            min_d = std::min(d_on_line, min_d);
            max_d = std::max(d_on_line, max_d);
            vpl_projection_distances.push_back(std::make_pair(cluster_cols[i], d_on_line));
        }

        std::vector<std::uint32_t> new_cluster1;
        std::vector<std::uint32_t> new_cluster2;

        std::sort(vpl_projection_distances.begin(), vpl_projection_distances.end(),
            [](const std::pair<std::uint32_t, float>& lhs, const std::pair<std::uint32_t, float>& rhs){
                return lhs.second < rhs.second;
            });
        
        for(std::uint32_t i = 0; i < vpl_projection_distances.size(); ++i){
            if(i < (vpl_projection_distances.size() / 2)){
                new_cluster1.push_back(vpl_projection_distances[i].first);
            }
            else{
                new_cluster2.push_back(vpl_projection_distances[i].first);
            }
        }

        //in case all the points are the same, which is possible in the case of full occlusion
        if(new_cluster1.size() == 0 || new_cluster2.size() == 0){
            std::uint32_t split_idx = std::max(new_cluster1.size(), new_cluster2.size()) / 2;

            if(new_cluster1.size() > new_cluster2.size()){
                new_cluster2.insert(new_cluster2.begin(), new_cluster1.begin() + split_idx, new_cluster1.end());
            }
            else{
                new_cluster1.insert(new_cluster1.begin(), new_cluster2.begin() + split_idx, new_cluster2.end());
            }
        }

        double new_cluster1_contrib = computeCost(cluster_contributions, nn, new_cluster1, norm_cache);
        double new_cluster2_contrib = computeCost(cluster_contributions, nn, new_cluster2, norm_cache);

        splitting_clusters[largest.first] = new_cluster1;
        splitting_clusters.push_back(new_cluster2);

        split_queue.push(std::make_pair(largest.first, new_cluster1_contrib));
        split_queue.push(std::make_pair(splitting_clusters.size() - 1, new_cluster2_contrib));
    }

    return splitting_clusters;
}

//representative sampling for lightslice
std::vector<VPL> sampleRepresentatives(const Eigen::MatrixXf& contributions, const std::vector<VPL>& vpls, 
    const std::vector<std::vector<std::uint32_t>>& clusters, std::mt19937& rng, float min_dist){
    
    std::vector<VPL> representatives;

    //might need to make this more deterministic
    for(std::uint32_t i = 0; i < clusters.size(); ++i){
        if(clusters[i].size() == 0){
            continue;
        }
        
        std::vector<float> cluster_dist(clusters[i].size());
        float total_vpl_power = 0.f;
        for(std::uint32_t j = 0; j < clusters[i].size(); ++j){
            cluster_dist[j] = contributions.col(clusters[i][j]).norm();
            total_vpl_power += vpls[clusters[i][j]].P.getLuminance();
        }

        if(total_vpl_power > 0.f){
            std::discrete_distribution<std::uint32_t> gen(cluster_dist.begin(), cluster_dist.end());

            std::uint32_t rep_idx = gen(rng);
            float rep_power = vpls[clusters[i][rep_idx]].P.getLuminance();

            representatives.push_back(vpls[clusters[i][rep_idx]]);
            representatives.back().P = representatives.back().P * total_vpl_power / rep_power;
        }
        else{
            representatives.push_back(vpls[clusters[i][0]]);
            representatives.back().P = Spectrum(0.f);
        }
    }

    updateVPLRadii(representatives, min_dist);

    return representatives;
}

std::unique_ptr<KDTNode<ReconstructionSample>> constructKDTree(Scene* scene, std::uint32_t size_threshold, 
    std::vector<ReconstructionSample>& samples, float min_dist, std::uint32_t min_slice_size, std::uint32_t spp){

    auto kdt_root = std::unique_ptr<KDTNode<ReconstructionSample>>(new KDTNode<ReconstructionSample>(&samples));

    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    samples.resize(film->getSize().y * film->getSize().x * spp);

    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for (std::int32_t y = 0; y < film->getSize().y; ++y) {
        //#pragma omp parallel for
        for (std::int32_t x = 0; x < film->getSize().x; ++x) {
            std::uint32_t cell_dim = sqrt(spp) + 0.5f;
            float cell_side_len = 1.f / cell_dim;

            for(std::uint32_t i = 0; i < spp; ++i){
                //calculates the position the ray intersects with
                Ray ray;

                float x_jitter = 0.5f * cell_side_len;
                float y_jitter = 0.5f * cell_side_len;
                float x_off = (i % cell_dim) * cell_side_len;
                float y_off = (i / cell_dim) * cell_side_len;

                Point2 sample_position(x + x_off + x_jitter, y + y_off + y_jitter);

                Point2 aperture_sample(0.5f, 0.5f);
                Float time_sample(0.5f);

                sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

                ReconstructionSample curr_sample;
                curr_sample.image_x = x;
                curr_sample.image_y = y;
                curr_sample.ray = ray;
                curr_sample.color = Spectrum(0.f);

                std::uint32_t num_bounces = 0;

                while(true){
                    curr_sample.intersected_scene = scene->rayIntersect(curr_sample.ray, curr_sample.its);
                    if(!curr_sample.intersected_scene){
                        break;
                    }

                    if((curr_sample.its.getBSDF()->getType() & BSDF::ESmooth) || curr_sample.its.isEmitter()
                        || ++num_bounces > 10){
                        //curr_sample.its.wi = -curr_sample.ray.d;
                        break;
                    }

                    BSDFSamplingRecord bsdf_sample_record(curr_sample.its, sampler);
                    bsdf_sample_record.typeMask = curr_sample.its.getBSDF()->isDielectric() ? 
                        BSDF::EDeltaTransmission | BSDF::ENull : BSDF::EDeltaReflection;
                    
                    curr_sample.its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());

                    curr_sample.ray = Ray(curr_sample.its.p, bsdf_sample_record.its.toWorld(bsdf_sample_record.wo), ray.time);
                }

                if(!curr_sample.intersected_scene){
                    if(scene->hasEnvironmentEmitter()){
                        curr_sample.color = scene->evalEnvironment(RayDifferential(curr_sample.ray));
                    }
                    else curr_sample.color = Spectrum(0.f);
                }
                else if(curr_sample.its.isEmitter()){
                    curr_sample.color = curr_sample.its.Le(-curr_sample.ray.d);
                }

                samples[y * film->getSize().x * spp + x * spp + i] = std::move(curr_sample);
            }
        }
    }

    //only samples that intersected the scene are recovered. The rest are sampled from the environment map if there is one
    kdt_root->sample_indices.reserve(samples.size());
    for(std::uint32_t i = 0; i < samples.size(); ++i){
        if(samples[i].intersected_scene){
            kdt_root->sample_indices.push_back(i);
        }
    }

    splitKDTree(kdt_root.get(), size_threshold, min_slice_size, min_dist * 0.05f);

    return kdt_root;
}

std::uint32_t computeUnoccludedSamples(KDTNode<ReconstructionSample>* slice, bool vsl, Scene* scene, 
    const std::vector<VPL>& vpls, float min_dist){
    
    std::uint32_t total_samples = 0;
    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));
    
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        auto& curr_sample = slice->sample(i);
        curr_sample.unoccluded_samples.resize(vpls.size());

        for(std::uint32_t j = 0; j < vpls.size(); ++j){
            const VPL& vpl = vpls[j];
            std::uint32_t samples_for_cell;
            curr_sample.unoccluded_samples[j] = sample(scene, sampler, curr_sample.its, curr_sample.ray,
                vpl, min_dist, false, 5, false, curr_sample.intersected_scene, false, vsl, samples_for_cell);
                total_samples += samples_for_cell;
        }
    }

    return total_samples;
}

void updateSliceWithMatData(const Eigen::MatrixXf& mat, KDTNode<ReconstructionSample>* slice, 
    Scene* scene, const std::vector<VPL>& vpls, const std::vector<std::uint32_t>& actual_vpl_indices, float min_dist){

    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    std::uint32_t total_samples = slice->sample_indices.size();
    for(std::uint32_t i = 0; i < total_samples; ++i){
        if(slice->sample(i).its.isEmitter()){
            Spectrum emitter_col = slice->sample(i).its.Le(-slice->sample(i).ray.d);
            slice->sample(i).color = emitter_col;
        }
        else{
            for(std::uint32_t j = 0; j < actual_vpl_indices.size(); ++j){
                float coeff = (mat(i, j) + 1.f) / 2.f;
                std::uint32_t idx = actual_vpl_indices[j];
                slice->sample(i).color += slice->sample(i).unoccluded_samples[idx] * coeff;
            }
        }
    }
}

void copySamplesToBuffer(std::uint8_t* output_image, const std::vector<ReconstructionSample>& samples, Vector2i image_size,
    std::uint32_t spp){
    std::unordered_map<std::uint32_t, Spectrum> output;

    for(std::uint32_t i = 0; i < samples.size(); ++i){
        std::uint32_t buffer_pos = samples[i].image_x + samples[i].image_y * image_size.x;
        if(output.find(buffer_pos) != output.end()){
            output[buffer_pos] += samples[i].color;
        }
        else output[buffer_pos] = samples[i].color;
    }

    for(auto iter = output.begin(); iter != output.end(); ++iter){
        Spectrum col = iter->second;
        col /= spp;

        float r, g, b;
        col.toSRGB(r, g, b);

        r = std::max(0.f, std::min(1.f, r));
        g = std::max(0.f, std::min(1.f, g));
        b = std::max(0.f, std::min(1.f, b));

        std::uint32_t buffer_pos = iter->first;

        output_image[buffer_pos*3] = r * 255 + 0.5f;
        output_image[buffer_pos*3 + 1] = g * 255 + 0.5f;
        output_image[buffer_pos*3 + 2] = b * 255 + 0.5f;
    }
}

//generates a set of indices based on the probabilities passed through the vector
std::vector<std::uint32_t> importanceSample(std::uint32_t num_samples, std::mt19937& rng, const std::vector<float>& probabilities){
    num_samples = std::min((size_t)num_samples, probabilities.size());

    if(num_samples == probabilities.size()){
        std::vector<std::uint32_t> all_indices(probabilities.size());
        std::iota(all_indices.begin(), all_indices.end(), 0);
        return all_indices;
    }

    std::discrete_distribution<std::uint32_t> gen(probabilities.begin(), probabilities.end());

    std::vector<std::uint32_t> sampled(num_samples);
    for(std::uint32_t i = 0; i < num_samples; ++i){
        sampled[i] = gen(rng);
    }

    return sampled;
}

//sparsely samples a single column for adaptive matrix recovery
std::vector<std::uint32_t> sampleCol(Scene* scene, KDTNode<ReconstructionSample>* slice, const std::vector<VPL>& vpls, 
    std::uint32_t col, float min_dist, std::uint32_t num_samples, std::mt19937& rng, Eigen::MatrixXf& mat, 
    const std::vector<std::uint32_t>& sample_set, bool resample,
    bool importance_sample, const std::vector<float>& probabilities){
    
    std::uint32_t num_rows = num_samples;
    std::uint32_t max_samples = slice->sample_indices.size();

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

    for(size_t i = 0; i < sampled_indices.size(); ++i){
        std::uint32_t vpl_index = col;
        std::uint32_t sample_index = sampled_indices[i];
        const VPL& vpl = vpls[vpl_index];
        ReconstructionSample& scene_sample = slice->sample(sample_index);

        bool los = sampleVisibility(scene, scene_sample.its, vpl, min_dist);
        mat(i, 0) = los ? 1.f : -1.f;
    }

    return sampled_indices;
}

//sample col as boolean
std::vector<std::uint32_t> sampleColB(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, std::uint32_t col, float min_dist, std::uint32_t num_samples, 
    std::mt19937& rng, std::unordered_map<std::uint32_t, std::uint8_t>& sampled_vals, const std::vector<float>& probabilities, 
    const std::vector<std::uint32_t>& sample_set, bool resample){

    std::uint32_t max_samples = slice->sample_indices.size();
    assert(num_samples <= max_samples);

    std::vector<std::uint32_t> sampled_indices = resample ? importanceSample(num_samples, rng, probabilities) : sample_set;

    for(size_t i = 0; i < sampled_indices.size(); ++i){
        const VPL& vpl = vpls[col];
        ReconstructionSample& scene_sample = slice->sample(sampled_indices[i]);

        if(sampled_vals.find(sampled_indices[i]) == sampled_vals.end()){
            sampled_vals[sampled_indices[i]] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;    
        }
    }

    return sampled_indices;
}

//sample col with leading indices for gaussian elimination
std::vector<std::uint32_t> sampleColBWithLeading(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, std::uint32_t col, float min_dist, std::uint32_t num_samples, 
    std::mt19937& rng, std::unordered_map<std::uint32_t, std::uint8_t>& sampled_vals, const std::vector<float>& leading_probabilities, 
    const std::vector<float>& non_leading_probabilities, const std::vector<std::uint32_t>& leading_indices, 
    float leading_perc, const std::vector<std::uint32_t>& sample_set, bool resample){

    std::uint32_t max_samples = slice->sample_indices.size();
    assert(num_samples <= max_samples);

    std::vector<std::uint32_t> sample_indices;

    if(resample){
        std::uint32_t num_leading_samples = std::min(std::uint32_t(leading_indices.size()), std::max(1u, std::uint32_t(leading_perc * num_samples)));
        if(num_leading_samples == leading_indices.size()){
            for(std::uint32_t i = 0; i < leading_indices.size(); ++i){
                sample_indices.push_back(leading_indices[i]);
            }
        }
        else{
            sample_indices = importanceSample(num_leading_samples, rng, leading_probabilities);
        }

        std::uint32_t remaining_samples = num_samples - num_leading_samples;
        std::vector<std::uint32_t> nonleading_sample_indices;
        if(remaining_samples > 0){
            nonleading_sample_indices = importanceSample(remaining_samples, rng, non_leading_probabilities);
        }
        sample_indices.insert(sample_indices.end(), nonleading_sample_indices.begin(), nonleading_sample_indices.end());
    }
    else{
        sample_indices = sample_set;
    }

    for(size_t i = 0; i < sample_indices.size(); ++i){
        const VPL& vpl = vpls[col];
        ReconstructionSample& scene_sample = slice->sample(sample_indices[i]);

        if(sampled_vals.find(sample_indices[i]) == sampled_vals.end()){
            sampled_vals[sample_indices[i]] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;    
        }
    }

    return sample_indices;
}

//gaussian elimination for matrix over gf2
std::vector<std::vector<std::uint8_t>> gf2elim(const std::vector<std::vector<std::uint8_t>>& basis, std::vector<std::uint32_t>& leading_pos){
    std::vector<std::vector<std::uint8_t>> reduced_basis = basis;
    
    if(basis.size() == 0){
        return reduced_basis;
    }

    //swap rows and columns because we perform elimination on a transposed matrix
    std::uint32_t rows = reduced_basis.size();
    std::uint32_t cols = reduced_basis[0].size();

    //actual gaussian elimination algorithm for gf2, should be noted that we are reducing the transpose of the above matrix, so all rows and columns 
    //are inherently changed.
    std::uint32_t curr_pivot_row = 0;
    std::uint32_t curr_pivot_col = 0;

    while(curr_pivot_row < rows && curr_pivot_col < cols){
        //this is to handle empty columns. we forward the column pivot to the next non-empty column of the submatrix
        for(; curr_pivot_col < cols; ++curr_pivot_col){
            std::uint8_t has_one = 0;

            for(std::uint32_t i = curr_pivot_row; i < rows; ++i){
                has_one |= reduced_basis[i][curr_pivot_col];
            }

            if(has_one){
                break;
            }
        }

        //finished
        if(curr_pivot_col >= cols){
            break;
        }

        //put first row with nonzero pivot element as the first
        int first_nonzero_idx = curr_pivot_row;
        for(std::uint32_t i = curr_pivot_row; i < reduced_basis.size(); ++i){
            if(reduced_basis[i][curr_pivot_col]){
                first_nonzero_idx = i;
                break;
            }
        }

        if(first_nonzero_idx != int(curr_pivot_row)){
            auto temp = reduced_basis[curr_pivot_row];
            reduced_basis[curr_pivot_row] = reduced_basis[first_nonzero_idx];
            reduced_basis[first_nonzero_idx] = temp;
        }

        //xor pivot row with all other rows in the matrix
        std::vector<std::uint8_t> aijn(cols - curr_pivot_col);
        for(std::uint32_t i = curr_pivot_col; i < cols; ++i){
            aijn[i - curr_pivot_col] = reduced_basis[curr_pivot_row][i];
        }

        std::vector<std::uint8_t> c(rows);
        for(std::uint32_t i = 0; i < rows; ++i){
            c[i] = reduced_basis[i][curr_pivot_col];
        }

        c[curr_pivot_row] = 0; //dont self xor pivot row

        for(std::uint32_t i = 0; i < c.size(); ++i){
            if(c[i]){
                for(std::uint32_t j = 0; j < aijn.size(); ++j){
                    reduced_basis[i][j + curr_pivot_col] ^= aijn[j];
                }
            }
        }

        curr_pivot_row++;
        curr_pivot_col++;
    }

    //drop all zero rows
    reduced_basis.resize(curr_pivot_row);

    //get all leading positions, will be needed when obtaining coefficients
    leading_pos.clear();

    for(std::uint32_t i = 0; i < reduced_basis.size(); ++i){
        int nonzero_pos = -1;
        for(std::uint32_t j = 0; j < reduced_basis[i].size(); ++j){
            if(reduced_basis[i][j]){
                nonzero_pos = j;
                break;
            }
        }

        assert(nonzero_pos >= 0);
        leading_pos.push_back(nonzero_pos);
    }

    return reduced_basis;
}

//basis is gaussian eliminated version of actual basis
bool gereconstruct(std::unordered_map<std::uint32_t, std::uint8_t>& sampled, const std::vector<std::vector<std::uint8_t>>& reduced_basis, 
    const std::vector<std::uint32_t>& leading_indices, std::uint32_t rows){

    std::vector<std::uint32_t> one_counts(rows, 0);

    for(std::uint32_t i = 0; i < leading_indices.size(); ++i){
        //consider
        std::uint32_t idx = leading_indices[i];

        //only consider basis if it's pivot has been sampled
        if(sampled.find(idx) != sampled.end()){
            bool even = (one_counts[idx] & 1) == 0;

            //check if need to consider basis, if yes update tally, this works because the matrix is at least in echelon form after
            //gaussian elimination, thus we know the column will be zero from now onwards
            if((sampled[idx] == 1 && even) || (sampled[idx] == 0 && !even)){
                for(std::uint32_t j = 0; j < one_counts.size(); ++j){
                    one_counts[j] += reduced_basis[i][j];
                }
            }
        }
    }

    for(std::uint32_t i = 0; i < one_counts.size(); ++i){
        one_counts[i] &= 1;
    }

    bool matching = true;
    for(auto iter = sampled.begin(); iter != sampled.end(); ++iter){
        std::uint32_t idx = iter->first;
        if(iter->second != one_counts[idx]){
            matching = false;
            break;
        }
    }

    if(matching){
        for(std::uint32_t i = 0; i < rows; ++i){
            sampled[i] = one_counts[i];
        }
    }


    return matching;
}

//adaptive matrix completion with gaussian elimination over gf2
std::uint32_t adaptiveMatrixReconstructionBGE(
    std::vector<std::uint8_t>& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc, float max_sample_perc, float verification_inc, 
    std::mt19937& rng, const std::vector<float>& col_estimations){

    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);

    std::uint32_t num_rows = slice->sample_indices.size();
    std::uint32_t num_cols = vpls.size();
    
    //re-allocate matrix if it is of the incorrect size
    if(mat.size() != slice->sample_indices.size() * vpls.size()){
        mat.resize(slice->sample_indices.size() * vpls.size(), 0);
    }

    std::uint32_t num_samples = num_rows * sample_perc + 0.5f;
    std::uint32_t num_verification_samples = verification_inc < std::numeric_limits<float>::epsilon() ? 
        0 : std::max(1, int(num_rows * verification_inc + 0.5f));
    //just in case, this shouldn't ever really happen
    if(num_samples == 0){
        num_samples = slice->sample_indices.size();
    }

    //we recover from the brightest to least bright because there will be more full samples initially, allow for better coverage of
    //higher energy sections
    std::vector<std::uint32_t> order(num_cols);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), 
        [&col_estimations](const std::uint32_t& lhs, const std::uint32_t& rhs){
            return col_estimations[lhs] > col_estimations[rhs];
        });


    std::vector<std::vector<std::uint8_t>> basis;
    std::vector<std::vector<std::uint8_t>> reduced_basis;
    std::vector<std::uint32_t> leading_indices;

    std::uint32_t total_samples = 0;
    std::uniform_real_distribution<float> gen(0.f, 1.f);
    std::vector<float> probabilities(num_rows, 1.f / num_rows);

    std::vector<float> leading_probabilities(num_rows, 0.f);
    std::vector<float> non_leading_probabilities(num_rows, 0.f);

    std::unordered_map<std::uint32_t, std::uint8_t> sample_omega;
    std::vector<std::uint8_t> col_to_add(num_rows);
    bool full_col_sampled = false;

    std::vector<std::uint32_t> sampled;

    //The actual adaptive matrix recovery algorithm
    for(std::uint32_t i = 0; i < order.size(); ++i){
        sample_omega.clear();
        std::uint32_t samples_for_col = 0;

       if(basis.size() > 0){
            sampled = sampleColBWithLeading(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega,
                leading_probabilities, non_leading_probabilities, leading_indices, 0.5f, sampled, true);

            full_col_sampled = false;

            samples_for_col = num_samples;

            if(sample_omega.size() == col_to_add.size()){
                for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                    col_to_add[j] = sample_omega[j];   
                }
            }
            else{
                std::unordered_map<std::uint32_t, std::uint8_t> reconstructed = sample_omega; 
                if(gereconstruct(reconstructed, reduced_basis, leading_indices, num_rows)){
                    for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                        col_to_add[j] = reconstructed[j];   
                    }

                    if(num_verification_samples > 0){
                        std::vector<std::uint32_t> ver_indices;

                        std::uniform_int_distribution<std::uint32_t> gen_row(0, num_rows - 1);
                        while(ver_indices.size() < num_verification_samples){
                            std::uint32_t row = gen_row(rng);
                            if(sample_omega.find(row) == sample_omega.end()){
                                ver_indices.push_back(row);
                            }
                        }


                        sampleColB(scene, slice, vpls, order[i], min_dist, num_verification_samples, rng, sample_omega, 
                            probabilities, ver_indices, false);

                        bool correct = true;
                        for(std::uint32_t j = 0; j < ver_indices.size(); ++j){
                            if(sample_omega[ver_indices[j]] != col_to_add[ver_indices[j]]){
                                correct = false;
                                break;
                            }
                        }

                        if(!correct){
                            sample_perc = std::min(max_sample_perc, sample_perc + verification_inc);
                            num_samples = num_rows * sample_perc + 0.5f;

                            sampleColB(scene, slice, vpls, order[i], min_dist, num_rows, rng, sample_omega, 
                                probabilities, sampled, true);
                            samples_for_col = num_rows;

                            for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                                col_to_add[j] = sample_omega[j];   
                            }
                            
                            full_col_sampled = true;
                        }
                    }
                }
                else{
                    sampleColB(scene, slice, vpls, order[i], min_dist, num_rows, rng, sample_omega, 
                        probabilities, sampled, true);
                    samples_for_col = num_rows;

                    for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                        col_to_add[j] = sample_omega[j];   
                    }
                    
                    full_col_sampled = true;
                }
            }
        }
        else{
            sampled = sampleColB(scene, slice, vpls, order[i], min_dist, num_rows, rng, sample_omega, 
                probabilities, sampled, true);
            samples_for_col = num_rows;

            for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                col_to_add[j] = sample_omega[j];   
            }
            
            full_col_sampled = true;
        }

        //probability and basis updates
        if(full_col_sampled){
            //basis update
            basis.push_back(col_to_add);
            reduced_basis = gf2elim(basis, leading_indices);

            //probability update
            std::vector<std::uint32_t> buckets;
            std::uint8_t last_sign = 2;
            for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                if(last_sign != col_to_add[j]){
                    buckets.push_back(j);
                    last_sign = col_to_add[j];
                }
            }

            buckets.push_back(col_to_add.size());
            float prob_per_bucket = 1.f / buckets.size();

            for(std::uint32_t j = 0; j < buckets.size(); ++j){
                std::uint32_t bucket_size = j == 0 ? buckets[0] : buckets[j] - buckets[j - 1];
                float curr_bucket_prob = prob_per_bucket / bucket_size;
                std::uint32_t bucket_start = j == 0 ? 0 : buckets[j - 1];
                std::uint32_t bucket_end = buckets[j];

                for(std::uint32_t k = bucket_start; k < bucket_end; ++k){
                    probabilities[k] += curr_bucket_prob;
                }
            }

            //set leading and non-leading probabilities
            non_leading_probabilities = probabilities;
            for(std::uint32_t j = 0; j < leading_indices.size(); ++j){
                std::uint32_t idx = leading_indices[j];
                non_leading_probabilities[idx] = 0.f;
                leading_probabilities[idx] = 1.f;
            }
        }

        std::uint32_t offset = order[i] * num_rows;
        std::copy(col_to_add.begin(), col_to_add.end(), mat.begin() + offset);
        total_samples += samples_for_col;
    }

    return total_samples;
}

//original adaptive matrix completion algorithm 
std::uint32_t adaptiveMatrixReconstruction(Eigen::MatrixXf& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc, std::mt19937& rng, bool importance_sample,
    const std::vector<float>& col_estimations){

    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);

    std::uint32_t total_rows = slice->sample_indices.size();
    std::uint32_t num_rows = total_rows;
    std::uint32_t num_cols = vpls.size();
    
    //re-allocate matrix if it is of the incorrect size
    if((size_t)mat.cols() != num_cols || (size_t)mat.rows() != num_rows){
        mat = Eigen::MatrixXf(num_rows, num_cols);
    }
    mat.setZero();

    std::uint32_t num_samples = total_rows * sample_perc + 0.5f;
    //just in case, this shouldn't ever really happen
    if(num_samples == 0){
        num_samples = total_rows;
    }

    //we recover from the brightest to least bright because there will be more full samples initially, allow for better coverage of
    //higher energy sections
    std::vector<std::uint32_t> order(num_cols);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), 
        [&col_estimations](const std::uint32_t& lhs, const std::uint32_t& rhs){
            return col_estimations[lhs] > col_estimations[rhs];
        });


    Eigen::MatrixXf reconstructed(num_rows, 1);
    std::vector<std::uint32_t> sampled;

    Eigen::MatrixXf q;
    Eigen::MatrixXf sample_omega;
    Eigen::MatrixXf q_omega;
    Eigen::MatrixXf q_omega_pseudoinverse;
    Eigen::MatrixXf qq_omega_pseudoinverse;

    std::uint32_t total_samples = 0;
    std::uniform_real_distribution<float> gen(0.f, 1.f);
    std::vector<float> probabilities(num_rows, 1.f / num_rows);

    bool full_col_sampled = false;

    //The actual adaptive matrix recovery algorithm
    for(std::uint32_t i = 0; i < order.size(); ++i){
        
        std::uint32_t samples_for_col = 0;

        //if the basis is not empty, we can try reproject, otherwise a full sample is required to populate the basis
        if(q.cols() > 0){
            std::uint32_t expected_omega_rows = num_samples;
            if(sample_omega.cols() != expected_omega_rows){
                sample_omega.resize(expected_omega_rows, 1);
            }

            //we may want to regenerate the sample indices for a variety of reasons, in which case the indices are generated and
            //the pseudoinverse is recalculated
            if(sample_omega.rows() != reconstructed.rows() && full_col_sampled){
                full_col_sampled = false;

                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    true, importance_sample, probabilities);

                q_omega.resize(expected_omega_rows, q.cols());

                for(std::uint32_t j = 0; j < sampled.size(); ++j){
                    q_omega.row(j) = q.row(sampled[j]);
                }

                auto svd = q_omega.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto sv = svd.singularValues();
                Eigen::MatrixXf singular_val_inv = Eigen::MatrixXf::Zero(sv.size(), sv.size());

                for(std::uint32_t j = 0; j < sv.size(); ++j){
                    singular_val_inv(j, j) = sv(j) < 1e-5f ? 0.f : 1.f / sv(j);
                }

                q_omega_pseudoinverse = svd.matrixV() * singular_val_inv * svd.matrixU().transpose();
                qq_omega_pseudoinverse = q * q_omega_pseudoinverse;
            }
                //no new direction was added so no need to regenerate sample indices
            else{
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    false, importance_sample, probabilities);
            }
            samples_for_col = num_samples;

            //no need to reconstruct if full sample
            reconstructed = sample_omega.rows() == reconstructed.rows() ? sample_omega : qq_omega_pseudoinverse * sample_omega;           
            
            float d = 0;
            for(std::uint32_t j = 0; j < sampled.size(); ++j){
                d += std::abs(reconstructed(sampled[j], 0) - sample_omega(j, 0));
            }

            //sampled values can't be reconstructed accurately so fully sample
            if(d > 1e-3f){
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, 
                    reconstructed, sampled, true, importance_sample, probabilities);
                samples_for_col = total_rows;

                q.conservativeResize(q.rows(), q.cols() + 1);
                q.col(q.cols() - 1) = reconstructed;
                
                full_col_sampled = true;
            }
        }
        else{
            sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, reconstructed, 
                sampled, true, importance_sample, probabilities);
            samples_for_col = total_rows;

            //only add to basis if vector isn't a zero vector, as it is otherwise meaningless
            //might want to change this to be above some epsilon instead
            if(reconstructed.norm() > 0.f){
                q = reconstructed;
                full_col_sampled = true;
            }
        }

        //probability update for importance sampling.
        if(full_col_sampled){
            std::vector<std::uint32_t> buckets;
            std::int32_t last_sign = 0;
            for(std::uint32_t j = 0; j < reconstructed.rows(); ++j){
                std::int32_t sign = reconstructed(j, 0) > 0.f ? 1 : -1;
                if(sign != last_sign){
                    buckets.push_back(j);
                    last_sign = sign;
                }
            }

            buckets.push_back(reconstructed.rows());
            float prob_per_bucket = 1.f / buckets.size();

            for(std::uint32_t j = 0; j < buckets.size(); ++j){
                std::uint32_t bucket_size = j == 0 ? buckets[0] : buckets[j] - buckets[j - 1];
                float curr_bucket_prob = prob_per_bucket / bucket_size;
                std::uint32_t bucket_start = j == 0 ? 0 : buckets[j - 1];
                std::uint32_t bucket_end = buckets[j];

                for(std::uint32_t k = bucket_start; k < bucket_end; ++k){
                    probabilities[k] += curr_bucket_prob;
                }
            }
        }
        mat.col(order[i]) = reconstructed.col(0);
        total_samples += samples_for_col;
    }

    return total_samples;
}

std::vector<std::pair<std::uint32_t, bool>> getMatchingCols(const std::vector<std::vector<std::uint8_t>>& cols, 
    const std::unordered_map<std::uint32_t, std::uint8_t>& sampled_vals){
    
    std::vector<std::pair<std::uint32_t, bool>> matching_cols;
    for(std::uint32_t i = 0; i < cols.size(); ++i){
        bool matching = true;
        bool opposite = true;

        for(auto iter = sampled_vals.begin(); iter != sampled_vals.end(); ++iter){
            if(iter->second == cols[i][iter->first]){
                opposite = false;
            }
            else{
                matching = false;
            }
        }

        if(matching || opposite){
            matching_cols.push_back(std::make_pair(i, matching));
        }
    }

    return matching_cols;
}


//adaptive matrix completion with binary matching
std::uint32_t adaptiveMatrixReconstructionB(
    std::vector<std::uint8_t>& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc, float max_sample_perc, float verification_inc, 
    std::mt19937& rng, const std::vector<float>& col_estimations){

    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);

    std::uint32_t num_rows = slice->sample_indices.size();
    std::uint32_t num_cols = vpls.size();
    
    //re-allocate matrix if it is of the incorrect size
    if(mat.size() != slice->sample_indices.size() * vpls.size()){
        mat.resize(slice->sample_indices.size() * vpls.size(), 0);
    }

    std::uint32_t num_samples = num_rows * sample_perc + 0.5f;
    std::uint32_t num_verification_samples = verification_inc < std::numeric_limits<float>::epsilon() ? 
        0 : std::max(1, int(num_rows * verification_inc + 0.5f));
    //just in case, this shouldn't ever really happen
    if(num_samples == 0){
        num_samples = slice->sample_indices.size();
    }

    //we recover from the brightest to least bright because there will be more full samples initially, allow for better coverage of
    //higher energy sections
    std::vector<std::uint32_t> order(num_cols);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), 
        [&col_estimations](const std::uint32_t& lhs, const std::uint32_t& rhs){
            return col_estimations[lhs] > col_estimations[rhs];
        });


    std::vector<std::vector<std::uint8_t>> basis;

    std::uint32_t total_samples = 0;
    std::uniform_real_distribution<float> gen(0.f, 1.f);
    std::vector<float> probabilities(num_rows, 1.f / num_rows);

    std::unordered_map<std::uint32_t, std::uint8_t> sample_omega;
    std::vector<std::uint8_t> col_to_add(num_rows);
    bool full_col_sampled = false;

    std::vector<std::uint32_t> sampled;

    //The actual adaptive matrix recovery algorithm
    for(std::uint32_t i = 0; i < order.size(); ++i){
        sample_omega.clear();
        std::uint32_t samples_for_col = 0;

       if(basis.size() > 0){
            sampled = sampleColB(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, 
                probabilities, sampled, true);
            full_col_sampled = false;

            samples_for_col = num_samples;
            
            std::vector<std::pair<std::uint32_t, bool>> matching_cols;

            if(sample_omega.size() == col_to_add.size()){
                for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                    col_to_add[j] = sample_omega[j];   
                }
            }
            else{
                matching_cols = getMatchingCols(basis, sample_omega);

                if(matching_cols.size() > 0){
                    std::uniform_int_distribution<std::uint32_t> select_col(0, matching_cols.size() - 1);
                    std::uint32_t sel = select_col(rng);
                    col_to_add = basis[matching_cols[sel].first];

                    //opposite direction so we flip
                    if(!matching_cols[sel].second){
                        for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                            col_to_add[j] = (col_to_add[j] + 1) % 2;
                        }
                    }

                    if(num_verification_samples > 0){
                        std::vector<std::uint32_t> ver_indices;

                        std::uniform_int_distribution<std::uint32_t> gen_row(0, num_rows - 1);
                        while(ver_indices.size() < num_verification_samples){
                            std::uint32_t row = gen_row(rng);
                            if(sample_omega.find(row) == sample_omega.end()){
                                ver_indices.push_back(row);
                            }
                        }

                        sampleColB(scene, slice, vpls, order[i], min_dist, num_verification_samples, rng, sample_omega, 
                            probabilities, ver_indices, false);

                        bool correct = true;
                        for(std::uint32_t j = 0; j < ver_indices.size(); ++j){
                            if(sample_omega[ver_indices[j]] != col_to_add[ver_indices[j]]){
                                correct = false;
                                break;
                            }
                        }

                        if(!correct){
                            sample_perc = std::min(max_sample_perc, sample_perc + verification_inc);
                            num_samples = num_rows * sample_perc + 0.5f;

                            sampleColB(scene, slice, vpls, order[i], min_dist, num_rows, rng, sample_omega, 
                                probabilities, sampled, true);
                            samples_for_col = num_rows;

                            for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                                col_to_add[j] = sample_omega[j];   
                            }

                            basis.push_back(col_to_add);
                            
                            full_col_sampled = true;
                        }
                    }
                }
                else{
                    sampleColB(scene, slice, vpls, order[i], min_dist, num_rows, rng, sample_omega, 
                        probabilities, sampled, true);
                    samples_for_col = num_rows;

                    for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                        col_to_add[j] = sample_omega[j];   
                    }
    
                    basis.push_back(col_to_add);
                    
                    full_col_sampled = true;
                }
            }
        }
        else{
            sampled = sampleColB(scene, slice, vpls, order[i], min_dist, num_rows, rng, sample_omega, 
                probabilities, sampled, true);
            samples_for_col = num_rows;

            for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                col_to_add[j] = sample_omega[j];   
            }

            basis.push_back(col_to_add);
            
            full_col_sampled = true;
        }

        //probability update for importance sampling.
        if(full_col_sampled){
            std::vector<std::uint32_t> buckets;
            std::uint8_t last_sign = 2;
            for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                if(last_sign != col_to_add[j]){
                    buckets.push_back(j);
                    last_sign = col_to_add[j];
                }
            }

            buckets.push_back(col_to_add.size());
            float prob_per_bucket = 1.f / buckets.size();

            for(std::uint32_t j = 0; j < buckets.size(); ++j){
                std::uint32_t bucket_size = j == 0 ? buckets[0] : buckets[j] - buckets[j - 1];
                float curr_bucket_prob = prob_per_bucket / bucket_size;
                std::uint32_t bucket_start = j == 0 ? 0 : buckets[j - 1];
                std::uint32_t bucket_end = buckets[j];

                for(std::uint32_t k = bucket_start; k < bucket_end; ++k){
                    probabilities[k] += curr_bucket_prob;
                }
            }
        }

        std::uint32_t offset = order[i] * num_rows;
        std::copy(col_to_add.begin(), col_to_add.end(), mat.begin() + offset);
        total_samples += samples_for_col;
    }

    return total_samples;
}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(const std::vector<VPL>& vpls,
    float sample_percentage, float min_dist, std::uint32_t slice_size,  
    bool adaptive_importance_sampling, bool show_slices, bool vsl, ClusteringStrategy clustering_strategy, 
    bool hw, bool bin_vis, std::uint32_t num_clusters, std::uint32_t samples_per_slice, float max_sample_perc, float ver_inc, bool ge) : 
        vpls_(vpls), 
        sample_percentage_(sample_percentage), 
        min_dist_(min_dist), 
        slice_size_(slice_size), 
        adaptive_importance_sampling_(adaptive_importance_sampling),
        vsl_(vsl),
        clustering_strategy_(clustering_strategy),
        hw_(hw),
        bin_vis_(bin_vis),
        num_clusters_(num_clusters),
        samples_per_slice_(samples_per_slice),
        sample_inc_(ver_inc), 
        max_sample_perc_(max_sample_perc),
        ge_(ge),
        cancel_(false){
}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other) : 
    vpls_(other.vpls_),
    sample_percentage_(other.sample_percentage_), 
    min_dist_(other.min_dist_), 
    slice_size_(other.slice_size_), 
    adaptive_importance_sampling_(other.adaptive_importance_sampling_),
    vsl_(other.vsl_),
    clustering_strategy_(other.clustering_strategy_),
    hw_(other.hw_),
    bin_vis_(other.bin_vis_),
    num_clusters_(other.num_clusters_),
    samples_per_slice_(other.samples_per_slice_),
    sample_inc_(other.sample_inc_), 
    max_sample_perc_(other.max_sample_perc_),
    ge_(other.ge_),
    cancel_(other.cancel_){
}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (MatrixReconstructionRenderer&& other){
    if(this != &other){
        vpls_ = other.vpls_;
        sample_percentage_ = other.sample_percentage_;
        min_dist_ = other.min_dist_;
        slice_size_ = other.slice_size_;
        adaptive_importance_sampling_ = other.adaptive_importance_sampling_;
        vsl_ = other.vsl_;
        clustering_strategy_ = other.clustering_strategy_;
        hw_ = other.hw_;
        bin_vis_ = other.bin_vis_;
        num_clusters_ = other.num_clusters_;
        samples_per_slice_ = other.samples_per_slice_;
        sample_inc_ = other.sample_inc_; 
        max_sample_perc_ = other.max_sample_perc_;
        ge_ = other.ge_;
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

struct AdaptiveConstructionParams{
    bool import_sample;
    bool bin_vis;
};

struct GeneralParams{
    float min_dist;
    bool adaptive_col;
    bool vsl;
    float sample_perc;
    float max_sample_perc;
    float sample_inc;
};

std::uint32_t getQuadrant(Vector3f v){
    std::uint32_t quadrant = 0;

    if(v.x < 0.f){
        quadrant |= 1;
    }

    if(v.y < 0.f){
        quadrant |= 2;
    }

    if(v.z < 0.f){
        quadrant |= 4;
    }

    return quadrant;
}

void recover(KDTNode<ReconstructionSample>* slice, std::mutex& stats_mutex, const std::vector<VPL>& vpls, 
    Scene* scene, const GeneralParams& general_params, const AdaptiveConstructionParams& ac_params,
    std::uint64_t& total_samples, std::uint64_t& performed_samples, std::mt19937& rng, bool ge){

    computeUnoccludedSamples(slice, general_params.vsl, scene, vpls, general_params.min_dist);

    Point3f slice_center_of_mass(0.f, 0.f, 0.f);
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        slice_center_of_mass += slice->sample(i).its.p;
    }
    slice_center_of_mass /= slice->sample_indices.size();

    std::vector<std::vector<VPL>> quadranted_vpls(8);
    std::vector<std::vector<std::uint32_t>> quadranted_vpl_indices(8);

    //further slice for importance sampling purposes
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        Vector3f pos_dir = normalize(slice_center_of_mass - vpls[i].its.p);
        Vector3f dir = vpls[i].type == EDirectionalEmitterVPL ? Vector3f(vpls[i].its.shFrame.n) : pos_dir;

        std::uint32_t quadrant = getQuadrant(dir);
        quadranted_vpls[quadrant].push_back(vpls[i]);
        quadranted_vpl_indices[quadrant].push_back(i);
    }

    std::uint32_t total_performed_samples = 0;

    for(std::uint32_t i = 0; i < quadranted_vpls.size(); ++i){
        if(quadranted_vpls[i].size() == 0){
            continue;
        }

        std::vector<float> cluster_contribs(quadranted_vpls[i].size(), 0.f);

        for(std::uint32_t j = 0; j < slice->sample_indices.size(); ++j){
            for(std::uint32_t k = 0; k < quadranted_vpls[i].size(); ++k){
                cluster_contribs[k] += slice->sample(j).unoccluded_samples[quadranted_vpl_indices[i][k]].getLuminance();
            }
        }

        std::uint32_t mat_rows = slice->sample_indices.size();
        Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(mat_rows, quadranted_vpls[i].size());
        std::uint32_t samples;

        if(!ac_params.bin_vis){
            samples = adaptiveMatrixReconstruction(mat, scene, slice, quadranted_vpls[i], general_params.min_dist, 
                general_params.sample_perc, rng, ac_params.import_sample, cluster_contribs);
        }
        else{
            std::vector<std::uint8_t> bin_vis(slice->sample_indices.size() * quadranted_vpls[i].size(), 0);
            if(ge){
                samples = adaptiveMatrixReconstructionBGE(bin_vis, scene, slice,
                    quadranted_vpls[i], general_params.min_dist, general_params.sample_perc, 
                    general_params.max_sample_perc, general_params.sample_inc, rng, 
                    cluster_contribs);
            }
            else{
                samples = adaptiveMatrixReconstructionB(bin_vis, scene, slice,
                    quadranted_vpls[i], general_params.min_dist, general_params.sample_perc, 
                    general_params.max_sample_perc, general_params.sample_inc, rng, 
                    cluster_contribs);
            }
            

            for(auto col = 0; col < mat.cols(); ++col){
                for(auto row = 0; row < mat.rows(); ++row){
                    mat(row, col) = bin_vis[col * slice->sample_indices.size() + row] == 0 ? -1.f : 1.f;
                }
            }
        }

        total_performed_samples += samples;
        updateSliceWithMatData(mat, slice, scene, vpls, quadranted_vpl_indices[i], general_params.min_dist);
    }

    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        slice->sample(i).unoccluded_samples.clear();
        slice->sample(i).unoccluded_samples.shrink_to_fit();
    }

    {
        std::lock_guard<std::mutex> lock(stats_mutex);
        performed_samples += total_performed_samples;
        total_samples += slice->sample_indices.size() * vpls.size();
    }
}

void sliceWorkerLS(std::vector<std::int32_t>& work, std::uint32_t thread_id, std::mutex& work_mutex, std::mutex& stats_mutex,
    const std::vector<KDTNode<ReconstructionSample>*>& slices, const std::vector<std::vector<std::uint32_t>>& cbsamp, 
    const Eigen::MatrixXf& contribs, const std::vector<std::vector<int>>& nn, const std::vector<VPL>& total_vpls, Scene* scene,
    const GeneralParams& general_params, const AdaptiveConstructionParams& ac_params,
    std::uint64_t& total_samples, std::uint64_t& performed_samples, std::uint32_t num_clusters, bool ge){

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * thread_id);
    std::uint32_t split_num_clusters = std::min(total_vpls.size(), num_clusters - cbsamp.size());

    while(true){
        std::int32_t slice_id;
        {
            std::lock_guard<std::mutex> lock(work_mutex);
            slice_id = work[thread_id];
        }

        if(slice_id == -2){
            break;
        }
        else if(slice_id == -1){
            continue;
        }

        std::vector<std::vector<std::uint32_t>> cbsplit = clusterVPLsBySplitting(cbsamp, contribs, nn[slice_id], split_num_clusters, rng);
        auto vpls = sampleRepresentatives(contribs, total_vpls, cbsplit, rng, general_params.min_dist);

        recover(slices[slice_id], stats_mutex, vpls, scene, general_params, ac_params, total_samples, performed_samples, rng, ge);

        {
            std::lock_guard<std::mutex> lock(work_mutex);
            work[thread_id] = -1;
        }
    }
}

void sliceWorkerMDLC(std::vector<std::int32_t>& work, std::uint32_t thread_id, std::mutex& work_mutex, std::mutex& stats_mutex,
    const std::vector<KDTNode<ReconstructionSample>*>& slices, LightTree* light_tree, const std::vector<VPL>& tot_vpls, Scene* scene,
    const GeneralParams& general_params, const AdaptiveConstructionParams& ac_params, std::uint64_t& total_samples, std::uint64_t& performed_samples, bool ge){

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * thread_id);

    while(true){
        std::int32_t slice_id;
        {
            std::lock_guard<std::mutex> lock(work_mutex);
            slice_id = work[thread_id];
        }

        if(slice_id == -2){
            break;
        }
        else if(slice_id == -1){
            continue;
        }

        std::vector<Intersection> slice_points(slices[slice_id]->sample_indices.size());

        for(std::uint32_t i = 0; i < slice_points.size(); ++i){
            slice_points[i] = slices[slice_id]->sample(i).its;
        }

        auto vpls = light_tree->getClusteringForPoints(scene, slice_points);
        updateVPLRadii(vpls, general_params.min_dist);

        recover(slices[slice_id], stats_mutex, vpls, scene, general_params, ac_params, total_samples, performed_samples, rng, ge);
        {
            std::lock_guard<std::mutex> lock(work_mutex);
            work[thread_id] = -1;
        }
    }
}

typedef std::pair<KDTNode<ReconstructionSample>*, std::uint32_t> HWWorkUnit;

void processBatch(const std::vector<HWWorkUnit>& batch, HWShader& hw_shader, 
    std::vector<std::vector<VPL>>& vpls, std::uint32_t cluster_size, float min_dist, bool vsl){
    
    std::vector<KDTNode<ReconstructionSample>*> slices;
    std::vector<std::vector<VPL>*> slice_vpls;
    for(std::uint32_t i = 0; i < batch.size(); ++i){
        slices.push_back(batch[i].first);
        slice_vpls.push_back(&vpls[batch[i].second]);
    }

    hw_shader.renderSlices(slices, slice_vpls, cluster_size, min_dist, vsl);
}

void unoccludedHWWorker(BlockingQueue<HWWorkUnit>& input, BlockingQueue<HWWorkUnit>& output, 
    std::vector<std::vector<VPL>>& vpls, float min_dist, bool vsl, std::uint32_t cluster_size, 
    std::uint32_t batch_size){
    
    HWShader hw_shader;

    HWWorkUnit work_unit;
    std::vector<HWWorkUnit> batch;
    
    while(input.pop(work_unit)){
        batch.push_back(work_unit);

        if(batch.size() >= batch_size){
            processBatch(batch, hw_shader, vpls, cluster_size, min_dist, vsl);
            for(std::uint32_t i = 0; i < batch.size(); ++i){
                output.push(batch[i]);
            }

            batch.clear();
        }
    }

    if(batch.size() > 0){
        processBatch(batch, hw_shader, vpls, cluster_size, min_dist, vsl);
        for(std::uint32_t i = 0; i < batch.size(); ++i){
            output.push(batch[i]);
        }

        batch.clear();
    }

    output.close();
}

std::tuple<std::uint64_t, std::uint64_t> recoverHW(KDTNode<ReconstructionSample>* slice, const std::vector<VPL>& vpls, Scene* scene,
    float sample_perc, float max_sample_perc, float sample_inc, float min_dist, std::mt19937& rng, bool importance_sample, bool bin_vis, bool ge){
    Point3f slice_com_pos(0.f, 0.f, 0.f);
    Vector3f slice_com_norm(0.f);
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        slice_com_pos += slice->sample(i).its.p;
        slice_com_norm += slice->sample(i).its.shFrame.n;
    }
    slice_com_pos /= slice->sample_indices.size();
    slice_com_norm /= slice_com_norm.length();

    std::vector<std::uint32_t> material_sample_positions;
    for(std::uint32_t i = 0; i < 5; ++i){
        material_sample_positions.push_back(rand() % slice->sample_indices.size());
    }

    std::vector<std::vector<VPL>> quadranted_vpls(8);
    std::vector<std::vector<std::uint32_t>> actual_vpl_index(8);

    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        Vector3f dir = vpls[i].type == EDirectionalEmitterVPL ? Vector3f(vpls[i].its.shFrame.n) : 
            slice_com_pos - vpls[i].its.p;
        std::uint32_t quadrant = getQuadrant(dir);
        quadranted_vpls[quadrant].push_back(vpls[i]);
        actual_vpl_index[quadrant].push_back(i);
    }

    std::uint32_t total_performed_samples = 0;

    slice->visibility_coefficients.resize(slice->sample_indices.size() * vpls.size(), 0.f);

    for(std::uint32_t i = 0; i < quadranted_vpls.size(); ++i){
        if(quadranted_vpls[i].size() == 0){
            continue;
        }

        std::vector<float> cluster_contribs(quadranted_vpls[i].size(), 0.f);
        for(std::uint32_t j = 0; j < quadranted_vpls[i].size(); ++j){
            float estimated_diffuse_power = 0.f;
            for(std::uint32_t k = 0; k < material_sample_positions.size(); ++k){
                const BSDF* bsdf = slice->sample(material_sample_positions[k]).its.getBSDF();
                auto& vpl = quadranted_vpls[i][j];
                Vector3f wi = vpl.type != EDirectionalEmitterVPL ? vpl.its.p - slice_com_pos : vpl.its.shFrame.n;
                float d = wi.length();
                wi /= d;

                Spectrum diffuse_col = bsdf->getDiffuseReflectance(slice->sample(material_sample_positions[k]).its);
                diffuse_col *= dot(wi, slice_com_norm) * vpl.P  / PI;

                if(vpl.type != EDirectionalEmitterVPL){
                    float atten = std::max(d, min_dist);
                    atten *= atten;
                    diffuse_col /= atten;
                }

                if(vpl.type == ESurfaceVPL){
                    if(vpl.emitter != nullptr){
                        DirectionSamplingRecord dir(-wi);
                        diffuse_col *= vpl.emitter->evalDirection(dir, vpl.psr);
                    }
                    else if(vpl.its.getBSDF() != nullptr){
                        BSDFSamplingRecord bsdf_sample_record(vpl.its, vpl.its.toLocal(-wi));
                        diffuse_col *= vpl.its.getBSDF()->eval(bsdf_sample_record);
                    }
                    else{
                        diffuse_col *= Spectrum(std::max(0.f, dot(vpl.its.shFrame.n, -wi))) / PI;
                    }
                }

                estimated_diffuse_power += diffuse_col.getLuminance();
            }

            cluster_contribs[j] = estimated_diffuse_power;
        }
        
        std::uint32_t samples;

        if(bin_vis){
            std::vector<std::uint8_t> bv(slice->sample_indices.size() * quadranted_vpls[i].size(), 0);
            if(ge){
                samples = adaptiveMatrixReconstructionBGE(bv, scene, slice, 
                    quadranted_vpls[i], min_dist, sample_perc, max_sample_perc, sample_inc, rng, 
                    cluster_contribs);
            }
            else{
                samples = adaptiveMatrixReconstructionB(bv, scene, slice, 
                    quadranted_vpls[i], min_dist, sample_perc, max_sample_perc, sample_inc, rng, 
                    cluster_contribs);
            }

            total_performed_samples += samples;

            Properties props("independent");
            ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
            sampler->configure();
            sampler->generate(Point2i(0));

            for(std::uint32_t j = 0; j < actual_vpl_index[i].size(); ++j){
                for(std::uint32_t k = 0; k < slice->sample_indices.size(); ++k){
                    std::uint32_t idx = actual_vpl_index[i][j] * slice->sample_indices.size() + k;
                    std::uint32_t bin_vis_idx = j * slice->sample_indices.size() + k;
                    slice->visibility_coefficients[idx] = bv[bin_vis_idx] == 0 ? 0.f : 1.f;
                }
            }
        }
        else{
            Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(slice->sample_indices.size(), quadranted_vpls[i].size());

            samples = adaptiveMatrixReconstruction(mat, scene, slice, quadranted_vpls[i], min_dist, 
                sample_perc, rng, importance_sample, cluster_contribs);

            total_performed_samples += samples;

            for(std::uint32_t j = 0; j < actual_vpl_index[i].size(); ++j){
                for(std::uint32_t k = 0; k < slice->sample_indices.size(); ++k){
                    std::uint32_t idx = actual_vpl_index[i][j] * slice->sample_indices.size() + k;
                    slice->visibility_coefficients[idx] = (mat(k, j) + 1.f) / 2.f;
                }
            }
        }
        
    }

    return std::make_tuple(slice->sample_indices.size() * vpls.size(), total_performed_samples);
}

void clusterWorkerMDLC(BlockingQueue<HWWorkUnit>& input, BlockingQueue<HWWorkUnit>& output,
    std::vector<std::vector<VPL>>& vpls, const std::vector<VPL>& tot_vpls, LightTree* light_tree, Scene* scene, float sample_perc, 
    float max_sample_perc, float sample_inc, float min_dist, std::uint32_t thread_id, std::uint32_t num_threads, 
    std::mutex& barrier_mutex, std::condition_variable& barrier, std::mutex& sample_update_mutex,
    std::uint64_t& total_samples, std::uint64_t& num_samples, bool bin_vis, bool importance_sample, bool ge){
    
    HWWorkUnit work_unit;

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * thread_id);

    while(input.pop(work_unit)){
        std::vector<Intersection> slice_points(work_unit.first->sample_indices.size());

        for(std::uint32_t i = 0; i < slice_points.size(); ++i){
            slice_points[i] = work_unit.first->sample(i).its;
        }

        vpls[work_unit.second] = light_tree->getClusteringForPoints(scene, slice_points);
        //not sure if this should change to a radius union approach instead
        updateVPLRadii(vpls[work_unit.second], min_dist);

        std::uint64_t slice_samples, num_slice_sampled;

        std::tie(slice_samples, num_slice_sampled) = recoverHW(work_unit.first, vpls[work_unit.second], scene, sample_perc,
            max_sample_perc, sample_inc, min_dist, rng, importance_sample, bin_vis, ge);
        
        {
            std::lock_guard<std::mutex> lock(sample_update_mutex);
            total_samples += slice_samples;
            num_samples += num_slice_sampled;
        }

        output.push(work_unit);
    }

    static std::uint32_t threads_completed = 0;
    threads_completed++;
    
    {
        std::unique_lock<std::mutex> lock(barrier_mutex);
        barrier.wait(lock, [num_threads](){return threads_completed >= num_threads;});
    }
    barrier.notify_all();

    output.close();
}

void clusterWorkerLS(BlockingQueue<HWWorkUnit>& input, BlockingQueue<HWWorkUnit>& output,
    std::vector<std::vector<VPL>>& vpls, const std::vector<std::vector<std::uint32_t>>& cbsamp, 
    const Eigen::MatrixXf& contribs, const std::vector<std::vector<int>>& nn, 
    const std::vector<VPL>& total_vpls, std::uint32_t num_clusters, Scene* scene, 
    float sample_perc, float max_sample_perc, float sample_inc, float min_dist, std::uint32_t thread_id, std::uint32_t num_threads, std::mutex& barrier_mutex,
    std::condition_variable& barrier, std::mutex& sample_update_mutex, std::uint64_t& total_samples, std::uint64_t& num_samples,
    bool bin_vis, bool importance_sample, bool ge){

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * thread_id);
    std::uint32_t split_num_clusters = std::min(total_vpls.size(), num_clusters - cbsamp.size());

    HWWorkUnit work_unit;
    while(input.pop(work_unit)){
        std::vector<std::vector<std::uint32_t>> cbsplit = 
            clusterVPLsBySplitting(cbsamp, contribs, nn[work_unit.second], split_num_clusters, rng);
        vpls[work_unit.second] = sampleRepresentatives(contribs, total_vpls, cbsplit, rng, min_dist);

        std::uint64_t slice_samples, num_slice_sampled;
        std::tie(slice_samples, num_slice_sampled) = recoverHW(work_unit.first, vpls[work_unit.second], scene, sample_perc,
            max_sample_perc, sample_inc, min_dist, rng, importance_sample, bin_vis, ge);

        {
            std::lock_guard<std::mutex> lock(sample_update_mutex);
            total_samples += slice_samples;
            num_samples += num_slice_sampled;
        }

        output.push(work_unit);
    }

    static std::uint32_t threads_completed = 0;
    threads_completed++;
    {
        std::unique_lock<std::mutex> lock(barrier_mutex);
        barrier.wait(lock, [num_threads](){return threads_completed >= num_threads;});
    }
    
    barrier.notify_all();

    output.close();
}

std::tuple<std::uint64_t, std::uint64_t> MatrixReconstructionRenderer::renderNonHW(Scene* scene, std::uint32_t spp, const RenderJob *job,
    std::vector<float>& timings, const std::vector<KDTNode<ReconstructionSample>*>& slices, 
    std::uint32_t samples_per_slice){
    
    std::uint64_t amount_sampled = 0;
    std::uint64_t total_samples = 0;

    AdaptiveConstructionParams ac_params;
    ac_params.import_sample = adaptive_importance_sampling_;
    ac_params.bin_vis = bin_vis_;
    
    GeneralParams general_params;
    general_params.min_dist = min_dist_;
    general_params.vsl = vsl_;
    general_params.sample_perc = sample_percentage_;
    general_params.max_sample_perc = max_sample_perc_;
    general_params.sample_inc = sample_inc_;

    std::uint32_t num_cores = 16;
    std::uint32_t scheduled_slices = num_cores;

    std::vector<std::int32_t> work(num_cores);
    std::iota(work.begin(), work.end(), 0);
    std::mutex work_mutex, stats_mutex;
    std::vector<std::thread> workers;
    
    Eigen::MatrixXf cluster_contributions;
    std::vector<std::vector<std::uint32_t>> clusters_by_sampling;
    std::vector<std::vector<int>> nearest_neighbours;
    std::vector<std::vector<float>> neighbour_distances;
    std::unique_ptr<LightTree> light_tree;

    if(clustering_strategy_ == ClusteringStrategy::LS){
        ProgressReporter lscc_pr("Calculating cluster contributions", 1, job);
        std::uint32_t num_clusters = std::min(std::uint32_t(vpls_.size()), std::uint32_t(num_clusters_ * 0.3f));
        cluster_contributions = calculateClusterContributions(vpls_, slices, scene, min_dist_,
            samples_per_slice, vsl_);
        lscc_pr.finish();
        std::cout << std::endl;

        ProgressReporter lscbs_pr("Cluster by sampling", 1, job);
        clusters_by_sampling = clusterVPLsBySampling(cluster_contributions,
            num_clusters, vpls_, min_dist_);
        lscbs_pr.finish();
        std::cout << std::endl;

        ProgressReporter lsann_pr("Finding cluster neighbours", 1, job);
        //approx nearest neighbours here
        flann::Matrix<float> slice_centroids(new float[slices.size() * 3], 
            slices.size(), 3);

        for(std::uint32_t i = 0; i < slices.size(); ++i){
            Point slice_centroid(0.f, 0.f, 0.f);
            for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
                slice_centroid += slices[i]->sample(j).its.p;
            }
            slice_centroid /= slices[i]->sample_indices.size();

            float* curr_centroid = (float*)slice_centroids[i];
            curr_centroid[0] = slice_centroid.x;
            curr_centroid[1] = slice_centroid.y;
            curr_centroid[2] = slice_centroid.z;
        }

        flann::Index<flann::L2<float>> index(slice_centroids, flann::KDTreeIndexParams(4));
        index.buildIndex();

        index.knnSearch(slice_centroids, nearest_neighbours, neighbour_distances, 10 * samples_per_slice, flann::SearchParams(128));

        lsann_pr.finish();
        std::cout << std::endl;

        for(std::uint32_t i = 0; i < num_cores; ++i){
            workers.emplace_back(sliceWorkerLS, std::ref(work), i, std::ref(work_mutex), std::ref(stats_mutex), std::ref(slices), 
            std::ref(clusters_by_sampling), std::ref(cluster_contributions), std::ref(nearest_neighbours), 
            std::ref(vpls_), scene, general_params, ac_params, std::ref(total_samples), std::ref(amount_sampled), num_clusters_, ge_);
        }
    }
    else{
        ProgressReporter lt_pr("Constructing light tree", 1, job);
        light_tree = std::unique_ptr<LightTree>(new LightTree(vpls_, min_dist_, num_clusters_, 0.f, true));
        lt_pr.finish();
        std::cout << std::endl;

        for(std::uint32_t i = 0; i < num_cores; ++i){
            workers.emplace_back(sliceWorkerMDLC, std::ref(work), i, std::ref(work_mutex), std::ref(stats_mutex), std::ref(slices), 
            light_tree.get(), std::ref(vpls_), scene, general_params, ac_params, std::ref(total_samples), 
            std::ref(amount_sampled), ge_);
        }
    }
    

    std::uint32_t finished_slices = 0;
    
    ProgressReporter progress_rs("Reconstructing slices", slices.size(), job);
    while(finished_slices < slices.size()){
        std::int32_t free_worker_id = -1;
        {
            std::lock_guard<std::mutex> lock(work_mutex);
            for(std::uint32_t i = 0; i < work.size(); ++i){
                if(work[i] == -1 && workers[i].joinable()){
                    free_worker_id = i;
                    break;
                }
            }
        }

        if(free_worker_id >= 0){
            progress_rs.update(finished_slices++);

            if(scheduled_slices < slices.size())
            {
                std::lock_guard<std::mutex> lock(work_mutex);
                work[free_worker_id] = scheduled_slices++;
            }
            else{
                {
                    std::lock_guard<std::mutex> lock(work_mutex);
                    work[free_worker_id] = -2;
                }
                workers[free_worker_id].join();
            }
        }
    }
    progress_rs.finish();
    std::cout << std::endl;

    float sample_perc = (float)amount_sampled / total_samples;
    timings.push_back(sample_perc);
    std::cout << "Sample percentage: " << sample_perc << std::endl;

    return std::make_tuple(total_samples, amount_sampled);
}

std::tuple<std::uint64_t, std::uint64_t> MatrixReconstructionRenderer::renderHW(Scene* scene, std::uint32_t spp, const RenderJob *job, 
    std::vector<float>& timings, const std::vector<KDTNode<ReconstructionSample>*>& slices, 
    std::uint32_t samples_per_slice, std::uint32_t slice_size){
    auto start = std::chrono::high_resolution_clock::now();
    std::uint32_t num_workers = 15; // for my computer, change to your cores to improve performance
    std::uint32_t batch_size = 500 / (std::max(slice_size / 1000u, 1u) * std::max(1u, num_clusters_ / 4000)); //num slices to be shaded at a time

    BlockingQueue<HWWorkUnit> to_cluster;
    BlockingQueue<HWWorkUnit> to_shade(batch_size * 2);
    BlockingQueue<HWWorkUnit> to_finish;

    for(std::uint32_t i = 0; i < slices.size(); ++i){
        to_cluster.push(std::make_pair(slices[i], i));
    }
    to_cluster.close();

    std::vector<std::vector<VPL>> vpls(slices.size());

    std::thread shader(unoccludedHWWorker, std::ref(to_shade), std::ref(to_finish), std::ref(vpls), 
        min_dist_, vsl_, num_clusters_, batch_size);

    std::vector<std::thread> clusterers;

    std::vector<std::vector<std::uint32_t>> cbsamp;
    Eigen::MatrixXf cluster_contributions;
    std::vector<std::vector<int>> nearest_neighbours;
    std::vector<std::vector<float>> neighbour_distances;
    std::unique_ptr<LightTree> light_tree;

    std::mutex barrier_mutex;
    std::condition_variable barrier;
    std::mutex sample_update_mutex;
    std::uint64_t total_samples = 0;
    std::uint64_t num_samples = 0;

    if(clustering_strategy_ == ClusteringStrategy::LS){
        std::uint32_t num_clusters = std::min(std::uint32_t(vpls_.size()), std::uint32_t(num_clusters_ * 0.3f));
        std::cout << "Clusters initialized" << std::endl;
        cluster_contributions = calculateClusterContributions(vpls_, slices, scene, min_dist_,
            samples_per_slice, vsl_);

        std::cout << "Contributions calculated" << std::endl;

        cbsamp = clusterVPLsBySampling(cluster_contributions, num_clusters, vpls_, min_dist_);

        std::cout << "VPLs clustered by sampling" << std::endl;

        flann::Matrix<float> slice_centroids(new float[slices.size() * 3], 
            slices.size(), 3);

        for(std::uint32_t i = 0; i < slices.size(); ++i){
            Point slice_centroid(0.f, 0.f, 0.f);
            for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
                slice_centroid += slices[i]->sample(j).its.p;
            }
            slice_centroid /= slices[i]->sample_indices.size();

            float* curr_centroid = (float*)slice_centroids[i];
            curr_centroid[0] = slice_centroid.x;
            curr_centroid[1] = slice_centroid.y;
            curr_centroid[2] = slice_centroid.z;
        }

        flann::Index<flann::L2<float>> index(slice_centroids, flann::KDTreeIndexParams(4));
        index.buildIndex();

        index.knnSearch(slice_centroids, nearest_neighbours, neighbour_distances, 10, flann::SearchParams(128));

        std::cout << "Nearest neighbours found" << std::endl;

        for(std::uint32_t i = 0; i < num_workers; ++i){
            clusterers.emplace_back(clusterWorkerLS, std::ref(to_cluster), std::ref(to_shade), std::ref(vpls),
                std::ref(cbsamp), std::ref(cluster_contributions), std::ref(nearest_neighbours), std::ref(vpls_),
                num_clusters_, scene, sample_percentage_, max_sample_perc_, sample_inc_, min_dist_, i,
                num_workers, std::ref(barrier_mutex), std::ref(barrier), std::ref(sample_update_mutex), std::ref(total_samples), std::ref(num_samples),
                bin_vis_, adaptive_importance_sampling_, ge_);
        }
        
    }
    else{
        std::cout << "Creating light tree" << std::endl;
        light_tree = std::unique_ptr<LightTree>(new LightTree(vpls_, min_dist_, num_clusters_, 0.f, true));
        for(std::uint32_t i = 0; i < num_workers; ++i){
            clusterers.emplace_back(clusterWorkerMDLC, std::ref(to_cluster), std::ref(to_shade), std::ref(vpls), std::ref(vpls_),
                light_tree.get(), scene, sample_percentage_, max_sample_perc_, sample_inc_, min_dist_, i,
                num_workers, std::ref(barrier_mutex), std::ref(barrier), std::ref(sample_update_mutex), 
                std::ref(total_samples), std::ref(num_samples), bin_vis_, adaptive_importance_sampling_, ge_);
        }
    }

    std::uint32_t slice_counter = 0;
    ProgressReporter progress("Reconstructing slices", slices.size(), job);
    HWWorkUnit unit;
    while(to_finish.pop(unit)){
        progress.update(slice_counter++);
    }
    progress.finish();

    for(std::uint32_t i = 0; i < clusterers.size(); ++i){
        clusterers[i].join();
    }

    shader.join();

    auto end = std::chrono::high_resolution_clock::now();
    timings.push_back(std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count());

    return std::make_tuple(total_samples, num_samples);
}

bool MatrixReconstructionRenderer::render(Scene* scene, std::uint32_t spp, const RenderJob *job){
    srand(time(0));
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
    std::uint32_t samples_per_slice = samples_per_slice_;

    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, size);
    std::uint8_t* output_image = output_bitmap->getUInt8Data();
    memset(output_image, 0, output_bitmap->getBytesPerPixel() * size.x * size.y);

    ProgressReporter kdt_pr("Constructing kd-tree", 1, job);
    auto kdt_root = constructKDTree(scene, slice_size_, samples_, scene->getAABB().getBSphere().radius, std::max(10u, samples_per_slice), spp);
    kdt_pr.finish();
    std::cout << std::endl;

    std::vector<KDTNode<ReconstructionSample>*> slices;
    getSlices(kdt_root.get(), slices);

    std::uint64_t total_samples, num_samples;
    std::vector<float> timings;

    if(hw_){
        std::tie(total_samples, num_samples) = renderHW(scene, spp, job, timings, slices, samples_per_slice, slice_size_);
    }
    else{
        std::tie(total_samples, num_samples) = renderNonHW(scene, spp, job, timings, slices, samples_per_slice);
    }

    copySamplesToBuffer(output_image, samples_, size, spp);

    film->setBitmap(output_bitmap);

    return !cancel_;
}

MTS_NAMESPACE_END