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
#include "filewriter.h"
#include "lighttree.h"
#include <mitsuba/core/statistics.h>
#include <set>
#include <unistd.h>

#include "common.h"
#include "hwshader.h"
#include "blockingqueue.hpp"

MTS_NAMESPACE_BEGIN

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

                float x_jitter = /*sampler->next1D()*/ 0.5f * cell_side_len;
                float y_jitter = /*sampler->next1D()*/ 0.5f * cell_side_len;
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
                curr_sample.fully_sampled_color = Spectrum(0.f);
                curr_sample.visibility_errors = 0;

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

//Used in the proximal gradient descent version of recovery. Matrix is sparsely populated uniformly with observations
//RGB assumed, which are dealt with in separate matrices
std::vector<std::uint32_t> calculateSparseSamples(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, Eigen::MatrixXf& rmat, Eigen::MatrixXf& gmat, Eigen::MatrixXf& bmat,
    std::uint32_t num_samples, float min_dist, bool vsl){
    assert(rmat.rows() * rmat.cols() > 0 && gmat.rows() * gmat.cols() > 0 && bmat.rows() * bmat.cols() > 0);

    std::uint32_t total_samples = slice->sample_indices.size() * vpls.size();
    num_samples = std::min(num_samples, total_samples);
    std::vector<std::uint32_t> indices(total_samples);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    indices.erase(indices.begin() + num_samples, indices.end());

    Properties props("independent");
	ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < indices.size(); ++i){
        std::uint32_t light_index = indices[i] % vpls.size();
        std::uint32_t sample_index = indices[i] / vpls.size();
        const VPL& vpl = vpls[light_index];
        ReconstructionSample& sample_to_compute = slice->sample(sample_index);

        std::uint32_t num_samples;
        Spectrum lightContribution = sample(scene, sampler, sample_to_compute.its, sample_to_compute.ray, vpl, min_dist, true,
            5, false, sample_to_compute.intersected_scene, false, vsl, num_samples);

        Float r, g, b;
        lightContribution.toLinearRGB(r, g, b);
        rmat(sample_index, light_index) = r;
        gmat(sample_index, light_index) = g;
        bmat(sample_index, light_index) = b;
    }

    return indices;
}

std::vector<std::uint32_t> calculateSparseSamples(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, Eigen::MatrixXf& vmat, std::uint32_t num_samples, float min_dist){
    assert((size_t)vmat.rows() == slice->sample_indices.size() && (size_t)vmat.cols() == vpls.size());

    std::uint32_t total_samples = slice->sample_indices.size() * vpls.size();
    num_samples = std::min(num_samples, total_samples);
    std::vector<std::uint32_t> indices(total_samples);
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    indices.erase(indices.begin() + num_samples, indices.end());

    Properties props("independent");
	ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t i = 0; i < indices.size(); ++i){
        std::uint32_t light_index = indices[i] % vpls.size();
        std::uint32_t sample_index = indices[i] / vpls.size();
        const VPL& vpl = vpls[light_index];
        ReconstructionSample& sample_to_compute = slice->sample(sample_index);

        bool los = sampleVisibility(scene, sample_to_compute.its, vpl, min_dist);
        vmat(sample_index, light_index) = los ? 1.f : -1.f;
    }

    return indices;
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

void updateSliceWithMatData(const Eigen::MatrixXf& rmat, const Eigen::MatrixXf& gmat,
    const Eigen::MatrixXf& bmat, KDTNode<ReconstructionSample>* slice){
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        for(std::uint32_t j = 0; j < rmat.cols(); ++j){
            Spectrum c;
            c.fromLinearRGB(rmat(i, j), gmat(i, j), bmat(i, j));
            slice->sample(i).color += c;
        }
    }
}

void updateSliceWithMatData(const Eigen::MatrixXf& mat, KDTNode<ReconstructionSample>* slice, 
    bool visibility_only, bool recover_transpose, bool vsl, Scene* scene, const std::vector<VPL>& vpls,
    const std::vector<std::uint32_t>& actual_vpl_indices, float min_dist, bool fully_compute_occlusion){

    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    std::uint32_t total_samples = slice->sample_indices.size();
    for(std::uint32_t i = 0; i < total_samples; ++i){
        if(visibility_only){
            if(slice->sample(i).its.isEmitter()){
                Spectrum emitter_col = slice->sample(i).its.Le(-slice->sample(i).ray.d);
                slice->sample(i).color = emitter_col;
                if(fully_compute_occlusion){
                    slice->sample(i).fully_sampled_color = emitter_col;
                }
            }
            else{
                for(std::uint32_t j = 0; j < actual_vpl_indices.size(); ++j){
                    float coeff = ((recover_transpose ? mat(j, i) : mat(i, j)) + 1.f) / 2.f;
                    //float coeff = mat(i, j) > 0.f ? 1.f : 0.f;
                    std::uint32_t idx = actual_vpl_indices[j];
                    slice->sample(i).color += slice->sample(i).unoccluded_samples[idx] * coeff;
                    if(fully_compute_occlusion){
                        bool los = sampleVisibility(scene, slice->sample(i).its, vpls[idx], min_dist);
                        float correct_coeff = los ? 1.f : 0.f;
                        slice->sample(i).fully_sampled_color += slice->sample(i).unoccluded_samples[idx] * correct_coeff;
                    }
                }
            }
        }
        else{
            for(std::uint32_t j = 0; j < actual_vpl_indices.size(); ++j){
                float r = recover_transpose ? mat(j * 3, i) : mat(i * 3, j);
                float g = recover_transpose ? mat(j * 3 + 1, i) : mat(i * 3 + 1, j);
                float b = recover_transpose ? mat(j * 3 + 2, i) : mat(i * 3 + 2, j);

                Spectrum c;
                c.fromLinearRGB(r, g, b);

                slice->sample(i).color += c;
            }
        }
    }
}

void copySamplesToBuffer(std::uint8_t* output_image, const std::vector<ReconstructionSample>& samples, Vector2i image_size,
    std::uint32_t spp){
    std::unordered_map<std::uint32_t, Spectrum> output;
    //#pragma omp parallel for
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

void constructSvdBuffers(std::vector<Spectrum>& col_ranks, const std::vector<KDTNode<ReconstructionSample>*>& slices, 
    std::uint32_t spp, std::uint32_t image_width){
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        auto slice = slices[i];

        for(std::uint32_t j = 0; j < slice->sample_indices.size(); ++j){
            std::uint32_t index = slice->sample(j).image_x + slice->sample(j).image_y * image_width;
            
            float r, g, b;
            std::tie(r, g, b) = floatToRGB(slice->rank_ratio2);
            Spectrum rank2_col;
            rank2_col.fromLinearRGB(r, g, b);
            col_ranks[index] += rank2_col / spp;
        }
    }
}

void constructStatBuffers(std::vector<Spectrum>& recovered, std::vector<Spectrum>& full, std::vector<std::vector<Spectrum>>& ranks,
    std::vector<Spectrum>& samples, std::vector<Spectrum>& slice_cols, std::vector<Spectrum>& visibility_errors, std::uint64_t& total_error,
    const std::vector<KDTNode<ReconstructionSample>*>& slices, std::uint32_t spp, std::uint32_t image_width,
    std::uint32_t columns){

    for(std::uint32_t i = 0; i < slices.size(); ++i){
        auto slice = slices[i];
        int y = rand() % 100 + 100;
        int u = rand() % 255;
        int v = rand() % 255;

        int sb = 1.164 * (y - 16) + 2.018 * (u - 128);
        int sg = 1.164 * (y - 16) - 0.813 * (v - 128) - 0.391 * (u - 128);
        int sr = 1.164 * (y - 16) + 1.596 * (v - 128);

        Spectrum slice_col;
        slice_col.fromSRGB(sr / 255.f, sg / 255.f, sb / 255.f);

        for(std::uint32_t j = 0; j < slice->sample_indices.size(); ++j){
            std::uint32_t index = slice->sample(j).image_x + slice->sample(j).image_y * image_width;

            recovered[index] += slice->sample(j).color / spp;
            full[index] += slice->sample(j).fully_sampled_color / spp;
            float r, g, b;
            
            for(std::uint32_t k = 0; k < 8; ++k){
                std::tie(r, g, b) = floatToRGB(slice->rank_ratio[k]);
                Spectrum rank_col;
                rank_col.fromLinearRGB(r, g, b);
                ranks[k][index] += rank_col / spp;
            }

            std::tie(r, g, b) = floatToRGB(slice->sample_ratio);
            Spectrum sample_ratio_col;
            sample_ratio_col.fromLinearRGB(r, g, b);
            samples[index] += sample_ratio_col / spp;

            slice_cols[index] += slice_col / spp;

            float ratio = float(slice->sample(j).visibility_errors) / float(columns);
            total_error += slice->sample(j).visibility_errors;
            std::tie(r, g, b) = floatToRGB(ratio);
            Spectrum visibility_err;
            visibility_err.fromLinearRGB(r, g, b);
            visibility_errors[index] += visibility_err / spp;
        }
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

    for(size_t i = 0; i < sampled_indices.size(); ++i){
        std::uint32_t vpl_index = recover_transpose ? sampled_indices[i] : col;
        std::uint32_t sample_index = recover_transpose ? col : sampled_indices[i];
        const VPL& vpl = vpls[vpl_index];
        ReconstructionSample& scene_sample = slice->sample(sample_index);

        if(visibility_only){
            bool los = sampleVisibility(scene, scene_sample.its, vpl, min_dist);
            mat(i, 0) = los ? 1.f : -1.f;
        }
        else{
            Properties props("independent");
            ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
            sampler->configure();
            sampler->generate(Point2i(0));

            std::uint32_t num_samples;
            Spectrum lightContribution = sample(scene, sampler, scene_sample.its, scene_sample.ray, vpl, 
                min_dist, true, 5, false, scene_sample.intersected_scene, false, vsl, num_samples);

            Float r, g, b;
            lightContribution.toLinearRGB(r, g, b);
            mat(i * 3, 0) = r;
            mat(i * 3 + 1, 0) = g;
            mat(i * 3 + 2, 0) = b;
        }
    }

    return sampled_indices;
}

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

std::vector<std::uint32_t> sampleColBWithLeading(Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, std::uint32_t col, float min_dist, std::uint32_t num_samples, 
    std::mt19937& rng, std::unordered_map<std::uint32_t, std::uint8_t>& sampled_vals, const std::vector<float>& leading_probabilities, 
    const std::vector<float>& non_leading_probabilities, const std::vector<std::uint32_t>& leading_indices, 
    const std::vector<std::uint32_t>& basis_indices, float leading_perc, const std::vector<std::uint32_t>& sample_set, bool resample){

    std::uint32_t max_samples = slice->sample_indices.size();
    assert(num_samples <= max_samples);

    std::vector<std::uint32_t> sample_indices;

    if(resample){
        std::uint32_t num_leading_samples = std::max(std::uint32_t(leading_indices.size()), std::min(1u, std::uint32_t(leading_perc * num_samples)));
        if(num_leading_samples == leading_indices.size()){
            for(std::uint32_t i = 0; i < leading_indices.size(); ++i){
                sample_indices.push_back(basis_indices[leading_indices[i]]);
            }
        }
        else{
            sample_indices = importanceSample(num_leading_samples, rng, leading_probabilities);
        }

        std::uint32_t remaining_samples = num_samples - num_leading_samples;
        std::vector<std::uint32_t> nonleading_sample_indices;
        if(remaining_samples > 0){
            nonleading_sample_indices = importanceSample(remaining_samples, rng, non_leading_probabilities)
        }

        sample_indices.insert(sample_indices.end(), nonleading_sample_indices.begin(), nonleading_sample_indices.end());
    }
    else{
        sample_indices = sample_set;
    }

    for(size_t i = 0; i < sampled_indices.size(); ++i){
        const VPL& vpl = vpls[col];
        ReconstructionSample& scene_sample = slice->sample(sampled_indices[i]);

        if(sampled_vals.find(sampled_indices[i]) == sampled_vals.end()){
            sampled_vals[sampled_indices[i]] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;    
        }
    }

    return sampled_indices;
}

std::vector<std::vector<std::uint8_t>> gf2elim(const std::vector<std::vector<std::uint8_t>>& basis, std::vector<std::uint32_t>& reduced_basis_rows, 
    std::vector<std::uint32_t>& leading_pos){
    
    reduced_basis_rows.clear();
    leading_pos.clear();
    std::vector<std::vector<std::uint8_t>> reduced_basis;
    
    if(basis.size() == 0){
        return reduced_basis;
    }

    reduced_basis_rows.clear();

    //find nonzero rows, we only reduce those
    for(std::uint32_t i = 0; i < basis[0].size(); ++i){
        std::uint8_t one_exists = 0;

        for(std::uint32_t j = 0; j < basis.size(); ++j){
            one_exists |= basis[j][i];
        }

        if(one_exists){
            reduced_basis_rows.push_back(i);
        }
    }

    if(reduced_basis_rows.size() == 0){
        return reduced_basis;
    }

    //swap rows and columns because we perform elimination on a transposed matrix
    std::uint32_t cols = basis.size();
    std::uint32_t rows = reduced_basis_rows.size();

    reduced_basis.resize(cols);

    //copy nonzero info to separate matrix
    for(std::uint32_t i = 0; i < cols; ++i){
        reduced_basis[i].resize(rows);
        for(std::uint32_t j = 0; j < rows; ++j){
            std::uint32_t idx = reduced_basis_rows[j];
            reduced_basis[i][j] = basis[i][idx];
        }
    }


    //actual gaussian elimination algorithm for gf2, should be noted that we are reducing the transpose of the above matrix, so all rows and columns 
    //are inherently changed
    std::uint32_t curr_pivot = 0;

    while(curr_pivot < rows && curr_pivot < cols){
        std::uint32_t first_nonzero_idx = curr_pivot;
        for(std::uint32_t i = curr_pivot; i < reduced_basis.size(); ++i){
            if(reduced_basis[i][curr_pivot]){
                first_nonzero_idx = i;
                break;
            }
        }

        //swap rows
        if(first_nonzero_idx != curr_pivot){
            auto temp = reduced_basis[curr_pivot];
            reduced_basis[curr_pivot] = reduced_basis[first_nonzero_idx];
            reduced_basis[first_nonzero_idx] = temp;
        }

        std::vector<std::uint8_t> aijn(cols - curr_pivot);
        for(std::uint32_t i = curr_pivot; i < cols; ++i){
            aijn[i - curr_pivot] = reduced_basis[curr_pivot][i];
        }

        std::vector<std::uint8_t> c(rows);
        for(std::uint32_t i = 0; i < rows; ++i){
            c[i] = reduced_basis[i][curr_pivot];
        }

        c[0] = 0; //dont self xor pivot

        //outer product
        std::vector<std::vector<std::uint8_t>> outer_prod(rows);

        for(std::uint32_t i = 0; i < c.size(); ++i){
            outer_prod.resize(aijn.size());
            for(std::uint32_t j = 0; j < aijn.size(); ++j){
                outer_prod[i][j] = c[i] & aijn[j];
            }
        }

        //xor
        for(std::uint32_t i = 0; i < outer_prod.size(); ++i){
            for(std::uint32_t j = 0; j < outer_prod[i].size(); ++j){
                reduced_basis[i][curr_pivot + j] ^= outer_prod[i][j];
            }
        }

        curr_pivot++;
    }

    //get all leading positions, will be needed when obtaining coefficients
    std::vector<std::vector<std::uint8_t>> nonzero_reduced;
    leading_pos.clear();

    for(std::uint32_t i = 0; i < reduced_basis.size(); ++i){
        int nonzero_pos = -1;
        for(std::uint32_t j = 0; j < reduced_basis[i].size(); ++j){
            if(reduced_basis[i][j]){
                nonzero_pos = j;
                break;
            }
        }

        if(nonzero_pos >= 0){
            nonzero_reduced.push_back(reduced_basis[i]);
            leading_pos.push_back(nonzero_pos);
        }
    }

    return nonzero_reduced;
}

//basis is gaussian eliminated version of actual basis after removing all zero rows, leading indices must be sorted
bool gereconstruct(std::unordered_map<std::uint32_t, std::uint8_t>& sampled, const std::vector<std::vector<std::uint8_t>>& reduced_basis, 
    const std::vector<std::uint32_t>& basis_indices, const std::vector<std::uint32_t>& leading_indices, std::uint32_t rows){

    //no basis columns, can't reconstruct
    if(reduced_basis.size() == 0){
        return false;
    }

    std::vector<std::uint32_t> one_counts(basis_indices.size(), 0);
    std::vector<std::uint32_t> basis_to_consider;

    for(std::uint32_t i = 0; i < leading_indices.size(); ++i){
        std::uint32_t actual_index = basis_indices[leading_indices[i]];

        //only consider basis if it's pivot has been sampled
        if(sampled.find(actual_index) != sampled.end()){
            bool even = (one_counts[leading_indices[i]] & 1) == 0;

            //check if need to consider basis, if yes update tally, this works because the matrix is at least in echelon form after
            //gaussian elimination, thus we know the column will be zero from now onwards
            if((sampled[actual_index] == 1 && even) || (sampled[actual_index] == 0) && !even){
                basis_to_consider.push_back(i);
                
                for(std::uint32_t j = 0; j < one_counts.size(); ++j){
                    one_counts[i] += reduced_basis[i][j];
                }
            }
        }
    }

    std::vector<std::uint32_t> reconstructed(rows, 0);

    for(std::uint32_t i = 0; i < one_counts.size(); ++i){
        std::uint32_t curr_idx = basis_indices[i];
        
        reconstructed[curr_idx] = one_counts[i] & 1;
    }

    bool matching = true;
    for(auto iter = sampled.begin(); iter != sampled.end(); ++iter){
        std::uint32_t idx = iter->first;
        if(iter->second != reconstructed[idx]){
            matching = false;
            break;
        }
    }

    if(matching){
        for(std::uint32_t i = 0; i < rows; ++i){
            sampled[i] = reconstructed[i];
        }
    }

    return matching;
}

std::uint32_t adaptiveMatrixReconstructionBGE(
    std::vector<std::uint8_t>& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc, float max_sample_perc, float verification_inc, 
    std::mt19937& rng, std::uint32_t& basis_rank, const std::vector<float>& col_estimations){

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
    std::vector<std::uint32_t> basis_indices;
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

    bool regenerate_sample_indices = false;

    //The actual adaptive matrix recovery algorithm
    for(std::uint32_t i = 0; i < order.size(); ++i){
        sample_omega.clear();
        std::uint32_t samples_for_col = 0;

       if(basis.size() > 0){
            sampled = sampleColBWithLeading(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega,
                leading_probabilities, non_leading_probabilities, leading_indices, basis_indices, 0.25f, sampled, regenerate_sample_indices);
            regenerate_sample_indices = false;

            full_col_sampled = false;

            samples_for_col = num_samples;

            if(sample_omega.size() == col_to_add.size()){
                for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                    col_to_add[j] = sample_omega[j];   
                }
            }
            else{
                std::unordered_map<std::uint32_t, std::uint8_t> reconstructed = sample_omega; 
                if(gereconstruct(reconstructed, reduced_basis, basis_indices, leading_indices, num_rows)){
                    for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                        col_to_add[j] = sample_omega[j];   
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
            reduced_basis = gf2elim(basis, basis_indices, leading_indices);

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
            for(std::uint32_t i = 0; i < leading_indices.size(); ++i){
                std::uint32_t idx = basis_indices[leading_indices[i]];
                non_leading_probabilities[idx] = 0.f;
                leading_probabilities[idx] = probabilities[idx];
            }
        }

        std::uint32_t offset = order[i] * num_rows;
        std::copy(col_to_add.begin(), col_to_add.end(), mat.begin() + offset);
        total_samples += samples_for_col;

        regenerate_sample_indices = true;
    }
    
    basis_rank = basis.size();

    return total_samples;
}

std::tuple<std::uint32_t, float, float, float, float> adaptiveMatrixReconstruction(Eigen::MatrixXf& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc,
    std::mt19937& rng, bool visibility_only, bool recover_transpose, bool importance_sample, bool force_resample,
    std::uint32_t& basis_rank, bool vsl, const std::vector<float>& col_estimations, bool gather_stats, bool show_svd,
    std::vector<float>& singular_values){

    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);

    std::uint32_t total_rows = recover_transpose ? vpls.size() : slice->sample_indices.size();
    std::uint32_t num_rows = visibility_only ? total_rows : total_rows * 3;
    std::uint32_t num_cols = recover_transpose ? slice->sample_indices.size() : vpls.size();
    
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

    float time_for_svd = 0.f;
    float time_for_reconstruct = 0.f;
    float time_for_sampling = 0.f;
    float time_other = 0.f;
    bool full_col_sampled = false;

    //The actual adaptive matrix recovery algorithm
    for(std::uint32_t i = 0; i < order.size(); ++i){
        
        std::uint32_t samples_for_col = 0;

        //if the basis is not empty, we can try reproject, otherwise a full sample is required to populate the basis
        if(q.cols() > 0){
            std::uint32_t expected_omega_rows = visibility_only ? num_samples : num_samples * 3;
            if(sample_omega.cols() != expected_omega_rows){
                sample_omega.resize(expected_omega_rows, 1);
            }

            //we may want to regenerate the sample indices for a variety of reasons, in which case the indices are generated and
            //the pseudoinverse is recalculated
            if(sample_omega.rows() != reconstructed.rows() && (full_col_sampled || force_resample)){
                full_col_sampled = false;
                auto sample_t = std::chrono::high_resolution_clock::now();
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    true, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
                auto sample_end = std::chrono::high_resolution_clock::now();
                time_for_sampling += std::chrono::duration_cast<std::chrono::duration<float>>(sample_end - sample_t).count();

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

                
                auto svd_start = std::chrono::high_resolution_clock::now();
                auto svd = q_omega.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto sv = svd.singularValues();
                Eigen::MatrixXf singular_val_inv = Eigen::MatrixXf::Zero(sv.size(), sv.size());

                for(std::uint32_t j = 0; j < sv.size(); ++j){
                    singular_val_inv(j, j) = sv(j) < 1e-5f ? 0.f : 1.f / sv(j);
                }

                q_omega_pseudoinverse = svd.matrixV() * singular_val_inv * svd.matrixU().transpose();
                auto svd_end = std::chrono::high_resolution_clock::now();
                time_for_svd += std::chrono::duration_cast<std::chrono::duration<float>>(svd_end - svd_start).count();
                auto other_start = std::chrono::high_resolution_clock::now();
                qq_omega_pseudoinverse = q * q_omega_pseudoinverse;
                auto other_end = std::chrono::high_resolution_clock::now();
                time_other += std::chrono::duration_cast<std::chrono::duration<float>>(other_end - other_start).count();
            }
                //no new direction was added so no need to regenerate sample indices
            else{
                auto sample_t = std::chrono::high_resolution_clock::now();
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, num_samples, rng, sample_omega, sampled, 
                    false, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
                auto sample_end = std::chrono::high_resolution_clock::now();
                time_for_sampling += std::chrono::duration_cast<std::chrono::duration<float>>(sample_end - sample_t).count();
            }
            samples_for_col = num_samples;

            auto recon_start = std::chrono::high_resolution_clock::now();
            //no need to reconstruct if full sample
            reconstructed = sample_omega.rows() == reconstructed.rows() ? sample_omega : qq_omega_pseudoinverse * sample_omega;           
            auto recon_end = std::chrono::high_resolution_clock::now();

            time_for_reconstruct += std::chrono::duration_cast<std::chrono::duration<float>>(recon_end - recon_start).count();
            
            float d = 0;
            for(std::uint32_t j = 0; j < sampled.size(); ++j){
                d += std::abs(reconstructed(sampled[j], 0) - sample_omega(j, 0));
            }

            float largest = col_estimations[order[0]];
            float curr = col_estimations[order[i]];
            float thresh = (largest - curr) / largest * float(sampled.size()) * 0.1f;
            thresh = 1e-3f;//std::max(1e-3f, thresh);
            //sampled values can't be reconstructed accurately so fully sample
            if(d > thresh){
                //std::cout << d << std::endl;
                auto sample_t = std::chrono::high_resolution_clock::now();
                sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, 
                    reconstructed, sampled, true, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
                samples_for_col = total_rows;
                auto sample_end = std::chrono::high_resolution_clock::now();
                time_for_sampling += std::chrono::duration_cast<std::chrono::duration<float>>(sample_end - sample_t).count();

                q.conservativeResize(q.rows(), q.cols() + 1);
                q.col(q.cols() - 1) = reconstructed;
                
                full_col_sampled = true;
            }
        }
        else{
            auto sample_t = std::chrono::high_resolution_clock::now();
            sampled = sampleCol(scene, slice, vpls, order[i], min_dist, total_rows, rng, reconstructed, 
                sampled, true, visibility_only, recover_transpose, importance_sample, probabilities, vsl);
            auto sample_end = std::chrono::high_resolution_clock::now();
            time_for_sampling += std::chrono::duration_cast<std::chrono::duration<float>>(sample_end - sample_t).count();
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

    if(show_svd){
        auto svd = mat.bdcSvd();

        for(basis_rank = 0; basis_rank < svd.singularValues().size(); ++basis_rank){
            if(svd.singularValues()(basis_rank) < 1e-3f){
                break;
            }
        }

        singular_values.resize(svd.singularValues().size());
        for(std::uint32_t i = 0; i < svd.singularValues().size(); ++i){
            singular_values[i] = svd.singularValues()(i);
        }
    }
    else basis_rank = q.cols();

    return std::make_tuple(total_samples, time_for_svd, time_for_reconstruct, time_for_sampling, time_other);
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

std::tuple<std::uint32_t, bool, std::vector<std::uint32_t>> computeMinHammingErr(const std::vector<std::vector<std::uint8_t>>& cols, 
    const std::unordered_map<std::uint32_t, std::uint8_t>& sampled_vals, std::mt19937& rng){
    
    std::uint32_t min_dist = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::int32_t> indices;

    for(std::uint32_t i = 0; i < cols.size(); ++i){
        std::uint32_t herr = 0;
        std::uint32_t oerr = 0;

        for(auto iter = sampled_vals.begin(); iter != sampled_vals.end(); ++iter){
            if(iter->second == cols[i][iter->first]){
                oerr++;
            }
            else{
                herr++;
            }
        }

        if(herr <= min_dist){
            if(herr < min_dist){
                indices.clear();
                min_dist = herr;
            }

            indices.push_back(i + 1);
        }

        if(oerr <= min_dist){
            if(oerr < min_dist){
                indices.clear();
                min_dist = oerr;
            }

            indices.push_back(-(i + 1));
        }
    }

    std::uniform_int_distribution<std::uint32_t> select_col(0, indices.size() - 1);
    std::uint32_t sel = select_col(rng);
    std::uint32_t basis_idx = std::abs(indices[sel]) - 1;
    bool flip = indices[sel] < 0;

    std::vector<std::uint32_t> error_indices;
    for(auto iter = sampled_vals.begin(); iter != sampled_vals.end(); ++iter){
        std::uint8_t val = flip ? (cols[basis_idx][iter->first] + 1) % 2 : cols[basis_idx][iter->first];
        if(val != iter->second){
            error_indices.push_back(iter->first);
        }
    }

    return std::make_tuple(basis_idx, flip, error_indices);
}

std::uint32_t recursiveComplete(Scene* scene, KDTNode<ReconstructionSample>* slice, float min_dist, const VPL& vpl, 
    OTN<ReconstructionSample>* curr_octreenode, std::unordered_map<std::uint32_t, std::uint8_t>& sample_omega, 
    const std::vector<std::uint8_t>& basis_col, bool flip_basis, const std::vector<std::uint32_t>& incorrect_indices, 
    const std::vector<std::uint32_t> sampled_indices){
    std::uint32_t samples_taken = 0;

    //no error in subsection that has been sampled
    if(incorrect_indices.size() == 0){
        if(sampled_indices.size() > 0){
            for(std::uint32_t i = 0; i < curr_octreenode->sample_indices.size(); ++i){
                std::uint32_t idx = curr_octreenode->sample_indices[i];
                if(sample_omega.find(idx) == sample_omega.end()){
                    sample_omega[idx] = flip_basis ? (basis_col[idx] + 1) % 2 : basis_col[idx];
                }
            }

        }
        //section has not been sampled, fully sample to be safe as this only happens when there is an error in one of its siblings
        else{
            for(std::uint32_t i = 0; i < curr_octreenode->sample_indices.size(); ++i){
                std::uint32_t idx = curr_octreenode->sample_indices[i];
                if(sample_omega.find(idx) == sample_omega.end()){
                    ReconstructionSample& scene_sample = slice->sample(idx);
                    sample_omega[idx] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;
                    samples_taken++;
                }
            }
        }

        return samples_taken;
    }

    //some error

    std::uint32_t num_children = 0;
    for(std::uint32_t i = 0; i < curr_octreenode->children.size(); ++i){
        if(curr_octreenode->children[i] != nullptr){
            num_children++;
        }
    }

    //recurse if has children
    if(num_children > 0){
        std::vector<std::vector<std::uint32_t>> children_incorrect_indices(8);
        std::vector<std::vector<std::uint32_t>> children_sampled_indices(8);
        Vector3f midpoint = (curr_octreenode->bb.first + curr_octreenode->bb.second) / 2.f;

        for(std::uint32_t i = 0; i < incorrect_indices.size(); ++i){
            Vector3f p(slice->sample(incorrect_indices[i]).its.p);
            std::uint32_t child_idx = curr_octreenode->getChildIndex(p, midpoint);
            children_incorrect_indices[child_idx].push_back(incorrect_indices[i]);
        }

        for(std::uint32_t i = 0; i < sampled_indices.size(); ++i){
            Vector3f p(slice->sample(sampled_indices[i]).its.p);
            std::uint32_t child_idx = curr_octreenode->getChildIndex(p, midpoint);
            children_sampled_indices[child_idx].push_back(sampled_indices[i]);
        }

        for(std::uint32_t i = 0; i < 8; ++i){
            if(curr_octreenode->children[i] != nullptr){
                samples_taken += recursiveComplete(scene, slice, min_dist, vpl, 
                    curr_octreenode->children[i].get(), sample_omega, basis_col, flip_basis, children_incorrect_indices[i], children_sampled_indices[i]);
            }
            else{
                if(children_incorrect_indices[i].size() > 0 || children_sampled_indices[i].size() > 0){
                    std::cout << "This should never happen" << std::endl;
                }
            }
        }
    }
    //no children but has errors, fully sample
    else{
        for(std::uint32_t i = 0; i < curr_octreenode->sample_indices.size(); ++i){
            std::uint32_t idx = curr_octreenode->sample_indices[i];
            ReconstructionSample& scene_sample = slice->sample(idx);

            if(sample_omega.find(idx) == sample_omega.end()){
                sample_omega[idx] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;
                samples_taken++;
            }
        }
    }

    return samples_taken;
    
}

std::uint32_t recursiveComplete(Scene* scene, KDTNode<ReconstructionSample>* slice, float min_dist, const VPL& vpl, 
    OTN<ReconstructionSample>* curr_octreenode, std::unordered_map<std::uint32_t, std::uint8_t>& sample_omega, 
    const std::vector<std::uint8_t>& basis_col, bool flip_basis, const std::vector<std::uint32_t>& incorrect_indices, 
    const std::unordered_set<std::uint32_t>& already_sampled, float validation_samples, bool& fully_sampled){
    //no error in subsection
    if(incorrect_indices.size() == 0){
        for(std::uint32_t i = 0; i < curr_octreenode->sample_indices.size(); ++i){
            std::uint32_t idx = curr_octreenode->sample_indices[i];
            if(sample_omega.find(idx) == sample_omega.end()){
                sample_omega[idx] = flip_basis ? (basis_col[idx] + 1) % 2 : basis_col[idx];
            }
        }

        fully_sampled = false;

        return 0;
    }

    std::uint32_t num_children = 0;
    for(std::uint32_t i = 0; i < curr_octreenode->children.size(); ++i){
        if(curr_octreenode->children[i] != nullptr){
            num_children++;
        }
    }

    std::uint32_t samples_taken = 0;

    //recurse if has children
    if(num_children > 0){
        std::vector<std::vector<std::uint32_t>> children_incorrect_indices(8); 
        Vector3f midpoint = (curr_octreenode->bb.first + curr_octreenode->bb.second) / 2.f;

        for(std::uint32_t i = 0; i < incorrect_indices.size(); ++i){
            Vector3f p(slice->sample(incorrect_indices[i]).its.p);
            std::uint32_t child_idx = curr_octreenode->getChildIndex(p, midpoint);
            children_incorrect_indices[child_idx].push_back(incorrect_indices[i]);
        }

        bool childbranch_fs = false;
        bool childbranch_ns = false;
        std::vector<std::uint32_t> valid_indices;

        for(std::uint32_t i = 0; i < 8; ++i){
            if(curr_octreenode->children[i] != nullptr){
                bool child_sampled;
                samples_taken = recursiveComplete(scene, slice, min_dist, vpl, 
                    curr_octreenode->children[i].get(), sample_omega, basis_col, flip_basis, children_incorrect_indices[i], already_sampled,
                    validation_samples, child_sampled);
                
                childbranch_fs |= child_sampled;
                childbranch_ns |= !child_sampled;

                if(!child_sampled){
                    for(std::uint32_t j = 0; j < curr_octreenode->children[i]->sample_indices.size(); ++j){
                        std::uint32_t idx = curr_octreenode->children[i]->sample_indices[j];
                        if(already_sampled.find(idx) == already_sampled.end()){
                            valid_indices.push_back(idx);
                        }
                    }
                }
            }
            else{
                if(children_incorrect_indices[i].size() > 0){
                    std::cout << "This should never happen" << std::endl;
                }
            }
        }

        if(valid_indices.size() == 0){
            fully_sampled = true;
        }
        //if some children have been fully sampled and others not, verify with sparse samples and if fail, then fully sample
        else if(childbranch_fs && childbranch_ns){
            std::random_shuffle(valid_indices.begin(), valid_indices.end());
            std::uint32_t num_validation_samples = std::max(1u, 
                std::min(std::uint32_t(valid_indices.size()), std::uint32_t(curr_octreenode->sample_indices.size() * validation_samples)));

            fully_sampled = false;
            for(std::uint32_t i = 0; i < num_validation_samples; ++i){
                std::uint32_t idx = valid_indices[i];

                ReconstructionSample& scene_sample = slice->sample(idx);
                std::uint8_t result = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;

                if(result != sample_omega[idx]){
                    fully_sampled = true;
                    break;
                }
            }

            if(fully_sampled){
                for(std::uint32_t i = 0; i < valid_indices.size(); ++i){
                    std::uint32_t idx = valid_indices[i];

                    if(already_sampled.find(idx) == already_sampled.end()){
                        ReconstructionSample& scene_sample = slice->sample(idx);
                        sample_omega[idx] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;
                        samples_taken++;
                    }
                }
            }
        }
    }
    //no children but has errors, fully sample
    else{
        for(std::uint32_t i = 0; i < curr_octreenode->sample_indices.size(); ++i){
            std::uint32_t idx = curr_octreenode->sample_indices[i];
            ReconstructionSample& scene_sample = slice->sample(idx);

            if(sample_omega.find(idx) == sample_omega.end()){
                sample_omega[idx] = sampleVisibility(scene, scene_sample.its, vpl, min_dist) ? 1 : 0;
                samples_taken++;
            }
        }

        fully_sampled = true;
    }

    return samples_taken;
    
}

std::uint32_t adaptiveMatrixReconstructionBRecursive(
    std::vector<std::uint8_t>& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc, float max_sample_perc, float verification_inc, 
    std::mt19937& rng, std::uint32_t& basis_rank, const std::vector<float>& col_estimations){

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

            if(sample_omega.size() == col_to_add.size()){
                for(std::uint32_t j = 0; j < col_to_add.size(); ++j){
                    col_to_add[j] = sample_omega[j];   
                }
            }
            else{
                std::uint32_t basis_index;
                bool flip;
                std::vector<std::uint32_t> error_indices;

                std::tie(basis_index, flip, error_indices) = computeMinHammingErr(basis, sample_omega, rng);

                if(error_indices.size() == 0){
                    //opposite direction so we flip
                    col_to_add = basis[basis_index];
                    if(flip){
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
                    samples_for_col += recursiveComplete(scene, slice, min_dist, vpls[order[i]], 
                        slice->octree_root.get(), sample_omega, basis[basis_index], flip, error_indices, 
                        sampled);

                    /*std::unordered_set<std::uint32_t> sampled_set(sampled.begin(), sampled.end());

                    bool sampled;
                    samples_for_col += recursiveComplete(scene, slice, min_dist, vpls[order[i]], 
                        slice->octree_root.get(), sample_omega, basis[basis_index], flip, error_indices, 
                        sampled_set, 0.05, sampled);*/

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
    
    basis_rank = basis.size();

    return total_samples;
}

std::uint32_t adaptiveMatrixReconstructionB(
    std::vector<std::uint8_t>& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc, float max_sample_perc, float verification_inc, 
    std::mt19937& rng, std::uint32_t& basis_rank, const std::vector<float>& col_estimations){

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
    
    basis_rank = basis.size();

    return total_samples;
}

//singular value thresholding algorithm for the non-adaptive version. Can use the RedSVD truncated svd, or just normal svd
void svt(Eigen::MatrixXf& reconstructed_matrix, const Eigen::MatrixXf& lighting_matrix, float step_size, 
    float tolerance, float tau, std::uint32_t max_iterations, const std::vector<std::uint32_t>& sampled_indices,
    bool truncated){

    std::uint32_t k0 = tau / (step_size * lighting_matrix.norm()) + 1.5f; //extra .5 for rounding in case of float error
    Eigen::MatrixXf y = step_size * (float)k0 * lighting_matrix;
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
    bool truncated, bool show_slices, bool vsl, bool gather_stat_images, bool show_svd,
    ClusteringStrategy clustering_strategy, float error_scale, bool hw, bool bin_vis, std::uint32_t num_clusters,
    std::uint32_t samples_per_slice, float max_sample_perc, float ver_inc, bool show_vmat) : 
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
        gather_stat_images_(gather_stat_images),
        show_svd_(show_svd),
        clustering_strategy_(clustering_strategy),
        error_scale_(error_scale),
        hw_(hw),
        bin_vis_(bin_vis),
        num_clusters_(num_clusters),
        samples_per_slice_(samples_per_slice),
        sample_inc_(ver_inc), 
        max_sample_perc_(max_sample_perc),
        show_vmat_(show_vmat),
        cancel_(false){
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
    gather_stat_images_(other.gather_stat_images_),
    show_svd_(other.show_svd_),
    clustering_strategy_(other.clustering_strategy_),
    error_scale_(other.error_scale_),
    hw_(other.hw_),
    bin_vis_(other.bin_vis_),
    num_clusters_(other.num_clusters_),
    samples_per_slice_(other.samples_per_slice_),
    sample_inc_(other.sample_inc_), 
    max_sample_perc_(other.max_sample_perc_),
    show_vmat_(other.show_vmat_),
    cancel_(other.cancel_){
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
        gather_stat_images_ = other.gather_stat_images_;
        show_svd_ = other.show_svd_;
        clustering_strategy_ = other.clustering_strategy_;
        error_scale_ = other.error_scale_;
        hw_ = other.hw_;
        bin_vis_ = other.bin_vis_;
        num_clusters_ = other.num_clusters_;
        samples_per_slice_ = other.samples_per_slice_;
        sample_inc_ = other.sample_inc_; 
        max_sample_perc_ = other.max_sample_perc_;
        show_vmat_ = other.show_vmat_;
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

struct AdaptiveConstructionParams{
    bool vis_only;
    bool rec_trans;
    bool import_sample;
    bool force_resample;
    bool bin_vis;
};

struct SVTParams{
    float tolerance;
    float tau;
    std::uint32_t max_iter;
    bool trunc;
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

std::tuple<float, float, float, float, float, float> recover(KDTNode<ReconstructionSample>* slice, std::mutex& stats_mutex, const std::vector<VPL>& vpls, 
    Scene* scene, const GeneralParams& general_params, const SVTParams& svt_params, const AdaptiveConstructionParams& ac_params,
    std::uint64_t& total_samples, std::uint64_t& performed_samples, bool gather_stat_images, bool show_svd, std::mt19937& rng){
    float sample_time, recovery_time = 0.f, svd_time = 0.f, reconstruct_time = 0.f, adaptive_sample_time = 0.f, other_time = 0.f;
    if(general_params.adaptive_col){
        auto start_t = std::chrono::high_resolution_clock::now();
        computeUnoccludedSamples(slice, general_params.vsl, scene, vpls, general_params.min_dist);

        auto sample_t = std::chrono::high_resolution_clock::now();

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

        slice->singular_values.resize(8);
        slice->rank_ratio.resize(8);

        for(std::uint32_t i = 0; i < quadranted_vpls.size(); ++i){
            if(quadranted_vpls[i].size() == 0){
                continue;
            }

            std::vector<float> cluster_contribs(quadranted_vpls[i].size(), 0.f);

            //change this to an estimation to the 6d centroid
            for(std::uint32_t j = 0; j < slice->sample_indices.size(); ++j){
                for(std::uint32_t k = 0; k < quadranted_vpls[i].size(); ++k){
                    cluster_contribs[k] += slice->sample(j).unoccluded_samples[quadranted_vpl_indices[i][k]].getLuminance();
                }
            }

            std::uint32_t mat_rows = ac_params.vis_only ? slice->sample_indices.size() : slice->sample_indices.size() * 3;
            Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(mat_rows, quadranted_vpls[i].size());
            std::uint32_t max_rank = std::min(mat.rows(), mat.cols());
            std::uint32_t basis_rank;
            std::uint32_t samples;
            auto recover_start_t = std::chrono::high_resolution_clock::now();

            if(!ac_params.bin_vis){
                float svd_t, recon_t, adaptive_sample_t, other_t;
                std::tie(samples, svd_t, recon_t, adaptive_sample_t, other_t) = 
                    adaptiveMatrixReconstruction(mat, scene, slice, quadranted_vpls[i], general_params.min_dist, 
                    general_params.sample_perc, rng, ac_params.vis_only, ac_params.rec_trans, ac_params.import_sample, 
                    ac_params.force_resample, basis_rank, general_params.vsl, cluster_contribs, gather_stat_images, show_svd,
                    slice->singular_values[i]);
                svd_time += svd_t;
                reconstruct_time += recon_t;
                adaptive_sample_time += adaptive_sample_t;
                other_time += other_t;
            }
            else{
                std::vector<std::uint8_t> bin_vis(slice->sample_indices.size() * quadranted_vpls[i].size(), 0);
                if(true){
                    samples = adaptiveMatrixReconstructionBGE(bin_vis, scene, slice,
                        quadranted_vpls[i], general_params.min_dist, general_params.sample_perc, 
                        general_params.max_sample_perc, general_params.sample_inc, rng, 
                        basis_rank, cluster_contribs);
                }
                else{
                    samples = adaptiveMatrixReconstructionB(bin_vis, scene, slice,
                        quadranted_vpls[i], general_params.min_dist, general_params.sample_perc, 
                        general_params.max_sample_perc, general_params.sample_inc, rng, 
                        basis_rank, cluster_contribs);
                }
                

                for(auto col = 0; col < mat.cols(); ++col){
                    for(auto row = 0; row < mat.rows(); ++row){
                        mat(row, col) = bin_vis[col * slice->sample_indices.size() + row] == 0 ? -1.f : 1.f;
                    }
                }
            }

            auto recover_t = std::chrono::high_resolution_clock::now();
            total_performed_samples += samples;
            slice->rank_ratio[i] = float(basis_rank) / max_rank; 
            updateSliceWithMatData(mat, slice, ac_params.vis_only, ac_params.rec_trans, general_params.vsl, 
                scene, vpls, quadranted_vpl_indices[i], general_params.min_dist, gather_stat_images);

            recovery_time += std::chrono::duration_cast<std::chrono::duration<float>>(recover_t - recover_start_t).count();
        }

        float col_rank = 0.f;
        if(show_svd){
            Eigen::MatrixXf m(slice->sample_indices.size() * 3, vpls.size());
            for(std::uint32_t j = 0; j < vpls.size(); ++j){
                for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
                    bool los = sampleVisibility(scene, slice->sample(i).its, vpls[j], general_params.min_dist);
                    float coeff = los ? 1.f : 0.f;
                    Spectrum c = slice->sample(i).unoccluded_samples[j] * coeff;
                    float r, g, b;
                    c.toLinearRGB(r, g, b);
                    
                    m(i * 3, j) = r;
                    m(i * 3 + 1, j) = g;
                    m(i * 3 + 2, j) = b;
                }
            }

            auto svd = m.bdcSvd();
            std::uint32_t r;
            for(r = 0; r < svd.singularValues().size(); ++r){
                if(svd.singularValues()(r) < 1e-3f){
                    break;
                }
            }

            slice->singular_values2.resize(svd.singularValues().size());

            for(std::uint32_t i = 0; i < svd.singularValues().size(); ++i){
                slice->singular_values2[i] = svd.singularValues()(i);
            }

            col_rank = float(r) / std::min(m.rows(), m.cols());
        }

        slice->rank_ratio2 = col_rank;

        for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
            slice->sample(i).unoccluded_samples.clear();
            slice->sample(i).unoccluded_samples.shrink_to_fit();
        }

        slice->sample_ratio = float(total_performed_samples) / (slice->sample_indices.size() * vpls.size());

        sample_time = std::chrono::duration_cast<std::chrono::duration<float>>(sample_t - start_t).count();

        {
            std::lock_guard<std::mutex> lock(stats_mutex);
            performed_samples += total_performed_samples;
            total_samples += slice->sample_indices.size() * vpls.size();
        }
    }
    else{
        std::uint32_t num_samples = slice->sample_indices.size() * vpls.size() * general_params.sample_perc;
        auto start_t = std::chrono::high_resolution_clock::now();
        auto sample_t = start_t;

        if(ac_params.vis_only){
            Eigen::MatrixXf vmat = Eigen::MatrixXf::Zero(slice->sample_indices.size(), vpls.size());
            auto indices = calculateSparseSamples(scene, slice, vpls, vmat, num_samples,
                general_params.min_dist);
            sample_t = std::chrono::high_resolution_clock::now();
            float step_size = 1.9f;//(1.2f * lighting_matrix.rows() * lighting_matrix.cols()) / (indices.size() * 3.f); 

            std::vector<std::uint32_t> actual_vpl_indices(vpls.size());
            std::iota(actual_vpl_indices.begin(), actual_vpl_indices.end(), 0);

            Eigen::MatrixXf reconstructed;
            svt(reconstructed, vmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
            updateSliceWithMatData(reconstructed, slice, true, false, general_params.vsl, 
                scene, vpls, actual_vpl_indices, general_params.min_dist, gather_stat_images);
        }
        else{
            Eigen::MatrixXf rmat = Eigen::MatrixXf::Zero(slice->sample_indices.size(), vpls.size());
            Eigen::MatrixXf gmat = Eigen::MatrixXf::Zero(slice->sample_indices.size(), vpls.size());
            Eigen::MatrixXf bmat = Eigen::MatrixXf::Zero(slice->sample_indices.size(), vpls.size());

            auto indices = calculateSparseSamples(scene, slice, vpls, rmat, gmat, bmat, num_samples, 
                general_params.min_dist, general_params.vsl);
            sample_t = std::chrono::high_resolution_clock::now();
            float step_size = 1.9f;//(1.2f * lighting_matrix.rows() * lighting_matrix.cols()) / (indices.size() * 3.f); 

            Eigen::MatrixXf reconstructed_r, reconstructed_b, reconstructed_g;
            svt(reconstructed_r, rmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
            svt(reconstructed_g, gmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
            svt(reconstructed_b, bmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
            updateSliceWithMatData(reconstructed_r, reconstructed_g, reconstructed_b, slice);
        }

        auto recover_t = std::chrono::high_resolution_clock::now();

        sample_time = std::chrono::duration_cast<std::chrono::duration<float>>(sample_t - start_t).count();
        recovery_time = std::chrono::duration_cast<std::chrono::duration<float>>(recover_t - sample_t).count();

        {
            std::lock_guard<std::mutex> lock(stats_mutex);
            performed_samples += num_samples;
            total_samples += slice->sample_indices.size() * vpls.size();
        }
    }
    return std::make_tuple(sample_time, recovery_time, svd_time, reconstruct_time, adaptive_sample_time, other_time);
}

void sliceWorkerLS(std::vector<std::int32_t>& work, std::uint32_t thread_id, std::mutex& work_mutex, std::mutex& stats_mutex,
    const std::vector<KDTNode<ReconstructionSample>*>& slices, const std::vector<std::vector<std::uint32_t>>& cbsamp, 
    const Eigen::MatrixXf& contribs, const std::vector<std::vector<int>>& nn, const std::vector<VPL>& total_vpls, Scene* scene,
    const GeneralParams& general_params, const SVTParams& svt_params, const AdaptiveConstructionParams& ac_params,
    std::uint64_t& total_samples, std::uint64_t& performed_samples, float& clustering_time, float& recovery_time,
    float& sample_time, float& svd_time, float& reconstruct_time, float& adaptive_sample_time, 
    float& other_time, bool gather_stat_images, bool show_svd, std::uint32_t num_clusters){

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

        auto start_t = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<std::uint32_t>> cbsplit = clusterVPLsBySplitting(cbsamp, contribs, nn[slice_id], split_num_clusters, rng);
        auto vpls = sampleRepresentatives(contribs, total_vpls, cbsplit, rng, general_params.min_dist);

        auto cluster_t = std::chrono::high_resolution_clock::now();

        float sample_t, recover_t, svd_t, reconstruct_t, asample_t, other_t;
        std::tie(sample_t, recover_t, svd_t, reconstruct_t, asample_t, other_t) = 
            recover(slices[slice_id], stats_mutex, vpls, scene, general_params, svt_params, 
            ac_params, total_samples, performed_samples, gather_stat_images, show_svd, rng);

        {
            std::lock_guard<std::mutex> lock(work_mutex);
            work[thread_id] = -1;
        }

        {
            std::lock_guard<std::mutex> lock(stats_mutex);
            clustering_time += std::chrono::duration_cast<std::chrono::microseconds>(cluster_t - start_t).count() / 1000000.f;
            recovery_time += recover_t;
            sample_time += sample_t;
            svd_time += svd_t;
            reconstruct_time += reconstruct_t;
            adaptive_sample_time += asample_t;
            other_time += other_t;
        }
    }
}

void sliceWorkerMDLC(std::vector<std::int32_t>& work, std::uint32_t thread_id, std::mutex& work_mutex, std::mutex& stats_mutex,
    const std::vector<KDTNode<ReconstructionSample>*>& slices, LightTree* light_tree, const std::vector<VPL>& tot_vpls, Scene* scene,
    const GeneralParams& general_params, const SVTParams& svt_params, const AdaptiveConstructionParams& ac_params,
    std::uint64_t& total_samples, std::uint64_t& performed_samples, float& clustering_time, float& recovery_time, 
    float& sample_time, float& svd_time, float& reconstruct_time, float& adaptive_sample_time, 
    float& other_time, bool gather_stat_images, bool show_svd){

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

        auto start_t = std::chrono::high_resolution_clock::now();

        std::vector<Intersection> slice_points(slices[slice_id]->sample_indices.size());

        for(std::uint32_t i = 0; i < slice_points.size(); ++i){
            slice_points[i] = slices[slice_id]->sample(i).its;
        }

        auto vpls = light_tree->getClusteringForPoints(scene, slice_points);
        updateVPLRadii(vpls, general_params.min_dist);

        auto cluster_t = std::chrono::high_resolution_clock::now();

        float sample_t, recover_t, svd_t, reconstruct_t, asample_t, other_t;
        std::tie(sample_t, recover_t, svd_t, reconstruct_t, asample_t, other_t) = recover(slices[slice_id], stats_mutex, vpls, scene, general_params, svt_params, ac_params,
            total_samples, performed_samples, gather_stat_images, show_svd, rng);
        {
            std::lock_guard<std::mutex> lock(work_mutex);
            work[thread_id] = -1;
        }

        {
            std::lock_guard<std::mutex> lock(stats_mutex);
            clustering_time += std::chrono::duration_cast<std::chrono::duration<double>>(cluster_t - start_t).count();
            recovery_time += recover_t;
            sample_time += sample_t;
            svd_time += svd_t;
            reconstruct_time += reconstruct_t;
            adaptive_sample_time += asample_t;
            other_time += other_t;
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
    bool gather_stat_images, bool show_svd, float sample_perc, float max_sample_perc, float sample_inc,
    float min_dist, std::mt19937& rng, bool importance_sample, bool bin_vis){
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

    slice->singular_values.resize(8);
    slice->rank_ratio.resize(8);

    slice->visibility_coefficients.resize(slice->sample_indices.size() * vpls.size(), 0.f);



    for(std::uint32_t i = 0; i < quadranted_vpls.size(); ++i){
        if(quadranted_vpls[i].size() == 0){
            continue;
        }

        //add support for specular materials later. 
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
        std::uint32_t basis_rank;

        if(bin_vis){
            std::vector<std::uint8_t> bv(slice->sample_indices.size() * quadranted_vpls[i].size(), 0);
            if(true){
                samples = adaptiveMatrixReconstructionBGE(bv, scene, slice, 
                    quadranted_vpls[i], min_dist, sample_perc, max_sample_perc, sample_inc, rng, 
                    basis_rank, cluster_contribs);
            }
            else{
                samples = adaptiveMatrixReconstructionB(bv, scene, slice, 
                    quadranted_vpls[i], min_dist, sample_perc, max_sample_perc, sample_inc, rng, 
                    basis_rank, cluster_contribs);
            }

            total_performed_samples += samples;
            slice->rank_ratio[i] = float(basis_rank) / std::min(slice->sample_indices.size(), quadranted_vpls[i].size());

            Properties props("independent");
            ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
            sampler->configure();
            sampler->generate(Point2i(0));

            for(std::uint32_t j = 0; j < actual_vpl_index[i].size(); ++j){
                for(std::uint32_t k = 0; k < slice->sample_indices.size(); ++k){
                    std::uint32_t idx = actual_vpl_index[i][j] * slice->sample_indices.size() + k;
                    std::uint32_t bin_vis_idx = j * slice->sample_indices.size() + k;
                    slice->visibility_coefficients[idx] = bv[bin_vis_idx] == 0 ? 0.f : 1.f;

                    if(gather_stat_images){
                        std::uint8_t actual_visible = 
                            sampleVisibility(scene, slice->sample(k).its, vpls[actual_vpl_index[i][j]], min_dist) ?
                            1 : 0;
                        if(actual_visible != bv[bin_vis_idx]){
                            slice->sample(k).visibility_errors++;
                        }

                        if(slice->sample(k).its.isEmitter()){
                            Spectrum emitter_col = slice->sample(k).its.Le(-slice->sample(i).ray.d);
                            slice->sample(k).fully_sampled_color = emitter_col;
                        }
                        else{
                            std::uint32_t num_samples;

                            if(actual_visible){
                                slice->sample(k).fully_sampled_color += sample(scene, sampler, 
                                slice->sample(k).its, slice->sample(k).ray, vpls[actual_vpl_index[i][j]], min_dist, 
                                false, 5, false, slice->sample(k).intersected_scene, false, false, num_samples);
                            }
                        }
                    }
                }
            }
        }
        else{
            Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(slice->sample_indices.size(), quadranted_vpls[i].size());

            float svd_t, recon_t, adaptive_sample_t, other_t;
            std::tie(samples, svd_t, recon_t, adaptive_sample_t, other_t) = 
                adaptiveMatrixReconstruction(mat, scene, slice, quadranted_vpls[i], min_dist, 
                sample_perc, rng, true, false, importance_sample, false, basis_rank, false, cluster_contribs, 
                false, false, slice->singular_values[i]);
            total_performed_samples += samples;
            slice->rank_ratio[i] = float(basis_rank) / std::min(slice->sample_indices.size(), quadranted_vpls[i].size());

            for(std::uint32_t j = 0; j < actual_vpl_index[i].size(); ++j){
                for(std::uint32_t k = 0; k < slice->sample_indices.size(); ++k){
                    std::uint32_t idx = actual_vpl_index[i][j] * slice->sample_indices.size() + k;
                    slice->visibility_coefficients[idx] = (mat(k, j) + 1.f) / 2.f;
                }
            }
        }
        
    }

    slice->sample_ratio = float(total_performed_samples) / (slice->sample_indices.size() * vpls.size());

    return std::make_tuple(slice->sample_indices.size() * vpls.size(), total_performed_samples);
}

void writeVisibilityToFile(const std::vector<float>& coefficients, std::uint32_t num_lights, std::uint32_t num_receivers, std::string file_name){
    std::ofstream file;
    file.open(file_name);

    for(std::uint32_t i = 0; i < num_lights; ++i){
        for(std::uint32_t j = 0; j < num_receivers; ++j){
            std::uint32_t idx = i * num_receivers + j;
            file << coefficients[idx];
            if(j != num_receivers - 1){
                file << ", ";
            }
            else file << std::endl;
        }
    }
    file.close();
}

float getVisibilityRank(const std::vector<float>& coefficients, std::uint32_t num_lights, std::uint32_t num_receivers){
    std::vector<std::vector<std::uint8_t>> basis;

    for(std::uint32_t i = 0; i < num_lights; ++i){
        std::vector<std::uint8_t> curr_col(num_receivers);

        for(std::uint32_t j = 0; j < num_receivers; ++j){
            std::uint32_t idx = i * num_receivers + j;
            curr_col[j] = coefficients[idx] > 0.5f ? 1 : 0;
        }

        int matching_col = -1;
        for(std::uint32_t j = 0; j < basis.size(); ++j){
            bool match = true;
            for(std::uint32_t k = 0; k < num_receivers; ++k){
                if(basis[j][k] != curr_col[k]){
                    match = false;
                    break;
                }
            }
            if(match){
                matching_col = j;
                break;
            }
        }

        if(matching_col < 0){
            basis.push_back(curr_col);
        }
    }

    return float(basis.size()) / num_lights;
}

void clusterWorkerMDLC(BlockingQueue<HWWorkUnit>& input, BlockingQueue<HWWorkUnit>& output,
    std::vector<std::vector<VPL>>& vpls, const std::vector<VPL>& tot_vpls, LightTree* light_tree, Scene* scene, bool gather_stats, 
    bool show_svd, float sample_perc, float max_sample_perc, float sample_inc, float min_dist, std::uint32_t thread_id, std::uint32_t num_threads, 
    std::mutex& barrier_mutex, std::condition_variable& barrier, std::mutex& sample_update_mutex,
    std::uint64_t& total_samples, std::uint64_t& num_samples, bool bin_vis, bool importance_sample, bool show_vmat){
    
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
        std::tie(slice_samples, num_slice_sampled) = recoverHW(work_unit.first, vpls[work_unit.second], scene, gather_stats, show_svd, sample_perc,
            max_sample_perc, sample_inc, min_dist, rng, importance_sample, bin_vis);

        if(gather_stats){
            fs::path scene_path = scene->getDestinationFile();
            std::string filename_and_path = (scene_path.parent_path() / std::string("vmat/vmat_" + std::to_string(work_unit.first->sample(0).image_x) + "_" + std::to_string(work_unit.first->sample(0).image_y))).string();

            writeVisibilityToFile(work_unit.first->visibility_coefficients, vpls[work_unit.second].size(), work_unit.first->sample_indices.size(), 
                filename_and_path);
        
            work_unit.first->rank_ratio2 = getVisibilityRank(work_unit.first->visibility_coefficients, 
                vpls[work_unit.second].size(), work_unit.first->sample_indices.size());
        }
        
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
    const std::vector<VPL>& total_vpls, std::uint32_t num_clusters, Scene* scene, bool gather_stats, bool show_svd, 
    float sample_perc, float max_sample_perc, float sample_inc, float min_dist, std::uint32_t thread_id, std::uint32_t num_threads, std::mutex& barrier_mutex,
    std::condition_variable& barrier, std::mutex& sample_update_mutex, std::uint64_t& total_samples, std::uint64_t& num_samples,
    bool bin_vis, bool importance_sample, bool show_vmat){

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * thread_id);
    std::uint32_t split_num_clusters = std::min(total_vpls.size(), num_clusters - cbsamp.size());

    HWWorkUnit work_unit;
    while(input.pop(work_unit)){
        std::vector<std::vector<std::uint32_t>> cbsplit = 
            clusterVPLsBySplitting(cbsamp, contribs, nn[work_unit.second], split_num_clusters, rng);
        vpls[work_unit.second] = sampleRepresentatives(contribs, total_vpls, cbsplit, rng, min_dist);

        std::uint64_t slice_samples, num_slice_sampled;
        std::tie(slice_samples, num_slice_sampled) = recoverHW(work_unit.first, vpls[work_unit.second], scene, gather_stats, show_svd, sample_perc,
            max_sample_perc, sample_inc, min_dist, rng, importance_sample, bin_vis);
        if(gather_stats){
            fs::path scene_path = scene->getDestinationFile();
            std::string filename_and_path = (scene_path.parent_path() / std::string("vmat/vmat_" + std::to_string(work_unit.first->sample(0).image_x) + "_" + std::to_string(work_unit.first->sample(0).image_y))).string();
            writeVisibilityToFile(work_unit.first->visibility_coefficients, vpls[work_unit.second].size(), work_unit.first->sample_indices.size(), 
                filename_and_path);
            
            work_unit.first->rank_ratio2 = getVisibilityRank(work_unit.first->visibility_coefficients, 
                vpls[work_unit.second].size(), work_unit.first->sample_indices.size());
        }

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
    ac_params.vis_only = visibility_only_;
    ac_params.rec_trans = adaptive_recover_transpose_;
    ac_params.import_sample = adaptive_importance_sampling_;
    ac_params.force_resample = adaptive_force_resample_;
    ac_params.bin_vis = bin_vis_;
    
    SVTParams svt_params;
    svt_params.tolerance = tolerance_;
    svt_params.tau = tau_;
    svt_params.max_iter = max_iterations_;
    svt_params.trunc = truncated_;
    
    GeneralParams general_params;
    general_params.min_dist = min_dist_;
    general_params.adaptive_col = adaptive_col_sampling_;
    general_params.vsl = vsl_;
    general_params.sample_perc = sample_percentage_;
    general_params.max_sample_perc = max_sample_perc_;
    general_params.sample_inc = sample_inc_;

    std::uint32_t num_cores = 16;//std::min(slices.size(), (size_t)std::thread::hardware_concurrency());
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

    float total_recovery_time = 0.f;
    float total_clustering_time = 0.f;
    float total_sample_time = 0.f;
    float total_svd_time = 0.f;
    float total_reconstruct_time = 0.f;
    float total_adaptive_sample_time = 0.f;
    float total_other_time = 0.f;

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
            std::ref(vpls_), scene, general_params, svt_params, ac_params, std::ref(total_samples), 
            std::ref(amount_sampled), std::ref(total_clustering_time), std::ref(total_recovery_time), 
            std::ref(total_sample_time), std::ref(total_svd_time), std::ref(total_reconstruct_time), 
            std::ref(total_adaptive_sample_time), std::ref(total_other_time), gather_stat_images_, show_svd_,
            num_clusters_);
        }
    }
    else{
        ProgressReporter lt_pr("Constructing light tree", 1, job);
        light_tree = std::unique_ptr<LightTree>(new LightTree(vpls_, min_dist_, num_clusters_, 0.f, true));
        lt_pr.finish();
        std::cout << std::endl;

        for(std::uint32_t i = 0; i < num_cores; ++i){
            workers.emplace_back(sliceWorkerMDLC, std::ref(work), i, std::ref(work_mutex), std::ref(stats_mutex), std::ref(slices), 
            light_tree.get(), std::ref(vpls_), scene, general_params, svt_params, ac_params, std::ref(total_samples), 
            std::ref(amount_sampled), std::ref(total_clustering_time), std::ref(total_recovery_time),
            std::ref(total_sample_time), std::ref(total_svd_time), std::ref(total_reconstruct_time), 
            std::ref(total_adaptive_sample_time), std::ref(total_other_time), gather_stat_images_, show_svd_);
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

    std::cout << "Total clustering time: " << total_clustering_time << std::endl;
    std::cout << "Total sampling time: " << total_sample_time << std::endl;
    std::cout << "Total recovery time: " << total_recovery_time << std::endl;
    std::cout << "Total svd time: " << total_svd_time << std::endl;
    std::cout << "Total reconstruction time: " << total_reconstruct_time << std::endl;
    std::cout << "Total adaptive sampling time: " << total_adaptive_sample_time << std::endl;
    std::cout << "Total other time: " << total_other_time << std::endl;

    timings.push_back(total_clustering_time);
    timings.push_back(total_sample_time);
    timings.push_back(total_recovery_time);
    timings.push_back(total_svd_time);
    timings.push_back(total_reconstruct_time);
    timings.push_back(total_adaptive_sample_time);
    timings.push_back(total_other_time);

    if(adaptive_col_sampling_){
        float sample_perc = (float)amount_sampled / total_samples;
        timings.push_back(sample_perc);
        std::cout << "Sample percentage: " << sample_perc << std::endl;
    }

    return std::make_tuple(total_samples, amount_sampled);
}

std::tuple<std::uint64_t, std::uint64_t> MatrixReconstructionRenderer::renderHW(Scene* scene, std::uint32_t spp, const RenderJob *job, 
    std::vector<float>& timings, const std::vector<KDTNode<ReconstructionSample>*>& slices, 
    std::uint32_t samples_per_slice, std::uint32_t slice_size){
    auto start = std::chrono::high_resolution_clock::now();
    std::uint32_t num_workers = 15;//std::max(size_t(1), std::min(slices.size(), size_t(std::thread::hardware_concurrency() / 2)));
    std::uint32_t batch_size = 500 / (std::max(slice_size / 1000u, 1u) * std::max(1u, num_clusters_ / 4000));//std::max(num_workers * 2, 64u);

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
                num_clusters_, scene, gather_stat_images_, show_svd_, sample_percentage_, max_sample_perc_, sample_inc_, min_dist_, i,
                num_workers, std::ref(barrier_mutex), std::ref(barrier), std::ref(sample_update_mutex), std::ref(total_samples), std::ref(num_samples),
                bin_vis_, adaptive_importance_sampling_, show_vmat_);
        }
        
    }
    else{
        std::cout << "Creating light tree" << std::endl;
        light_tree = std::unique_ptr<LightTree>(new LightTree(vpls_, min_dist_, num_clusters_, 0.f, true));
        for(std::uint32_t i = 0; i < num_workers; ++i){
            clusterers.emplace_back(clusterWorkerMDLC, std::ref(to_cluster), std::ref(to_shade), std::ref(vpls), std::ref(vpls_),
                light_tree.get(), scene, gather_stat_images_, show_svd_, sample_percentage_, max_sample_perc_, sample_inc_, min_dist_, i,
                num_workers, std::ref(barrier_mutex), std::ref(barrier), std::ref(sample_update_mutex), 
                std::ref(total_samples), std::ref(num_samples), bin_vis_, adaptive_importance_sampling_, show_vmat_);
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

    std::uint32_t image_size = size.x * size.y;
    
    if(gather_stat_images_){
        std::vector<Spectrum> recovered(image_size, Spectrum(0.f));
        std::vector<Spectrum> full_sample(image_size, Spectrum(0.f));
        std::vector<std::vector<Spectrum>> rank;
        for(std::uint32_t i = 0; i < 8; ++i){
            rank.push_back(std::vector<Spectrum>(image_size, Spectrum(0.f)));
        }
        std::vector<Spectrum> samples(image_size, Spectrum(0.f));
        std::vector<Spectrum> slice_cols(image_size, Spectrum(0.f));
        std::vector<Spectrum> visibility_errors(image_size, Spectrum(0.f));
        std::uint64_t total_verr = 0;

        constructStatBuffers(recovered, full_sample, rank, samples, slice_cols, visibility_errors, total_verr, 
            slices, spp, size.x, num_clusters_);
        for(std::uint32_t i = 0; i < rank.size(); ++i){
            writeOutputImage(scene, "rank_" + std::to_string(i), size.x, size.y, true, rank[i]);
        }
        
        writeOutputImage(scene, "sample_ratio", size.x, size.y, true, samples);
        writeOutputImage(scene, "slices", size.x, size.y, true, slice_cols);
        writeOutputImage(scene, "full_sample", size.x, size.y, true, full_sample);
        writeOutputImage(scene, "visibility_errors", size.x, size.y, true, visibility_errors);

        if(hw_){
            std::vector<Spectrum> col_rank(image_size, Spectrum(0.f));
            constructSvdBuffers(col_rank, slices, spp, size.x);
            writeOutputImage(scene, "visibility_rank", size.x, size.y, true, col_rank);
        }
        
        /*float err = writeOutputErrorImage(scene, "error", size.x, size.y, true, recovered, full_sample, error_scale_);
        std::cout << "Total error: " << err << std::endl;
        timings.push_back(err);

        float vis_err_ratio = float(total_verr) / total_samples;
        timings.push_back(vis_err_ratio);*/
    }

    std::vector<float> sample_rate;
    float sr = (float)num_samples / total_samples;
    sample_rate.push_back(sr);
    std::cout << "Sample rate: " << sr << std::endl;

    writeOutputData(scene, "timings", false, timings, ',');
    writeOutputData(scene, "samplerates", false, sample_rate, ',');

    if(show_svd_){
        std::vector<Spectrum> col_rank(image_size, Spectrum(0.f));
        constructSvdBuffers(col_rank, slices, spp, size.x);

        for(std::uint32_t i = 0; i < slices.size(); ++i){
            for(std::uint32_t j = 0; j < slices[i]->singular_values.size(); ++j){
                writeOutputData(scene, "vis_sv_" + std::to_string(j), false, slices[i]->singular_values[j], ',');
            }
            
            writeOutputData(scene, "col_sv", false, slices[i]->singular_values2, ',');
        }

        writeOutputImage(scene, "coloured_rank", size.x, size.y, true, col_rank);
    }

    film->setBitmap(output_bitmap);

    return !cancel_;
}

MTS_NAMESPACE_END