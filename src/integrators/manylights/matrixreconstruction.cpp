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
#include <unordered_map>
#include <thread>
#include "filewriter.h"

#include "common.h"

MTS_NAMESPACE_BEGIN

void updateVPLRadii(std::vector<VPL>& vpls, float min_dist){
	//radius calculation, 11 to account 10 + 1 for adding a node's self in nearest neighbours
	std::uint32_t num_neighbours = std::min((std::uint32_t)vpls.size(), 11u);

	std::uint32_t num_sl = 0;
	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(vpls[i].type == ESurfaceVPL){
			num_sl++;
		}
	} 

    flann::Matrix<float> dataset(new float[num_sl * 3], num_sl, 3);
	std::uint32_t curr_light = 0;
	std::vector<std::uint32_t> pointlight_indices;
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(vpls[i].type == ESurfaceVPL){
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
    float min_dist, std::vector<std::uint32_t>& sampled_slice_rows, std::uint32_t samples_per_slice,
    bool vsl){

    assert(slices.size() > 0 && vpls.size() > 0);

    Eigen::MatrixXf contributions = Eigen::MatrixXf::Zero(slices.size() * 3 * samples_per_slice, vpls.size());

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        std::set<std::uint32_t> generated;
        while(generated.size() < samples_per_slice){
            generated.insert(rand() % slices[i]->sample_indices.size());
        }

        std::uint32_t pos = 0;
        for(auto idx = generated.begin(); idx != generated.end(); ++idx){
            sampled_slice_rows[i * samples_per_slice + pos] = *idx;
            auto& curr_sample = slices[i]->sample(*idx); 
            for(std::uint32_t j = 0; j < vpls.size(); ++j){
                std::uint32_t num_samples;
                Spectrum c = sample(scene, sampler, curr_sample.its, curr_sample.ray, vpls[j], min_dist, true, 5, false, 
                    curr_sample.intersected_scene, false, vsl, num_samples);

                float r, g, b;
                c.toLinearRGB(r, g, b);

                contributions(i * 3 * samples_per_slice + pos, j) = r;
                contributions(i * 3 * samples_per_slice + pos + 1, j) = g;
                contributions(i * 3 * samples_per_slice + pos + 2, j) = b;
            }
            pos++;
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
            if(min > dist_mat(j, i)){
                min = dist_mat(j, i);
                min_idx = j;
            }
        }

        clusters[min_idx].push_back(nonzero_norm_cols[i]);
    }

    if(zero_norm_cols.size() > 0){
        clusters.push_back(zero_norm_cols);
    }

    return clusters;
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
    
    /*float total = 0.f;
    for(std::uint32_t i = 0; i < cluster.size(); ++i){
        for(std::uint32_t j = i; j < cluster.size(); ++j){
            total += (contrib_mat.col(cluster[j]) - contrib_mat.col(cluster[i])).norm() * norm_cache[cluster[i]] * 
                norm_cache[cluster[j]];
        }
    }

    return total;*/

    /*float total_norm = 0.f;
    Eigen::MatrixXf variance = Eigen::MatrixXf::Zero(contrib_mat.rows(), 1);
    Eigen::MatrixXf mean = Eigen::MatrixXf::Zero(contrib_mat.rows(), 1);

    for(std::uint32_t i = 0; i < cluster.size(); ++i){
        total_norm += norm_cache[cluster[i]];
        Eigen::MatrixXf old_mean = mean;
        mean = mean + (contrib_mat.col(cluster[i]) - mean) / float(i + 1);
        variance = variance + 
            (contrib_mat.col(cluster[i]) - mean).cwiseProduct((contrib_mat.col(cluster[i]) - old_mean));
    }
    variance /= float(cluster.size() - 1);

    return variance.norm() * total_norm;*/

}

std::vector<std::vector<std::uint32_t>> clusterVPLsBySplitting(const std::vector<std::vector<std::uint32_t>>& sampling_clusters, 
    const Eigen::MatrixXf& cluster_contributions, const std::vector<int>& nn, std::uint32_t num_clusters,
     std::mt19937& rng, std::vector<float>& norm_cache){
    
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

    norm_cache.reserve(cluster_contributions.cols());
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

        /*float midpoint = (max_d + min_d) / 2.f;

        for(std::uint32_t i = 0; i < vpl_projection_distances.size(); ++i){
            if(vpl_projection_distances[i].second < midpoint){
                new_cluster1.push_back(vpl_projection_distances[i].first);
            }
            else{
                new_cluster2.push_back(vpl_projection_distances[i].first);
            }
        }*/

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
        
        float tot_cluster_power = 0.f;
        float total_vpl_power = 0.f;
        for(std::uint32_t j = 0; j < clusters[i].size(); ++j){
            tot_cluster_power += contributions.col(clusters[i][j]).norm();
            float r, g, b;
            vpls[clusters[i][j]].P.toLinearRGB(r, g, b);
            total_vpl_power += sqrt(r * r + g * g + b * b);
        }

        std::uniform_real_distribution<float> gen(0, tot_cluster_power);
        std::uint32_t rep_idx = 0;
        float v = 0.f;
        float selection = gen(rng);
        float selected_norm = 0.f;

        for(std::uint32_t j = 0; j < clusters[i].size(); ++j){
            /*selected_norm = contributions.col(clusters[i][j]).norm();
            v += selected_norm;
            rep_idx = j;

            if(v > selection){
                break;
            }*/

            float curr_norm = contributions.col(clusters[i][j]).norm();
            if(curr_norm > selected_norm){
                selected_norm = curr_norm;
                rep_idx = j;
            }
        }

        representatives.push_back(vpls[clusters[i][rep_idx]]);
        float r, g, b;
        representatives.back().P.toLinearRGB(r, g, b);
        float rep_power = sqrt(r * r + g * g + b * b);
        representatives.back().P = representatives.back().P * total_vpl_power / rep_power;
    }

    //updateVPLRadii(representatives, min_dist);

    return representatives;
}

std::unique_ptr<KDTNode<ReconstructionSample>> constructKDTree(Scene* scene, std::uint32_t size_threshold, 
    std::vector<ReconstructionSample>& samples, float min_dist, std::uint32_t min_slice_size, std::uint32_t spp){

    auto kdt_root = std::unique_ptr<KDTNode<ReconstructionSample>>(new KDTNode<ReconstructionSample>(&samples));

    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    samples.resize(film->getSize().y * film->getSize().x * spp);

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for (std::int32_t y = 0; y < film->getSize().y; ++y) {
        #pragma omp parallel for
        for (std::int32_t x = 0; x < film->getSize().x; ++x) {
            std::uint32_t cell_dim = sqrt(spp) + 0.5f;
            float cell_side_len = 1.f / cell_dim;

            for(std::uint32_t i = 0; i < spp; ++i){
                //calculates the position the ray intersects with
                Ray ray;

                float x_jitter = 0.5f;//sampler->next1D() / cell_dim;
                float y_jitter = 0.5f;//sampler->next1D() / cell_dim;
                float x_off = (i % spp) * cell_side_len;
                float y_off = (i / spp) * cell_side_len;

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

                    if(curr_sample.its.getBSDF()->getType() & BSDF::ESmooth || curr_sample.its.isEmitter()){
                        break;
                    }

                    if(++num_bounces > 5){
                        break;
                    }

                    BSDFSamplingRecord bsdf_sample_record(curr_sample.its, sampler);
                    //bsdf_sample_record.typeMask = BSDF::EReflection;
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

    splitKDTree(kdt_root.get(), size_threshold, min_slice_size, min_dist);

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
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
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
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
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
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));
    
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        auto& curr_sample = slice->sample(i);
        curr_sample.unoccluded_samples.resize(vpls.size());

        for(std::uint32_t j = 0; j < vpls.size(); ++j){
            std::uint32_t samples_for_cell;
            curr_sample.unoccluded_samples[j] = sample(scene, sampler, curr_sample.its, curr_sample.ray,
                vpls[j], min_dist, false, 5, false, curr_sample.intersected_scene, false, vsl, samples_for_cell);
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
    float min_dist, bool fully_compute_occlusion){
    
    Sampler *sampler = nullptr;

    Properties props("independent");
    sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
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
                for(std::uint32_t j = 0; j < vpls.size(); ++j){
                    float coeff = ((recover_transpose ? mat(j, i) : mat(i, j)) + 1.f) / 2.f;
                    slice->sample(i).color += slice->sample(i).unoccluded_samples[j] * coeff;
                    if(fully_compute_occlusion){
                        bool los = sampleVisibility(scene, slice->sample(i).its, vpls[j], min_dist);
                        float correct_coeff = los ? 1.f : 0.f;
                        slice->sample(i).fully_sampled_color += slice->sample(i).unoccluded_samples[j] * correct_coeff;
                    }
                }
            }
        }
        else{
            for(std::uint32_t j = 0; j < vpls.size(); ++j){
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

void constructStatBuffers(std::vector<Spectrum>& recovered, std::vector<Spectrum>& full, std::vector<Spectrum>& ranks,
    std::vector<Spectrum>& samples, std::vector<Spectrum>& slice_cols, 
    const std::vector<KDTNode<ReconstructionSample>*>& slices, std::uint32_t spp, std::uint32_t image_width){
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
            
            std::tie(r, g, b) = floatToRGB(slice->rank_ratio);
            Spectrum rank_col;
            rank_col.fromLinearRGB(r, g, b);
            ranks[index] += rank_col / spp;

            std::tie(r, g, b) = floatToRGB(slice->sample_ratio);
            Spectrum sample_ratio_col;
            sample_ratio_col.fromLinearRGB(r, g, b);
            samples[index] += sample_ratio_col / spp;

            slice_cols[index] += slice_col / spp;
        }
    }
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

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

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

std::uint32_t adaptiveMatrixReconstruction(Eigen::MatrixXf& mat, Scene* scene, KDTNode<ReconstructionSample>* slice, 
    const std::vector<VPL>& vpls, float min_dist, float sample_perc,
    std::mt19937& rng, bool visibility_only, bool recover_transpose, bool importance_sample, bool force_resample,
    std::uint32_t& basis_rank, bool vsl, const std::vector<float>& col_estimations){

    assert(sample_perc > 0.f && slice->sample_indices.size() > 0 && vpls.size() > 0);
    std::random_shuffle(slice->sample_indices.begin(), slice->sample_indices.end());

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
    Eigen::MatrixXf q_normalized;

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

                /*q_normalized = q;

                for(std::uint32_t j = 0; j < q_omega.cols(); ++j){
                    float l2 = q_omega.col(j).norm();
                    if(l2 < std::numeric_limits<float>::epsilon()){
                        q_omega.col(j).setZero();
                        q_normalized.col(j).setZero();
                    }
                    else{
                        q_omega.col(j).normalize();

                        for(std::uint32_t k = j + 1; k < q_omega.cols(); ++k){
                            q_omega.col(k) = q_omega.col(k) - q_omega.col(j).dot(q_omega.col(k)) * q_omega.col(j);
                        }
                        q_normalized.col(j) /= l2;
                    }
                }

                q_omega.transposeInPlace();*/

                auto svd = q_omega.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
                auto sv = svd.singularValues();
                Eigen::MatrixXf singular_val_inv = Eigen::MatrixXf::Zero(sv.size(), sv.size());

                for(std::uint32_t j = 0; j < sv.size(); ++j){
                    singular_val_inv(j, j) = sv(j) < 1e-5f ? 0.f : 1.f / sv(j);
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
            //reconstructed = sample_omega.rows() == reconstructed.rows() ? sample_omega : q_normalized * q_omega * sample_omega;            

            float d = 0;
            for(std::uint32_t j = 0; j < sampled.size(); ++j){
                d += std::abs(reconstructed(sampled[j], 0) - sample_omega(j, 0));
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
            if(d > 1e-3f){
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
                q_normalized = q = reconstructed;
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
    bool truncated, bool show_slices, bool vsl, bool gather_stat_images) : 
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
};

void sliceWorker(std::vector<std::int32_t>& work, std::uint32_t thread_id, std::mutex& work_mutex, std::mutex& stats_mutex,
    const std::vector<KDTNode<ReconstructionSample>*>& slices, const std::vector<std::vector<std::uint32_t>>& cbsamp, 
    const Eigen::MatrixXf& contribs, const std::vector<std::vector<int>>& nn, const std::vector<VPL>& total_vpls, Scene* scene,
    const GeneralParams& general_params, const SVTParams& svt_params, const AdaptiveConstructionParams& ac_params,
    std::uint32_t& total_samples, std::uint32_t& performed_samples, bool gather_stat_images){

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count() * thread_id);
    std::uint32_t split_num_clusters = std::min(total_vpls.size(), 1000 - cbsamp.size());

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

        std::vector<float> norm_cache;
        std::vector<std::vector<std::uint32_t>> cbsplit = clusterVPLsBySplitting(cbsamp, contribs, nn[slice_id], split_num_clusters, rng, norm_cache);

        auto vpls = sampleRepresentatives(contribs, total_vpls, cbsplit, rng, general_params.min_dist);

        if(general_params.adaptive_col){
            std::uint32_t total_slice_samples = 
                computeUnoccludedSamples(slices[slice_id], general_params.vsl, scene, vpls, general_params.min_dist);

            std::uint32_t base = slices[slice_id]->sample_indices.size() * vpls.size();

            //std::cout << total_slice_samples << " " << base << " " <<
            //    (float)total_slice_samples / base << std::endl;

            std::vector<float> cluster_contribs(vpls.size(), 0.f);
            for(std::uint32_t j = 0; j < slices[slice_id]->sample_indices.size(); ++j){
                for(std::uint32_t k = 0; k < vpls.size(); ++k){
                    cluster_contribs[k] += slices[slice_id]->sample(j).unoccluded_samples[k].getLuminance();
                }
            }

            std::uint32_t mat_rows = ac_params.vis_only ? slices[slice_id]->sample_indices.size() : slices[slice_id]->sample_indices.size() * 3;
            Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(mat_rows, vpls.size());
            std::uint32_t max_rank = std::min(mat.rows(), mat.cols());
            std::uint32_t basis_rank;
            std::uint32_t samples = adaptiveMatrixReconstruction(mat, scene, slices[slice_id], vpls, general_params.min_dist, 
                general_params.sample_perc, rng, ac_params.vis_only, ac_params.rec_trans, ac_params.import_sample, 
                ac_params.force_resample, basis_rank, general_params.vsl, cluster_contribs);
            
            slices[slice_id]->sample_ratio = float(samples) / (mat.rows() * mat.cols());
            slices[slice_id]->rank_ratio = float(basis_rank) / max_rank;

            updateSliceWithMatData(mat, slices[slice_id], ac_params.vis_only, ac_params.rec_trans, general_params.vsl, 
                scene, vpls, general_params.min_dist, gather_stat_images);

            {
                std::lock_guard<std::mutex> lock(stats_mutex);
                performed_samples += samples;
                total_samples += slices[slice_id]->sample_indices.size() * vpls.size();
            }
        }
        else{
            std::uint32_t num_samples = slices[slice_id]->sample_indices.size() * vpls.size() * general_params.sample_perc;

            if(ac_params.vis_only){
                Eigen::MatrixXf vmat = Eigen::MatrixXf::Zero(slices[slice_id]->sample_indices.size(), vpls.size());
                auto indices = calculateSparseSamples(scene, slices[slice_id], vpls, vmat, num_samples,
                    general_params.min_dist);

                float step_size = 1.9f;//(1.2f * lighting_matrix.rows() * lighting_matrix.cols()) / (indices.size() * 3.f); 

                Eigen::MatrixXf reconstructed;
                svt(reconstructed, vmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
                updateSliceWithMatData(reconstructed, slices[slice_id], true, false, general_params.vsl, 
                    scene, vpls, general_params.min_dist, gather_stat_images);
            }
            else{
                Eigen::MatrixXf rmat = Eigen::MatrixXf::Zero(slices[slice_id]->sample_indices.size(), vpls.size());
                Eigen::MatrixXf gmat = Eigen::MatrixXf::Zero(slices[slice_id]->sample_indices.size(), vpls.size());
                Eigen::MatrixXf bmat = Eigen::MatrixXf::Zero(slices[slice_id]->sample_indices.size(), vpls.size());

                auto indices = calculateSparseSamples(scene, slices[slice_id], vpls, rmat, gmat, bmat, num_samples, 
                    general_params.min_dist, general_params.vsl);

                float step_size = 1.9f;//(1.2f * lighting_matrix.rows() * lighting_matrix.cols()) / (indices.size() * 3.f); 

                Eigen::MatrixXf reconstructed_r, reconstructed_b, reconstructed_g;
                svt(reconstructed_r, rmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
                svt(reconstructed_g, gmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
                svt(reconstructed_b, bmat, step_size, svt_params.tolerance, svt_params.tau, svt_params.max_iter, indices, svt_params.trunc);
                updateSliceWithMatData(reconstructed_r, reconstructed_g, reconstructed_b, slices[slice_id]);
            }

            {
                std::lock_guard<std::mutex> lock(stats_mutex);
                performed_samples += num_samples;
                total_samples += slices[slice_id]->sample_indices.size() * vpls.size();
            }
        }

        {
            std::lock_guard<std::mutex> lock(work_mutex);
            work[thread_id] = -1;
        }
    }
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

    std::uint32_t samples_per_slice = 5;

    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, size);
    std::uint8_t* output_image = output_bitmap->getUInt8Data();
    memset(output_image, 0, output_bitmap->getBytesPerPixel() * size.x * size.y);

    std::cout << "constructing kd tree" << std::endl;
    auto kdt_root = constructKDTree(scene, slice_size_, samples_, min_dist_, samples_per_slice, spp);
    std::cout << "getting slices" << std::endl;

    std::vector<KDTNode<ReconstructionSample>*> slices;
    getSlices(kdt_root.get(), slices);

    std::uint32_t amount_sampled = 0;
    std::uint32_t total_samples = 0;

    std::cout << "calculating contributions" << std::endl;
    std::uint32_t num_clusters = std::min(std::uint32_t(vpls_.size()), std::uint32_t(1000 * 0.3f));
    std::vector<std::uint32_t> sampled_slice_rows(slices.size() * samples_per_slice);
    Eigen::MatrixXf cluster_contributions = calculateClusterContributions(vpls_, slices, scene, min_dist_,
        sampled_slice_rows, samples_per_slice, vsl_);
    std::cout << "cluster by sampling" << std::endl;
    std::vector<std::vector<std::uint32_t>> clusters_by_sampling = clusterVPLsBySampling(cluster_contributions,
        num_clusters, vpls_, min_dist_);

    std::cout << "approx nearest neighbour computation" << std::endl;
    //approx nearest neighbours here
    flann::Matrix<float> data_points(new float[sampled_slice_rows.size() * 3], 
        sampled_slice_rows.size(), 3);
    flann::Matrix<float> slice_centroids(new float[slices.size() * 3], 
        slices.size(), 3);

    for(std::uint32_t i = 0; i < sampled_slice_rows.size(); ++i){
        float* curr_col = (float*)data_points[i];

        std::uint32_t slice_idx = i / samples_per_slice;
        auto& sample = slices[slice_idx]->sample(sampled_slice_rows[i]);

        curr_col[0] = sample.its.p.x;
        curr_col[1] = sample.its.p.y;
        curr_col[2] = sample.its.p.z;
    }

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

    flann::Index<flann::L2<float>> index(data_points, flann::KDTreeIndexParams(4));
    index.buildIndex();

    std::vector<std::vector<int>> nearest_neighbours;
    std::vector<std::vector<float>> neighbour_distances;

    index.knnSearch(slice_centroids, nearest_neighbours, neighbour_distances, 10 * samples_per_slice, flann::SearchParams(128));

    AdaptiveConstructionParams ac_params;
    ac_params.vis_only = visibility_only_;
    ac_params.rec_trans = adaptive_recover_transpose_;
    ac_params.import_sample = adaptive_importance_sampling_;
    ac_params.force_resample = adaptive_force_resample_;
    
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

    std::uint32_t num_cores = std::min(slices.size(), (size_t)std::thread::hardware_concurrency());
    std::uint32_t scheduled_slices = num_cores;

    std::vector<std::int32_t> work(num_cores);
    std::iota(work.begin(), work.end(), 0);
    std::mutex work_mutex, stats_mutex;
    std::vector<std::thread> workers;

    for(std::uint32_t i = 0; i < num_cores; ++i){
        workers.emplace_back(sliceWorker, std::ref(work), i, std::ref(work_mutex), std::ref(stats_mutex), std::ref(slices), 
        std::ref(clusters_by_sampling), std::ref(cluster_contributions), std::ref(nearest_neighbours), 
        std::ref(vpls_), scene, general_params, svt_params, ac_params, std::ref(total_samples), 
        std::ref(amount_sampled), gather_stat_images_);
    }

    std::uint32_t finished_slices = 0;
    
    std::cout << "reconstructing slice 0/" << slices.size();
    fflush(stdout);
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
            printf("\rreconstructing slice %d/%d     ", ++finished_slices, int(slices.size()));
            fflush(stdout);

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
    
    copySamplesToBuffer(output_image, samples_, size, spp);

    if(gather_stat_images_){
        std::uint32_t image_size = size.x * size.y;
        std::vector<Spectrum> recovered(image_size, Spectrum(0.f));
        std::vector<Spectrum> full_sample(image_size, Spectrum(0.f));
        std::vector<Spectrum> rank(image_size, Spectrum(0.f));
        std::vector<Spectrum> samples(image_size, Spectrum(0.f));
        std::vector<Spectrum> slice_cols(image_size, Spectrum(0.f));

        constructStatBuffers(recovered, full_sample, rank, samples, slice_cols, slices, spp, size.x);
        writeOutputImage(scene, "rank", size.x, size.y, true, rank);
        writeOutputImage(scene, "sample_ratio", size.x, size.y, true, samples);
        writeOutputImage(scene, "slices", size.x, size.y, true, slice_cols);
        writeOutputErrorImage(scene, "error", size.x, size.y, true, recovered, full_sample, 1.f);
    }

    if(adaptive_col_sampling_){
        float sample_perc = (float)amount_sampled / total_samples;
        std::cout << "Sample percentage: " << sample_perc << std::endl;
    }
    
    film->setBitmap(output_bitmap);

    return !cancel_;
}

MTS_NAMESPACE_END