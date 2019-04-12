#include "rowcolumnsampling.h"

#include <mitsuba/core/plugin.h>

#include <random>
#include <queue>
#include <utility>
#include <chrono>
#include <algorithm>

#include "common.h"

MTS_NAMESPACE_BEGIN

std::vector<std::pair<std::uint32_t, std::uint32_t>> subsampleRows(std::uint32_t rows, std::tuple<std::uint32_t, std::uint32_t> resolution){
    std::vector<std::pair<std::uint32_t, std::uint32_t>> row_indices;

    //find factors to try get a decent bucket size for the requested number of buckets.   
    typedef std::pair<std::int32_t, std::int32_t> FactorPair;

    auto comparator = [](FactorPair l, FactorPair r){
        return abs(l.first - l.second) > abs(r.first - r.second);
    };

    std::priority_queue<FactorPair, std::vector<FactorPair>, decltype(comparator)> factors(comparator);

    for(int i = 1; i < (int)(sqrt((float)rows) + 1.f); ++i){
        if(rows % i == 0){
            factors.push(std::make_pair(i, rows / i));
        }
    }

    auto dims = factors.top();
    bool rows_more = std::get<0>(resolution) > std::get<1>(resolution);
    std::uint32_t bucket_rows = rows_more ? std::max(dims.first, dims.second) : std::min(dims.first, dims.second);
    std::uint32_t bucket_cols = rows_more ? std::min(dims.first, dims.second) : std::max(dims.first, dims.second);

    std::uint32_t bucket_width = std::get<0>(resolution) / bucket_rows;
    std::uint32_t bucket_height = std::get<1>(resolution) / bucket_cols;

    bucket_rows += std::get<0>(resolution) % bucket_rows == 0 ? 0 : 1;
    bucket_cols += std::get<1>(resolution) % bucket_cols == 0 ? 0 : 1;

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);
    
    for(std::uint32_t i = 0; i < bucket_rows; ++i){
        for(std::uint32_t j = 0; j < bucket_cols; ++j){
            std::uint32_t x_start = j * bucket_width;
            std::uint32_t y_start = i * bucket_height;

            std::uniform_int_distribution<std::uint32_t> gen_y(0, 
                std::min(bucket_height, std::get<1>(resolution) - y_start) - 1);

            std::uniform_int_distribution<std::uint32_t> gen_x(0, 
                std::min(bucket_width, std::get<0>(resolution) - x_start) - 1);

            int x = gen_x(rng) + x_start;
            int y = gen_y(rng) + y_start;

            row_indices.push_back(std::make_pair(x, y));
        }
    }

    return row_indices;
}

const double PI = 3.14159265359;

std::vector<float> calculateLightContributions(const std::vector<VPL>& vpls, 
                                    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& rows,
                                    const Scene* scene, float min_dist){
    std::vector<float> contributions(vpls.size(), 0.f);

    const Sensor* sensor = scene->getSensor();

    //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
    Point2 aperture_sample(0.5f, 0.5f);
    Float time_sample(0.5f);

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    for(size_t j = 0; j < rows.size(); ++j){
        Point2 sample_position(rows[j].first + 0.5f, rows[j].second + 0.5f);

        Ray ray;
        sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

        Intersection its;
        bool intersected;

        for(size_t i = 0; i < vpls.size(); ++i){
            Spectrum s = sample(const_cast<Scene*>(scene), sampler, its, ray, vpls[i], min_dist, true, 
                    10, i == 0, intersected, true, true);

            if(!intersected || its.isEmitter()){
                break;
            }

            contributions[i] += s.getLuminance();
        }
    }

    return contributions;
}

std::vector<VPL> calculateClustering(std::vector<VPL> vpls, std::vector<float> contributions, 
    std::uint32_t num_clusters, float min_dist){

    std::vector<VPL> output_clusters;

    if(num_clusters == 0){
        return output_clusters;
    }

    if(vpls.size() <= num_clusters){
        return vpls;
    }
    
    std::uint32_t clusters_by_sampling = num_clusters > 1 ? num_clusters * 2 / 3 : num_clusters;

    auto comparator = [](const std::pair<float, std::uint32_t>& l, const std::pair<float, std::uint32_t>& r){
        return l.first > r.first;
    };

    typedef std::priority_queue<std::pair<float, std::uint32_t>, std::vector<std::pair<float, 
        std::uint32_t>>, decltype(comparator)> DistancePQueue;

    std::vector<std::vector<VPL>> clusters;
    std::vector<std::pair<float, DistancePQueue*>> distance_sqr_from_clusters;

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    //cluster by sampling
    std::uniform_int_distribution<std::uint32_t> first_pick(0, vpls.size() - 1);

    std::uint32_t pick_idx = first_pick(rng);
    clusters.emplace_back(1, vpls[pick_idx]);
    clusters.back()[0].emitterScale = contributions[pick_idx];

    vpls.erase(vpls.begin() + pick_idx);
    contributions.erase(contributions.begin() + pick_idx);

    //this is necessary because we can't put the actual priority queues with custom comparators in a vector
    //we want to remove from due to copying assignment for lambdas being deleted
    //NB!!!! this vector must never be changed after creation so pointers remain the same
    std::vector<DistancePQueue> actual_priority_queues(vpls.size(), DistancePQueue(comparator));

    for(size_t i = 0; i < vpls.size(); ++i){
        float dsqr = (clusters.back()[0].its.p - vpls[i].its.p).lengthSquared() + 
            dot(clusters.back()[0].its.shFrame.n, vpls[i].its.shFrame.n) * min_dist;
        distance_sqr_from_clusters.emplace_back(0.f, &actual_priority_queues[i]);
        distance_sqr_from_clusters.back().first = dsqr;
        distance_sqr_from_clusters.back().second->push(std::make_pair(dsqr, 0));
    }

    while(clusters.size() < clusters_by_sampling){
        float total_dist = 0.f;
        for(size_t i = 0; i < distance_sqr_from_clusters.size(); ++i){
            total_dist += distance_sqr_from_clusters[i].first;
        }

        std::uniform_real_distribution<float> gen_pick(0.f, total_dist);
        float pick = gen_pick(rng);
        pick_idx = distance_sqr_from_clusters.size() - 1;

        for(size_t i = 0; i < distance_sqr_from_clusters.size(); ++i){
            pick -= distance_sqr_from_clusters[i].first;
            if(pick <= 0){
                pick_idx = i;
                break;
            }
        }

        clusters.emplace_back(1, vpls[pick_idx]);
        clusters.back()[0].emitterScale = contributions[pick_idx];

        vpls.erase(vpls.begin() + pick_idx);
        contributions.erase(contributions.begin() + pick_idx);
        distance_sqr_from_clusters.erase(distance_sqr_from_clusters.begin() + pick_idx);

        for(size_t i = 0; i < distance_sqr_from_clusters.size(); ++i){
            float dsqr = (clusters.back()[0].its.p - vpls[i].its.p).lengthSquared() + 
                dot(clusters.back()[0].its.shFrame.n, vpls[i].its.shFrame.n) * min_dist;
            distance_sqr_from_clusters[i].first += dsqr;
            distance_sqr_from_clusters[i].second->push(std::make_pair(dsqr, clusters.size() - 1));
        }
    }

    for(size_t i = 0; i < distance_sqr_from_clusters.size(); ++i){
        std::uint32_t cluster_idx = distance_sqr_from_clusters[i].second->top().second;
        clusters[cluster_idx].push_back(vpls[i]);
        clusters[cluster_idx].back().emitterScale = contributions[i];
    }

    //cluster by splitting
    typedef std::pair<std::vector<VPL>, float> ClusterContributionPair;

    //it's not possible to split size 1 clusters, so we put them at the back of the queue regardless of their contribution
    auto cluster_comparator = [](const ClusterContributionPair& l, const ClusterContributionPair& r){
        if(l.first.size() == 1){
            return true;
        }

        if(r.first.size() == 1){
            return false;
        }

        return l.second < r.second;
    };

    std::priority_queue<ClusterContributionPair, std::vector<ClusterContributionPair>, 
        decltype(cluster_comparator)> split_queue(cluster_comparator);

    for(size_t i = 0; i < clusters.size(); ++i){
        float total_contrib = 0.f;
        for(size_t j = 0; j < clusters[i].size(); ++j){
            total_contrib += clusters[i][j].emitterScale;
        }

        //it's ok to move this since we won't be using it any longer in the vector
        split_queue.push(std::make_pair(std::move(clusters[i]), total_contrib));
    }

    std::uniform_real_distribution<float> gen(0, 1.f);

    auto sort_comparator = [](const std::pair<VPL, float>& l, const std::pair<VPL, float>& r){
        return l.second < r.second;
    };

    while(split_queue.size() < num_clusters){
        ClusterContributionPair largest_cluster(split_queue.top());

        if(largest_cluster.first.size() == 1){
            std::cerr << "I should not be here" << std::endl;
            exit(0);
        }

        split_queue.pop();

        Vector3f random_line(gen(rng), gen(rng), gen(rng));
        Vector3f random_line2(gen(rng), gen(rng), gen(rng));

        //important to discard vectors with a larger norm than 1 as otherwise there will be higher densities
        //in the corners
        while(sqrt(random_line.lengthSquared() + random_line2.lengthSquared()) > 1.f || 
            sqrt(random_line.lengthSquared() + random_line2.lengthSquared()) == 0.f){
            random_line.x = gen(rng);
            random_line.y = gen(rng);
            random_line.z = gen(rng);
            random_line2.x = gen(rng);
            random_line2.y = gen(rng);
            random_line2.z = gen(rng);
        }

        std::vector<std::pair<VPL, float>> ordered_projections;

        //project points in cluster to 1D
        for(size_t i = 0; i < largest_cluster.first.size(); ++i){
            Vector3f vectorized_p(largest_cluster.first[i].its.p[0],largest_cluster.first[i].its.p[1], largest_cluster.first[i].its.p[2]);
            float projection_dist = (dot(random_line, vectorized_p) + dot(largest_cluster.first[i].its.geoFrame.n, random_line2)) /
                (dot(random_line, random_line) + dot(random_line2, random_line2));
            ordered_projections.push_back(std::make_pair(largest_cluster.first[i], projection_dist));
        }

        std::sort(ordered_projections.begin(), ordered_projections.end(), sort_comparator);

        float contrib_half = largest_cluster.second / 2.f;
        
        std::uint8_t idx = 0;
        ClusterContributionPair split_clusters[2];

        for(size_t i = 0; i < ordered_projections.size(); ++i){
            if(idx == 0 && (contrib_half <= 0.f || i == ordered_projections.size() - 2)){
                idx = 1;
            }

            split_clusters[idx].first.push_back(ordered_projections[i].first);
            split_clusters[idx].second += ordered_projections[i].first.emitterScale;

            contrib_half -= ordered_projections[i].first.emitterScale;
        }

        for(std::uint8_t i = 0; i < 2; ++i){
            split_queue.push(std::move(split_clusters[i]));
        }
    }

    //representative picking and emission scaling
    while(split_queue.size() > 0){
        const ClusterContributionPair& cluster = split_queue.top();

        float total_contrib = 0.f;
        float total_lum = 0.f;
        for(size_t i = 0; i < cluster.first.size(); ++i){
            total_contrib += cluster.first[i].emitterScale;
            total_lum += cluster.first[i].P.getLuminance();
        }

        std::uniform_real_distribution<float> gen_representative(0.f, total_contrib);
        float representative = gen_representative(rng);

        std::uint32_t representative_idx = cluster.first.size() - 1;

        for(size_t i = 0; i < cluster.first.size(); ++i){
            representative -= cluster.first[i].emitterScale;
            if(representative <= 0.f){
                representative_idx = i;
                break;
            }
        }

        output_clusters.push_back(cluster.first[representative_idx]);
        output_clusters.back().emitterScale = 1.f;
        output_clusters.back().P = output_clusters.back().P / output_clusters.back().P.getLuminance() * total_lum;

        split_queue.pop();
    }

    return output_clusters;
}

RowColumnSampling::RowColumnSampling(const std::vector<VPL>& vpls, std::uint32_t rows, std::uint32_t cols,
        std::tuple<std::uint32_t, std::uint32_t> resolution, const Scene* scene, float min_dist) : min_dist_(min_dist){

    auto indices = subsampleRows(rows, resolution);
    auto contributions = calculateLightContributions(vpls, indices, scene, min_dist_);
    clustering_ = calculateClustering(vpls, contributions, cols, min_dist);
}

RowColumnSampling::RowColumnSampling(const RowColumnSampling& other) : 
    clustering_(other.clustering_), min_dist_(other.min_dist_){
}

RowColumnSampling::RowColumnSampling(RowColumnSampling&& other) : 
    clustering_(std::move(other.clustering_)), min_dist_(other.min_dist_){

}

RowColumnSampling& RowColumnSampling::operator = (const RowColumnSampling& other){
    if(&other != this){
        clustering_ = other.clustering_;
        min_dist_ = other.min_dist_;
    }

    return *this;
}

RowColumnSampling& RowColumnSampling::operator = (RowColumnSampling&& other){
    if(&other != this){
        clustering_ = std::move(other.clustering_);
        min_dist_ = other.min_dist_;
    }

    return *this;
}

RowColumnSampling::~RowColumnSampling(){

}

std::vector<VPL> RowColumnSampling::getClusteringForPoint(const Intersection& its){
    return clustering_;
}

MTS_NAMESPACE_END