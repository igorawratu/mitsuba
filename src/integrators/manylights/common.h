#ifndef COMMON_H
#define COMMON_H

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/vpl.h>
#include <cstdint>
#include <tuple>
#include <eigen3/Eigen/Dense>
#include <flann/flann.hpp>
#include <mutex>

#include "definitions.h"

MTS_NAMESPACE_BEGIN

float calculateError(Scene* scene, const std::vector<VPL>& vpls, float min_dist, std::uint8_t* image_buffer);

std::tuple<float, float, float> floatToRGB(float v);

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf> partialSvd(const Eigen::MatrixXf& mat, std::uint32_t num_singular_values);

Eigen::MatrixXf softThreshRank(const Eigen::MatrixXf& mat, float theta, const std::uint32_t initial, 
    const std::uint32_t step_size, std::mutex& block_mutex);

Spectrum sample(Scene* scene, Sampler* sampler, const Intersection& its, const VPL& vpl, float min_dist, bool check_occlusion);

Eigen::MatrixXf computeMoorePenroseInverse(const Eigen::MatrixXf& m);

const std::array<std::function<bool(float, const Intersection&)>, 6> searchers {
    [](float value, const Intersection& its){
        return value < its.p.x;
    },
    [](float value, const Intersection& its){
        return value < its.p.y;
    },
    [](float value, const Intersection& its){
        return value < its.p.z;
    },
    [](float value, const Intersection& its){
        return value < its.geoFrame.n.x;
    },
    [](float value, const Intersection& its){
        return value < its.geoFrame.n.y;
    },
    [](float value, const Intersection& its){
        return value < its.geoFrame.n.z;
    }
};

template<class Sample>
struct KDTNode{
    std::vector<std::uint32_t> sample_indices;
    std::vector<Sample>* samples;
    std::unique_ptr<KDTNode> left;
    std::unique_ptr<KDTNode> right;
    std::vector<std::vector<int>> nearest_neighbours;
    std::vector<std::vector<float>> neighbour_distances;

    KDTNode(std::vector<Sample>* sample_set) : sample_indices(), samples(sample_set), 
        left(nullptr), right(nullptr){

    }

    void Split(float norm_scale, std::uint32_t size_threshold){
        float maxf = std::numeric_limits<float>::max();
        float minf = std::numeric_limits<float>::min();
        Vector3f min_pos(maxf, maxf, maxf), max_pos(minf, minf, minf);
        Vector3f min_normal(maxf, maxf, maxf), max_normal(minf, minf, minf);

        for(size_t i = 0; i < sample_indices.size(); ++i){
            auto& curr_sample = (*samples)[sample_indices[i]];

            min_pos.x = std::min(curr_sample.its.p.x, min_pos.x);
            min_pos.y = std::min(curr_sample.its.p.y, min_pos.y);
            max_pos.z = std::max(curr_sample.its.p.z, max_pos.z);
            min_pos.z = std::min(curr_sample.its.p.z, min_pos.z);
            max_pos.x = std::max(curr_sample.its.p.x, max_pos.x);
            max_pos.y = std::max(curr_sample.its.p.y, max_pos.y);

            min_normal.x = std::min(curr_sample.its.geoFrame.n.x, min_normal.x);
            min_normal.y = std::min(curr_sample.its.geoFrame.n.y, min_normal.y);
            min_normal.z = std::min(curr_sample.its.geoFrame.n.z, min_normal.z);
            max_normal.x = std::max(curr_sample.its.geoFrame.n.x, max_normal.x);
            max_normal.y = std::max(curr_sample.its.geoFrame.n.y, max_normal.y);
            max_normal.z = std::max(curr_sample.its.geoFrame.n.z, max_normal.z);
        }

        std::array<std::pair<std::uint8_t, float>, 6> ranges;
        ranges[0] = std::make_pair(0, max_pos.x - min_pos.x);
        ranges[1] = std::make_pair(1, max_pos.y - min_pos.y);
        ranges[2] = std::make_pair(2, max_pos.z - min_pos.z);
        ranges[3] = std::make_pair(3, (max_normal.x - min_normal.x) * norm_scale);
        ranges[4] = std::make_pair(4, (max_normal.y - min_normal.y) * norm_scale);
        ranges[5] = std::make_pair(5, (max_normal.z - min_normal.z) * norm_scale);

        std::array<float, 6> midpoints;
        midpoints[0] = (max_pos.x + min_pos.x) / 2.f;
        midpoints[1] = (max_pos.y + min_pos.y) / 2.f;
        midpoints[2] = (max_pos.z + min_pos.z) / 2.f;
        midpoints[3] = (max_normal.x + min_normal.x) / 2.f;
        midpoints[4] = (max_normal.y + min_normal.y) / 2.f;
        midpoints[5] = (max_normal.z + min_normal.z) / 2.f;

        std::sort(ranges.begin(), ranges.end(), 
            [](const std::pair<std::uint8_t, float>& lhs, const std::pair<std::uint8_t, float>& rhs){
                return lhs.second > rhs.second;
            });

        left = std::unique_ptr<KDTNode>(new KDTNode(samples));
        right = std::unique_ptr<KDTNode>(new KDTNode(samples));

        for(std::uint32_t i = 0; i < sample_indices.size(); ++i){
            if(searchers[ranges[0].first](midpoints[ranges[0].first], sample(i).its)){
                right->sample_indices.push_back(sample_indices[i]);
            }
            else{
                left->sample_indices.push_back(sample_indices[i]);
            }
        }

        sample_indices.clear();
        nearest_neighbours.clear();
        neighbour_distances.clear();

        sample_indices.shrink_to_fit();
        nearest_neighbours.shrink_to_fit();
        neighbour_distances.shrink_to_fit();
    }

    Sample& sample(std::uint32_t index){
        assert(samples != nullptr && index < sample_indices.size());

        return (*samples)[sample_indices[index]];
    }

    void buildNearestNeighbours(std::uint32_t nn = 0){
        assert(nn < sample_indices.size());

        if(nn == 0){
            nn = std::min((size_t)100, sample_indices.size() / 10 + 1);
        }

        flann::Matrix<float> query_points(new float[sample_indices.size() * 3], sample_indices.size(), 3);
        constructFlannSampleMatrix(query_points, *(samples), sample_indices);
        flann::Index<flann::L2<float>> index(query_points, flann::KDTreeIndexParams(4));
        index.buildIndex();
        index.knnSearch(query_points, nearest_neighbours, neighbour_distances, nn, flann::SearchParams(128));

        for(std::uint32_t i = 0; i < nearest_neighbours.size(); ++i){
            auto iter = std::find(nearest_neighbours[i].begin(), nearest_neighbours[i].end(), i);
            if(iter != nearest_neighbours[i].end()){
                nearest_neighbours[i].pop_back();
                neighbour_distances[i].pop_back();
            }
            else{
                std::uint32_t pos = iter - nearest_neighbours[i].begin();
                nearest_neighbours[i].erase(iter);
                neighbour_distances[i].erase(neighbour_distances[i].begin() + pos);
            }
        }

        delete query_points.ptr();
    }

private:
    void constructFlannSampleMatrix(flann::Matrix<float>& mat, const std::vector<Sample>& samples, 
        const std::vector<std::uint32_t>& indices){

        assert(indices.size() > 0 && mat.rows == indices.size() && mat.cols == 3);

        for(std::uint32_t i = 0; i < indices.size(); ++i){
            float* curr_pos = (float*)mat[i];

            curr_pos[0] = samples[indices[i]].its.p.x;
            curr_pos[1] = samples[indices[i]].its.p.y;
            curr_pos[2] = samples[indices[i]].its.p.z;
        }
    }
};

template<class Sample>
void splitKDTree(KDTNode<Sample>* node, std::uint32_t size_threshold, float min_dist){
    if(node == nullptr || node->sample_indices.size() < size_threshold){
        return;
    }

    if(node->left == nullptr && node->right == nullptr){
        node->Split(min_dist / 10.f, size_threshold);
        splitKDTree(node->left.get(), size_threshold, min_dist);
        splitKDTree(node->right.get(), size_threshold, min_dist);
    }
}

template<class Sample>
void getSlices(KDTNode<Sample>* curr, std::vector<KDTNode<Sample>*>& slices){
    if(curr == nullptr){
        return;
    }

    if(curr->left == nullptr && curr->right == nullptr){
        slices.push_back(curr);
    }

    getSlices(curr->left.get(), slices);
    getSlices(curr->right.get(), slices);
}

MTS_NAMESPACE_END

#endif