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
#include "arpaca.hpp"
#include "RedSVD-h.hpp"

MTS_NAMESPACE_BEGIN

std::tuple<float, float, float> floatToRGB(float v);

template<typename MatrixType>
std::tuple<Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>, 
    Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic>> 
    partialSvd(const MatrixType& mat, typename MatrixType::Index num_singular_values){
    
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;

    Index max_rank = std::min(mat.rows(), mat.cols());
    num_singular_values = std::min(num_singular_values, max_rank);

    RedSVD::RedSVD<DenseMatrix> svd;
    svd.compute(mat, num_singular_values);

    DenseMatrix singular_values = DenseMatrix::Zero(mat.rows(), mat.cols());
    DenseMatrix u = DenseMatrix::Zero(mat.rows(), mat.rows());
    DenseMatrix v = DenseMatrix::Zero(mat.cols(), mat.cols());

    const DenseMatrix uvectors = svd.matrixU();
    const DenseMatrix vvectors = svd.matrixV();
    const ScalarVector svs = svd.singularValues();

    for(Index i = 0; i < svs.size(); ++i){
        singular_values(i, i) = svs(i);
        u.col(i) = uvectors.col(i);
        v.col(i) = vvectors.col(i);
    }

    return std::make_tuple(u, singular_values, v);
}

template<typename MatrixType>
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic> 
    softThreshRankNoTrunc(const MatrixType& mat, typename MatrixType::Scalar theta){
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;

    auto svd = mat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto singular_values = svd.singularValues();
    DenseMatrix diagonal_singular = DenseMatrix::Zero(svd.matrixU().rows(), svd.matrixV().rows());
    
    for(Index j = 0; j < singular_values.rows(); ++j){
        diagonal_singular(j, j) = std::max(Scalar(0), singular_values(j, 0) - theta);
    }

    return svd.matrixU() * diagonal_singular * svd.matrixV().transpose();
}

template<typename MatrixType>
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic> 
    softThreshRank(const MatrixType& mat, float theta, const typename MatrixType::Index initial, 
    const typename MatrixType::Index step_size){
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ScalarVector;

    Index max_rank = std::min(mat.rows(), mat.cols());
    Index curr_step_size = std::min(max_rank, initial);

    RedSVD::RedSVD<DenseMatrix> svd;

    while(true){
        Index actual_dim_to_compute = std::min(max_rank, curr_step_size + max_rank / 10);
        svd.compute(mat, actual_dim_to_compute);
        if(curr_step_size == max_rank || svd.singularValues()(curr_step_size - 1) < theta){
            break;
        }

        curr_step_size = std::min(max_rank, curr_step_size + step_size);
    }

    DenseMatrix singular_values = DenseMatrix::Zero(mat.rows(), mat.cols());
    DenseMatrix u = DenseMatrix::Zero(mat.rows(), mat.rows());
    DenseMatrix v = DenseMatrix::Zero(mat.cols(), mat.cols());

    const DenseMatrix uvectors = svd.matrixU();
    const DenseMatrix vvectors = svd.matrixV();
    const ScalarVector svs = svd.singularValues();

    for(Index i = 0; i < curr_step_size; ++i){
        singular_values(i, i) = std::max((Scalar)0, svs(i) - theta);
        u.col(i) = uvectors.col(i);
        v.col(i) = vvectors.col(i);
    }

    return u * singular_values * v.transpose();
}

template<typename MatrixType>
Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, Eigen::Dynamic> 
    computeMoorePenroseInverse(const MatrixType& m){
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> DenseMatrix;

    auto svd = m.jacobiSvd(Eigen::ComputeFullV | Eigen::ComputeFullU);
    auto singular_values = svd.singularValues();
    DenseMatrix svmat = Eigen::MatrixXf::Zero(m.cols(), m.rows());
    for(Index i = 0; i < singular_values.size(); ++i){
        if(fabs(singular_values(i)) > 1e-10){
            svmat(i, i) = 1.0 / singular_values(i);
        }
    }
    return svd.matrixV() * svmat * svd.matrixU().adjoint();
}


Spectrum sample(Scene* scene, Sampler* sampler, Intersection& its, const Ray& ray, const VPL& vpl, 
    float min_dist, bool check_occlusion, std::uint32_t max_specular_bounces, 
    bool perform_ray_intersection, bool& intersected, bool show_emitter, bool vsl, std::uint32_t& samples_taken);

bool sampleVisibility(Scene* scene, const Intersection& its, const VPL& vpl, float min_dist);

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
    std::vector<float> rank_ratio;
    float sample_ratio;
    float rank_ratio2;
    std::vector<std::vector<float>> singular_values;
    std::vector<float> singular_values2;
    std::vector<float> visibility_coefficients;

    KDTNode(std::vector<Sample>* sample_set) : sample_indices(), samples(sample_set), 
        left(nullptr), right(nullptr), rank_ratio(0.f), sample_ratio(0.f){

    }

    void Split(float norm_scale, std::uint32_t min_size, std::uint32_t max_size){
        float maxf = std::numeric_limits<float>::max();
        Vector3f min_pos(maxf, maxf, maxf), max_pos(-maxf, -maxf, -maxf);
        Vector3f min_normal(maxf, maxf, maxf), max_normal(-maxf, -maxf, -maxf);

        std::array<float, 6> means;
        for(std::uint8_t i = 0; i < 6; ++i){
            means[i] = 0.f;
        }

        for(size_t i = 0; i < sample_indices.size(); ++i){
            auto& curr_sample = (*samples)[sample_indices[i]];

            means[0] += curr_sample.its.p.x;
            means[1] += curr_sample.its.p.y;
            means[2] += curr_sample.its.p.z;

            means[3] += curr_sample.its.geoFrame.n.x;
            means[4] += curr_sample.its.geoFrame.n.y;
            means[5] += curr_sample.its.geoFrame.n.z;

            min_pos.x = std::min(curr_sample.its.p.x, min_pos.x);
            min_pos.y = std::min(curr_sample.its.p.y, min_pos.y);
            min_pos.z = std::min(curr_sample.its.p.z, min_pos.z);
            max_pos.x = std::max(curr_sample.its.p.x, max_pos.x);
            max_pos.y = std::max(curr_sample.its.p.y, max_pos.y);
            max_pos.z = std::max(curr_sample.its.p.z, max_pos.z);

            min_normal.x = std::min(curr_sample.its.geoFrame.n.x, min_normal.x);
            min_normal.y = std::min(curr_sample.its.geoFrame.n.y, min_normal.y);
            min_normal.z = std::min(curr_sample.its.geoFrame.n.z, min_normal.z);
            max_normal.x = std::max(curr_sample.its.geoFrame.n.x, max_normal.x);
            max_normal.y = std::max(curr_sample.its.geoFrame.n.y, max_normal.y);
            max_normal.z = std::max(curr_sample.its.geoFrame.n.z, max_normal.z);
        }
        for(std::uint8_t i = 0; i < 6; ++i){
            means[i] /= sample_indices.size();
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

        if(left->sample_indices.size() < min_size || right->sample_indices.size() < min_size){
            if(sample_indices.size() > max_size * 2){
                std::uint32_t midpoint = sample_indices.size() / 2;
                left->sample_indices.clear();
                left->sample_indices.insert(left->sample_indices.end(), sample_indices.begin(), sample_indices.begin() + midpoint);
                right->sample_indices.clear();
                right->sample_indices.insert(right->sample_indices.end(), sample_indices.begin() + midpoint, sample_indices.end());
            }
            else{
                left.release();
                right.release();
                return;
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
void splitKDTree(KDTNode<Sample>* node, std::uint32_t size_threshold, std::uint32_t min_slice_size,
    float min_dist){
    if(node == nullptr || node->sample_indices.size() < size_threshold){
        return;
    }

    if(node->left == nullptr && node->right == nullptr){
        node->Split(min_dist / 10.f, min_slice_size, size_threshold);

        if(node->left != nullptr && node->right != nullptr){
            splitKDTree(node->left.get(), size_threshold, min_slice_size, min_dist);
            splitKDTree(node->right.get(), size_threshold, min_slice_size, min_dist);
        }
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

std::uint64_t upperPo2(std::uint64_t v);

MTS_NAMESPACE_END

#endif