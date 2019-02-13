#include "matrixseparationrenderer.h"

#include <vector>
#include <limits>
#include <utility>
#include <memory>
#include <random>
#include <tuple>
#include <cmath>
#include <chrono>
#include <iostream>
#include <functional>
#include <numeric>

#include <mitsuba/core/plugin.h>

#include <eigen3/Eigen/Dense>
#include <flann/flann.hpp>

#include "definitions.h"
#include "common.h"

#include "arpaca.hpp"

MTS_NAMESPACE_BEGIN

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
                
void constructFlannSampleMatrix(flann::Matrix<float>& mat, const std::vector<RowSample>& samples, 
    const std::vector<std::uint32_t>& indices){

    assert(indices.size() > 0 && mat.rows == indices.size() && mat.cols == 3);

    for(std::uint32_t i = 0; i < indices.size(); ++i){
        float* curr_pos = (float*)mat[i];

        curr_pos[0] = samples[indices[i]].its.p.x;
        curr_pos[1] = samples[indices[i]].its.p.y;
        curr_pos[2] = samples[indices[i]].its.p.z;
    }
}

struct KDTNode{
    std::vector<std::uint32_t> sample_indices;
    std::vector<RowSample>* samples;
    std::unique_ptr<KDTNode> left;
    std::unique_ptr<KDTNode> right;
    std::vector<std::vector<int>> nearest_neighbours;
    std::vector<std::vector<float>> neighbour_distances;

    KDTNode(std::vector<RowSample>* sample_set) : sample_indices(), samples(sample_set), 
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

    RowSample& sample(std::uint32_t index){
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
};

void splitKDTree(KDTNode* node, std::uint32_t size_threshold, float min_dist){
    if(node == nullptr || node->sample_indices.size() < size_threshold){
        return;
    }

    if(node->left == nullptr && node->right == nullptr){
        node->Split(min_dist / 10.f, size_threshold);
        splitKDTree(node->left.get(), size_threshold, min_dist);
        splitKDTree(node->right.get(), size_threshold, min_dist);
    }
}

void getSlices(KDTNode* curr, std::vector<KDTNode*>& slices){
    if(curr == nullptr){
        return;
    }

    if(curr->left == nullptr && curr->right == nullptr){
        slices.push_back(curr);
    }

    getSlices(curr->left.get(), slices);
    getSlices(curr->right.get(), slices);
}

std::unique_ptr<KDTNode> constructKDTree(Scene* scene, const std::vector<VPL>& vpls, 
    std::uint32_t size_threshold, float min_dist, std::vector<RowSample>& samples){

    auto kdt_root = std::unique_ptr<KDTNode>(new KDTNode(&samples));

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

            RowSample curr_sample(x, y, vpls.size(), intersected, its);;

            if(its.isEmitter()){
                curr_sample.emitter_color = its.Le(-ray.d);
            }

            for (std::uint32_t i = 0; i < vpls.size(); ++i) {
                curr_sample.col_samples[i] = sample(scene, sampler, its, vpls[i], min_dist, false);
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

std::vector<VISIBILITY> knnPredictor(KDTNode* slice, std::uint32_t neighbours, std::uint32_t col, std::mt19937& rng){
    assert(slice != nullptr && neighbours > 0 && slice->sample_indices.size() > 0
        && col < slice->sample(0).visibility.size());

    std::uint32_t num_neighbours = std::min(neighbours, (std::uint32_t)slice->nearest_neighbours.size());
    std::vector<VISIBILITY> output(slice->sample_indices.size());

    if(num_neighbours == 0){
        for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
            if(slice->sample(slice->sample_indices[i]).visibility[col] != VISIBLE || 
                slice->sample(slice->sample_indices[i]).visibility[col] != NOT_VISIBLE){
                output[i] = P_VISIBLE;
            }
            else output[i] = slice->sample(slice->sample_indices[i]).visibility[col];
        }

        return output;
    }

    std::uniform_real_distribution<float> gen(0.f, 1.f);

    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        if(slice->sample(i).visibility[col] == VISIBLE || slice->sample(i).visibility[col] == NOT_VISIBLE){
            output[i] = slice->sample(i).visibility[col];
        }
        else{
            std::vector<std::uint32_t> sampled_neighbours_indices;
            for(std::uint32_t j = 0; j < slice->nearest_neighbours[i].size(); ++j){
                if(slice->sample(slice->nearest_neighbours[i][j]).visibility[col] == VISIBLE || 
                    slice->sample(slice->nearest_neighbours[i][j]).visibility[col] == NOT_VISIBLE){
                    sampled_neighbours_indices.push_back(slice->nearest_neighbours[i][j]);
                    if(sampled_neighbours_indices.size() == num_neighbours){
                        break;
                    }
                }
            }

            if(sampled_neighbours_indices.size() > 0){
                std::uint32_t idx = std::min(sampled_neighbours_indices.size() - 1, 
                    (size_t)(gen(rng) * sampled_neighbours_indices.size()));
                output[i] = slice->sample(sampled_neighbours_indices[idx]).visibility[col] == NOT_VISIBLE ? 
                    P_NOT_VISIBLE : P_VISIBLE;
            }
            else output[i] = gen(rng) < 0.5f ? P_NOT_VISIBLE : P_VISIBLE;
        }
    }

    return output;
}

std::vector<VISIBILITY> linearPredictor(KDTNode* slice, std::uint32_t col, float min_dist, std::mt19937& rng){
    assert(slice != nullptr && slice->sample_indices.size() > 0);

    Eigen::MatrixXf q(slice->sample_indices.size(), 4);

    std::vector<std::uint32_t> sampled;
    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        if(slice->sample(i).visibility[col] == VISIBLE || slice->sample(i).visibility[col] == NOT_VISIBLE){
            sampled.push_back(i);
        }
        q(i, 0) = slice->sample(i).its.p.x;
        q(i, 1) = slice->sample(i).its.p.y;
        q(i, 2) = slice->sample(i).its.p.z;
        q(i, 3) = 1;
    }

    Eigen::MatrixXf x(sampled.size(), 4);
    Eigen::MatrixXf y(sampled.size(), 1);

    for(std::uint32_t i = 0; i < sampled.size(); ++i){
        x(i, 0) = slice->sample(sampled[i]).its.p.x;
        x(i, 1) = slice->sample(sampled[i]).its.p.y;
        x(i, 2) = slice->sample(sampled[i]).its.p.z;
        x(i, 3) = 1;

        y(i, 0) = slice->sample(sampled[i]).visibility[col] == VISIBLE ? 1 : -1;
    }

    //might be able to speed this up by calculating pseudoinverse via svd method instead
    Eigen::MatrixXf w = (x.transpose() * x).inverse() * x.transpose() * y;

    Eigen::MatrixXf predictions = q * w;

    std::vector<VISIBILITY> output(slice->sample_indices.size());

    std::uniform_real_distribution<float> gen(0.f, 1.f);

    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        if(slice->sample(i).visibility[col] == VISIBLE || slice->sample(i).visibility[col] == NOT_VISIBLE){
            output[i] = slice->sample(i).visibility[col];
        }
        else{
            float pv0 = std::min(1.f, std::max(0.f, 0.5f - predictions(i, 0) * 0.5f));
            output[i] = gen(rng) < pv0 ? P_NOT_VISIBLE : P_VISIBLE;
        }
    }

    return output;
}

std::vector<VISIBILITY> naiveBayes(KDTNode* slice, std::uint32_t col, std::uint32_t neighbours, 
    const std::vector<VPL>& vpls, std::mt19937& rng, const std::vector<std::vector<int>>& nearest_vpls){

    assert(slice != nullptr && slice->sample_indices.size() > 0 && slice->sample(0).visibility.size() > 0 && 
        vpls.size() > col);
    
    std::uint32_t num_neighbours = std::min(neighbours, (std::uint32_t)nearest_vpls.size());
    std::vector<VISIBILITY> output(slice->sample_indices.size());
    if(num_neighbours == 0){
        std::fill(output.begin(), output.end(), P_VISIBLE);

        return output;
    }

    float tot_vis = 0.f;
    float neighbor_vis = 0.f;
    float tot_sampled = 0.f;

    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        for(std::uint32_t j = 0; j < num_neighbours; ++j){
            VISIBILITY v = slice->sample(i).visibility[nearest_vpls[col][j]];
            if(v == NOT_VISIBLE || v == VISIBLE){
                if(v == VISIBLE){
                    tot_vis += 1.f;
                    neighbor_vis += 1.f;
                }
                tot_sampled += 1.f;
            }
        }

        VISIBILITY v = slice->sample(i).visibility[col];
        if(v == NOT_VISIBLE || v == VISIBLE){
            if(v == VISIBLE){
                tot_vis += 1.f;
            }
            tot_sampled += 1.f;
        }
    }

    float pv = tot_vis / tot_sampled;
    float pj_v = tot_vis == 0.f ? 1.f : (tot_vis - neighbor_vis) / tot_vis;
    float pij = 1.f / (slice->sample_indices.size() * (num_neighbours + 1));

    std::uniform_real_distribution<float> gen(0.f, 1.f);

    for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
        VISIBILITY v = slice->sample(i).visibility[col];
        if(v == VISIBLE || v == NOT_VISIBLE){
            output[i] = v;
            continue;
        }
        
        float pi_v = 0.f;

        if(tot_vis == 0.f){
            pi_v = 1.f;
        }
        else{
            for(std::uint32_t k = 0; k < slice->nearest_neighbours[i].size(); ++k){
                for(std::uint32_t j = 0; j < num_neighbours; ++j){
                    VISIBILITY v = slice->sample(slice->nearest_neighbours[i][k]).visibility[nearest_vpls[col][j]];
                    if(v == VISIBLE){
                        pi_v += 1.f;
                    }   
                }
            }
            pi_v /= (float)slice->nearest_neighbours[i].size();
            pi_v /= tot_vis;
        }
        
        float pv_ij = (pi_v * pj_v * pv) / pij;

        //std::cout << pv_ij << " " << pv << " " << pi_v << " " << pj_v << std::endl;

        output[i] = gen(rng) < pv_ij ? P_VISIBLE : P_NOT_VISIBLE;
    }

    return output;
}

void performInitialVisibilitySamples(KDTNode* slice, float sample_percentage, std::mt19937& rng, 
    Scene* scene, const std::vector<VPL>& vpls){

    assert(slice != nullptr);

    std::uint32_t num_samples = sample_percentage * slice->sample_indices.size();
    
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        std::vector<std::uint32_t> sample_set(slice->sample_indices.size());
        for(std::uint32_t j = 0; j < slice->sample_indices.size(); ++j){
            sample_set[j] = j;
        }

        std::vector<std::uint32_t> to_sample;

        while(to_sample.size() < num_samples){
            std::uniform_int_distribution<std::uint32_t> gen(0, sample_set.size() - 1);
            auto pos = gen(rng);
            to_sample.push_back(sample_set[pos]);
            sample_set[pos] = sample_set.back();
            sample_set.pop_back();
        }

        for(std::uint32_t j = 0; j < to_sample.size(); ++j){
            Point ray_origin = slice->sample(to_sample[j]).its.p;
            Ray shadow_ray(ray_origin, normalize(vpls[i].its.p - ray_origin), 0.f);

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            slice->sample(to_sample[j]).visibility[i] = VISIBLE;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if(abs((ray_origin - vpls[i].its.p).length() - t) > 0.0001f ){
                    slice->sample(to_sample[j]).visibility[i] = NOT_VISIBLE;
                }
            }
        }
    }
}

std::set<std::uint32_t> sampleAndPredictVisibility(KDTNode* slice, float sample_percentage, std::mt19937& rng, Scene* scene,
    const std::vector<VPL>& vpls, const std::set<std::uint32_t>& active_columns, float error_threshold, float min_dist,
    std::uint32_t prediction_mask, const std::vector<std::vector<int>>& nearest_vpls){

    std::set<std::uint32_t> new_active_cols;

    std::uint32_t num_samples = sample_percentage * slice->sample_indices.size();
    std::vector<std::vector<VISIBILITY>> new_predictions(active_columns.size());

    for(auto iter = active_columns.begin(); iter != active_columns.end(); iter++){
        std::vector<std::uint32_t> unsampled;
        for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
            auto& curr_samples = slice->sample(i);
            if(curr_samples.visibility[*iter] != VISIBLE && curr_samples.visibility[*iter] != NOT_VISIBLE){
                unsampled.push_back(i);
            }
        }

        std::vector<std::vector<VISIBILITY>> predictions;
        
        if(prediction_mask & 1){
            predictions.push_back(linearPredictor(slice, *iter, min_dist, rng));
        }

        if(prediction_mask & 2){
            predictions.push_back(naiveBayes(slice, *iter, 3, vpls, rng, nearest_vpls));
        }
        
        if(prediction_mask & 4){
            predictions.push_back(knnPredictor(slice, 3, *iter, rng));
        }

        std::vector<std::uint32_t> to_sample;
        
        if(num_samples < unsampled.size()){
            while(to_sample.size() < num_samples){
                std::uniform_int_distribution<std::uint32_t> gen(0, unsampled.size() - 1);
                auto pos = gen(rng);
                to_sample.push_back(unsampled[pos]);
                unsampled[pos] = unsampled.back();
                unsampled.pop_back();
            }
        }
        else{
            to_sample = unsampled;
        }

        std::vector<VISIBILITY> samples(to_sample.size());

        for(std::uint32_t i = 0; i < to_sample.size(); ++i){
            Point ray_origin = slice->sample(to_sample[i]).its.p;
            Ray shadow_ray(ray_origin, normalize(vpls[*iter].its.p - ray_origin), 0.f);

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            samples[i] = VISIBLE;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if(abs((ray_origin - vpls[*iter].its.p).length() - t) > 0.0001f ){
                    samples[i] = NOT_VISIBLE;
                }
            }
        }

        std::uint32_t smallest_idx = 0;
        std::uint32_t smallest_error = std::numeric_limits<std::uint32_t>::max();

        for(std::uint32_t i = 0; i < predictions.size(); ++i){
            std::uint32_t curr_error = 0;
            for(std::uint32_t j = 0; j < to_sample.size(); ++j){
                if((samples[j] == VISIBLE && predictions[i][to_sample[j]] != P_VISIBLE) || 
                    (samples[j] == NOT_VISIBLE && predictions[i][to_sample[j]] != P_NOT_VISIBLE)){
                    curr_error++;
                }

                predictions[i][to_sample[j]] = samples[j];
            }

            if(curr_error < smallest_error){
                smallest_idx = i;
                smallest_error = curr_error;
            }
        }

        for(std::uint32_t i = 0; i < slice->sample_indices.size(); ++i){
            slice->sample(i).visibility[*iter] = predictions[smallest_idx][i];
            slice->sample(i).predictors[*iter] = smallest_idx;
        }

        float normalized_error = (float)smallest_error / to_sample.size();
        //std::cout << normalized_error << std::endl;
        if(normalized_error > error_threshold){
            new_active_cols.insert(*iter);
        }
    }

    return new_active_cols;
}

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> sliceToMatrices(KDTNode* node){
    assert(node->sample_indices.size() > 0 && node->sample(0).col_samples.size() > 0);

    Eigen::MatrixXf slice_mat(node->sample_indices.size(), node->sample(0).col_samples.size() * 3);
    Eigen::MatrixXf error_mat(node->sample_indices.size(), node->sample(0).col_samples.size() * 3);

    for(std::uint32_t row = 0; row < node->sample_indices.size(); ++row){
        for(std::uint32_t col = 0; col < node->sample(row).col_samples.size(); ++col){
            float r, g , b;
            node->sample(row).col_samples[col].toLinearRGB(r, g, b);
            if(node->sample(row).visibility[col] == VISIBLE || node->sample(row).visibility[col] == P_VISIBLE){
                slice_mat(row, col * 3) = r;
                slice_mat(row, col * 3 + 1) = g;
                slice_mat(row, col * 3 + 2) = b;
            }
            else{
                slice_mat(row, col * 3) = slice_mat(row, col * 3 + 1) = slice_mat(row, col * 3 + 2) = 0.f;
            }
            
            error_mat(row, col * 3) = r > 1e-4 ? 2.f / r : 0.f;
            error_mat(row, col * 3 + 1) = g > 1e-4 ? 2.f / g : 0.f;
            error_mat(row, col * 3 + 2) = b > 1e-4 ? 2.f / b : 0.f;
        }
    }

    return std::make_tuple(slice_mat, error_mat);
}

void matrixToSlice(KDTNode* node, const Eigen::MatrixXf& mat){
    assert(node->sample_indices.size() == (std::uint32_t)mat.rows() && 
        (node->sample(0).col_samples.size() * 3) == (std::uint32_t)mat.cols());

    for(std::uint32_t row = 0; row < node->sample_indices.size(); ++row){
        for(std::uint32_t col = 0; col < node->sample(row).col_samples.size(); ++col){
            float r = mat(row, col * 3);
            float g = mat(row, col * 3 + 1);
            float b = mat(row, col * 3 + 2);

            node->sample(row).col_samples[col].fromLinearRGB(r, g, b);
        }
    }
}

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> separate(const Eigen::MatrixXf& mat, const Eigen::MatrixXf& two_q,
    const std::uint32_t max_iterations, float beta, float step, float theta, const Eigen::MatrixXf& sampled,
    bool show_rank){

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::normal_distribution<float> dist(0.f, 0.1f);

    Eigen::MatrixXf y = Eigen::MatrixXf::Identity(mat.cols(), mat.cols());

    Eigen::MatrixXf x;
    Eigen::MatrixXf z = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf lambda = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf b = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf bvt = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf lambda_over_beta = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf lambda_converge = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf pi_converge = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf xy;

    Eigen::MatrixXf h = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf pi = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf rescaled_h = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    float c = mat.rows() * mat.cols();

    for(std::uint32_t i = 0; i < max_iterations; ++i){
        //u and v update
        lambda_over_beta = lambda / beta;
        b = mat - z - lambda_over_beta;

        bvt = b * y.transpose();
        auto qr = bvt.colPivHouseholderQr();
        Eigen::VectorXf diag = qr.matrixR().diagonal();
        std::uint32_t rank_estimate = 0;
        for(; rank_estimate < diag.size(); ++rank_estimate){
            if(fabs(diag(rank_estimate)) < 1e-5){
                break;
            }
        }

        x = qr.householderQ();

        rank_estimate = qr.rank();
        x.conservativeResize(x.rows(), std::min((std::uint32_t)diag.size(), rank_estimate + 1));
        y = x.transpose() * b;

        //z update
        xy = x * y;
        Eigen::MatrixXf pi_over_beta = pi / beta;

        for(std::uint32_t j = 0; j < z.rows(); ++j){
            for(std::uint32_t k = 0; k < z.cols(); ++k){
                //float val = mat(j, k) - xy(j, k) - lambda_over_beta(j, k);
                float val = mat(j, k) - xy(j, k) + h(j, k) - lambda_over_beta(j, k) - pi_over_beta(j, k);
                if(val > 0){
                    z(j, k) = 0.5f * std::max(0.f, val - 1.f / beta);
                }
                else{
                    z(j, k) = 0.5f * std::min(0.f, val + 1.f / beta);
                }
            }
        }

        for(std::uint32_t j = 0; j < h.rows(); ++j){
            for(std::uint32_t k = 0; k < h.cols(); ++k){
                h(j, k) = (1.f - sampled(j, k)) * (z(j, k) + pi_over_beta(j, k));
                rescaled_h(j, k) = h(j, k) * two_q(j, k) - 1.f;
            }
        }
        float rhfrobsq = rescaled_h.norm();
        //std::cout << rhfrobsq << " " << pi_over_beta.norm() << " " << two_q.norm() << " " << h.norm() << std::endl;
        rhfrobsq *= rhfrobsq;
        float gamma = sqrt(c / rhfrobsq);

        for(std::uint32_t j = 0; j < h.rows(); ++j){
            for(std::uint32_t k = 0; k < h.cols(); ++k){
                float inv_twoq = fabs(two_q(j, k)) < 1e-5 ? 0.f : 1.f / two_q(j, k);
                h(j, k) = gamma * (h(j, k) - inv_twoq) + inv_twoq;
            }
        }

        lambda_converge = xy + z - mat;
        pi_converge = z - h;
        float lc = lambda_converge.norm() / mat.norm();
        float pc = pi_converge.norm();
        float lvt = (lambda * y.transpose()).norm();
        float utl = (x.transpose() * lambda).norm();
        //std::cout << i << " " << lc << " " << pc << " " << z.lpNorm<1>() << " " << h.lpNorm<1>() << " " <<
        //    lambda.norm() << " " << pi.norm() << " " << x.cols() << std::endl;
        if(lc < theta && pc < theta){
            //break;
        }

        lambda = lambda + beta * lambda_converge;
        pi = pi + beta * pi_converge;
    }

    for(int i = 0; i < xy.rows(); ++i){
        for(int j = 0; j < xy.cols(); ++j){
            xy(i, j) = std::max(0.f, xy(i, j));
        }
    }

    if(show_rank){
        float val = (float)x.cols() / (float)std::min((int)x.rows(), (int)y.cols());
        float r, g, b;
        std::tie(r, g, b) = floatToRGB(val);
        r /= (float)(xy.cols() / 3);
        g /= (float)(xy.cols() / 3);
        b /= (float)(xy.cols() / 3);

        for(int i = 0; i < xy.rows(); ++i){
            for(int j = 0; j < xy.cols() / 3; ++j){
                xy(i, j * 3) = r;
                xy(i, j * 3 + 1) = g;
                xy(i, j * 3 + 2) = b;
            }
        }
    }

    return std::make_tuple(xy, mat - xy);
}

float gaussian(float x, float mu, float sigma){
    float a = (x - mu) / sigma;
    return std::exp(-0.5f * a * a );
}

Eigen::VectorXf createGaussianKernel(std::uint8_t size){
    Eigen::VectorXf kernel(size);

    float total = 0.f;
    float rad = (float)size / 2.f;
    float sigma = rad / 2.f;

    for(std::uint8_t idx = 0; idx < size; ++idx){
        kernel(idx) = gaussian((float)idx, rad, sigma);
        total += kernel(idx);
    }

    return kernel / total;
}

//performs 2d convolution
void convolve(Eigen::MatrixXf& mat, const Eigen::VectorXf& row_kernel, const Eigen::VectorXf& col_kernel){
    Eigen::MatrixXf temp(mat.rows(), mat.cols());

    for(std::uint32_t row = 0; row < mat.rows(); ++row){
        for(std::uint32_t col = 0; col < mat.cols(); ++col){
            for(std::uint32_t kernel_idx = 0; kernel_idx < row_kernel.size(); ++kernel_idx){
                int uncorrected_idx = kernel_idx - row_kernel.size() / 2 + row;
                int sample_idx = std::min(std::max(uncorrected_idx, (int)mat.rows() - 1), 0);
                temp(row, col) = row_kernel(kernel_idx) * mat(sample_idx, col);
            }
        }
    }

    for(std::uint32_t row = 0; row < mat.rows(); ++row){
        for(std::uint32_t col = 0; col < mat.cols(); ++col){
            for(std::uint32_t kernel_idx = 0; kernel_idx < col_kernel.size(); ++kernel_idx){
                int uncorrected_idx = kernel_idx - col_kernel.size() / 2 + row;
                int sample_idx = std::min(std::max(uncorrected_idx, (int)mat.cols() - 1), 0);
                mat(row, col) = col_kernel(kernel_idx) * temp(row, sample_idx);
            }
        }
    }
}

void reincorporateDenseHighRank(Eigen::MatrixXf& low_rank, const Eigen::MatrixXf& sparse, KDTNode* slice,
    float sparsity_threshold, const Eigen::VectorXf& kernel, float density_threshold, Scene* scene,
    const std::vector<VPL> vpls, float min_dist){

    Eigen::MatrixXf discrete_sparse = Eigen::MatrixXf::Zero(sparse.rows(), sparse.cols());
    for(std::uint32_t row = 0; row < discrete_sparse.rows(); ++row){
        for(std::uint32_t col = 0; col < discrete_sparse.cols(); ++col){
            if(sparse(row, col) > sparsity_threshold){
                discrete_sparse(row, col) = 1;
            }
        }
    }

    convolve(discrete_sparse, kernel, kernel);

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for(std::uint32_t row = 0; row < discrete_sparse.rows(); ++row){
        for(std::uint32_t light = 0; light < discrete_sparse.cols() / 3; ++light){
            bool requires_direct_sample = discrete_sparse(row, light * 3) > density_threshold ||
                discrete_sparse(row, light * 3 + 1) > density_threshold  ||
                discrete_sparse(row, light * 3 + 2) > density_threshold;

            if(requires_direct_sample){
                Spectrum col = sample(scene, sampler, slice->sample(row).its, vpls[light], min_dist, true);

                float r, g, b;
                col.toLinearRGB(r, g, b);

                low_rank(row, light * 3) = r;
                low_rank(row, light * 3 + 1) = g;
                low_rank(row, light * 3 + 2) = b;
            }
        }
    }
}



MatrixSeparationRenderer::MatrixSeparationRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, 
        float min_dist, float sample_percentage, float error_threshold, float reincorporation_density_threshold,
        std::uint32_t slice_size, std::uint32_t max_prediction_iterations, std::uint32_t max_separation_iterations,
        std::uint32_t show_slices, std::uint32_t only_directsamples, bool separate, bool show_error, bool show_sparse,
        std::uint32_t predictor_mask, bool show_rank, bool show_predictors) : 
        clusterer_(std::move(clusterer)), min_dist_(min_dist), sample_percentage_(sample_percentage),
        error_threshold_(error_threshold), reincorporation_density_threshold_(reincorporation_density_threshold),
        slice_size_(slice_size), max_prediction_iterations_(max_prediction_iterations), 
        max_separation_iterations_(max_separation_iterations), show_slices_(show_slices > 0),
        show_only_directsamples_(only_directsamples > 0), cancel_(false), separate_(separate), show_error_(show_error),
        show_sparse_(show_sparse), predictor_mask_(predictor_mask), show_rank_(show_rank), show_predictors_(show_predictors){

}

MatrixSeparationRenderer::MatrixSeparationRenderer(MatrixSeparationRenderer&& other) : clusterer_(std::move(other.clusterer_)),
        min_dist_(other.min_dist_), sample_percentage_(other.sample_percentage_), error_threshold_(other.error_threshold_),
        reincorporation_density_threshold_(other.reincorporation_density_threshold_), slice_size_(other.slice_size_),
        max_prediction_iterations_(other.max_prediction_iterations_),
        max_separation_iterations_(other.max_separation_iterations_), show_slices_(other.show_slices_), 
        show_only_directsamples_(other.show_only_directsamples_), cancel_(false), samples_(std::move(other.samples_)),
        separate_(other.separate_), show_error_(other.show_error_), show_sparse_(other.show_sparse_),
        predictor_mask_(other.predictor_mask_), show_rank_(other.show_rank_), show_predictors_(other.show_predictors_){
    
}

MatrixSeparationRenderer& MatrixSeparationRenderer::operator = (MatrixSeparationRenderer&& other){
    if(this != &other){
        clusterer_ = std::move(other.clusterer_);
        min_dist_ = other.min_dist_;
        sample_percentage_ = other.sample_percentage_;
        error_threshold_ = other.error_threshold_;
        reincorporation_density_threshold_ = other.reincorporation_density_threshold_;
        slice_size_ = other.slice_size_;
        max_prediction_iterations_ = other.max_prediction_iterations_;
        max_separation_iterations_ = other.max_separation_iterations_;
        show_slices_ = other.show_slices_;
        show_only_directsamples_ = other.show_only_directsamples_;
        cancel_ = false;
        samples_ = std::move(other.samples_);
        separate_ = other.separate_;
        show_error_ = other.show_error_;
        show_sparse_ = other.show_sparse_;
        predictor_mask_ = other.predictor_mask_;
        show_rank_ = other.show_rank_;
        show_predictors_ = other.show_predictors_;
    }
    
    return *this;
}

MatrixSeparationRenderer::~MatrixSeparationRenderer(){

}

void calculateNearestVPLNeighbours(const std::vector<VPL>& vpls, std::vector<std::vector<int>>& indices, 
    std::vector<std::vector<float>>& distances, std::uint32_t num_neighbours){
    assert(num_neighbours > 1 && vpls.size() > 0);

    num_neighbours = std::min((std::uint32_t)vpls.size(), num_neighbours + 1);

    flann::Matrix<float> dataset(new float[vpls.size() * 3], vpls.size(), 3);
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        float* curr = (float*)dataset[i];
        curr[0] = vpls[i].its.p.x;
        curr[1] = vpls[i].its.p.y;
        curr[2] = vpls[i].its.p.z;
    }

    flann::Index<flann::L2<float>> index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();
    index.knnSearch(dataset, indices, distances, num_neighbours, flann::SearchParams(128));

    for(std::uint32_t i = 0; i < indices.size(); ++i){
        auto iter = std::find(indices[i].begin(), indices[i].end(), i);
        if(iter != indices[i].end()){
            indices[i].pop_back();
            distances[i].pop_back();
        }
        else{
            std::uint32_t pos = iter - indices[i].begin();
            indices[i].erase(iter);
            distances[i].erase(distances[i].begin() + pos);
        }
    }

    delete dataset.ptr();
}

bool MatrixSeparationRenderer::render(Scene* scene){
    {
        std::lock_guard<std::mutex> lock(cancel_lock_);
        cancel_ = false;
    }

    Intersection its;
    std::vector<VPL> vpls = clusterer_->getClusteringForPoint(its);
    std::vector<std::vector<int>> vpl_nearest_neighbours;
    std::vector<std::vector<float>> vpl_nearest_distances;

    calculateNearestVPLNeighbours(vpls, vpl_nearest_neighbours, vpl_nearest_distances, 10);
    
    std::cout << "constructing kd tree" << std::endl;
    auto root = constructKDTree(scene, vpls, slice_size_, min_dist_, samples_);
    std::vector<KDTNode*> slices;

    getSlices(root.get(), slices);

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        slices[i]->buildNearestNeighbours();
    }

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

    std::cout << "processing slices..." << std::endl;
    if(!show_slices_){
        #pragma omp parallel for
        for(std::uint32_t i = 0; i < slices.size(); ++i){
            std::mt19937 rng(seed * i);

            performInitialVisibilitySamples(slices[i], sample_percentage_, rng, scene, vpls);

            std::set<std::uint32_t> active_cols;
            for(std::uint32_t j = 0; j < vpls.size(); ++j){
                active_cols.insert(j);
            }

            std::uint32_t iter = 0;

            while(active_cols.size() > 0 && iter < max_prediction_iterations_){
                active_cols = sampleAndPredictVisibility(slices[i], sample_percentage_, rng, scene, vpls,  
                    active_cols, error_threshold_, min_dist_, predictor_mask_, vpl_nearest_neighbours);
                iter++;
            }

            Eigen::MatrixXf d, two_q;
            std::tie(d, two_q) = sliceToMatrices(slices[i]);

            Eigen::MatrixXf sampled = Eigen::MatrixXf::Zero(d.rows(), d.cols());

            for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
                for(std::uint32_t k = 0; k < slices[i]->sample(j).visibility.size(); ++k){
                    if(slices[i]->sample(j).visibility[k] == VISIBLE || slices[i]->sample(j).visibility[k] == NOT_VISIBLE){
                        sampled(j, k) = 1.f;
                    }
                }
            }

            if(separate_){
                float beta = 5000.f / d.norm();
                float step = 0.1f;
                float theta = 0.0001f;
                
                Eigen::MatrixXf l, s;
                std::tie(l, s) = separate(d, two_q, max_separation_iterations_, beta, step, theta, sampled, show_rank_);

                auto gaussian_kernel = createGaussianKernel(7);

                //reincorporateDenseHighRank(l, s, slices[i], 0.00001f, gaussian_kernel, reincorporation_density_threshold_, 
                //    scene, vpls, min_dist_);

                matrixToSlice(slices[i], show_sparse_ ? s : l);
            }
            
        }
    }

    std::cout << "generating final image" << std::endl;
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    auto size = film->getSize();
    if(size.x == 0 || size.y == 0){
        return true;
    }

    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, size);
    std::uint8_t* output_image = output_bitmap->getUInt8Data();
    memset(output_image, 0, output_bitmap->getBytesPerPixel() * size.x * size.y);

    std::uint32_t fully_sampled = 0;
    std::uint32_t total = 0;

    //slice iteration is mainly for slice colour
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        int y = rand() % 255;
        int u = rand() % 255;
        int v = ((float)i / slices.size()) * 255.f;

        int sb = 1.164 * (y - 16) + 2.018 * (u - 128);
        int sg = 1.164 * (y - 16) - 0.813 * (v - 128) - 0.391 * (u - 128);
        int sr = 1.164 * (y - 16) + 1.596 * (v - 128);

        total += slices[i]->sample_indices.size() > 0 ? 
                slices[i]->sample_indices.size() * slices[i]->sample(0).col_samples.size() :
                0;

        for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
            std::uint32_t offset = (slices[i]->sample(j).image_x + slices[i]->sample(j).image_y * 
                output_bitmap->getSize().x) * output_bitmap->getBytesPerPixel();

            if(show_slices_){
                output_image[offset] = sr;
                output_image[offset + 1] = sg;
                output_image[offset + 2] = sb;

                continue;
            }

            if(show_predictors_){
                Spectrum r, g, b;
                r.fromLinearRGB(1.f, 0.f, 0.f);
                g.fromLinearRGB(0.f, 1.f, 0.f);
                b.fromLinearRGB(0.f, 0.f, 1.f);

                float rtot = 0.f, gtot = 0.f, btot = 0.f;
                for(std::uint32_t k = 0; k < slices[i]->sample(j).predictors.size(); ++k){
                    if(slices[i]->sample(j).predictors[k] == 0){
                        rtot += 1.f;
                    }
                    else if(slices[i]->sample(j).predictors[k] == 1){
                        gtot += 1.f;
                    }
                    else if(slices[i]->sample(j).predictors[k] == 2){
                        btot += 1.f;
                    }
                }

                rtot /= slices[i]->sample(j).predictors.size();
                gtot /= slices[i]->sample(j).predictors.size();
                btot /= slices[i]->sample(j).predictors.size();

                Spectrum predictor_col = rtot * r + gtot * g + btot * b;

                float pr, pg, pb;
                predictor_col.toSRGB(pr, pg, pb);

                output_image[offset] = pr * 255.f + 0.5f;
                output_image[offset + 1] = pg * 255.f + 0.5f;
                output_image[offset + 2] = pb * 255.f + 0.5f;

                continue;
            }

            Spectrum s(0.f);
            
            if(slices[i]->sample(j).its.isEmitter()){
                s = slices[i]->sample(j).emitter_color;
            }
            else{
                for(std::uint32_t k = 0; k < slices[i]->sample(j).col_samples.size(); ++k){
                    if(slices[i]->sample(j).visibility[k] == NOT_VISIBLE || slices[i]->sample(j).visibility[k] == VISIBLE){
                        fully_sampled++;
                    }

                    if(!separate_){
                        if(slices[i]->sample(j).visibility[k] == NOT_VISIBLE ||
                            (slices[i]->sample(j).visibility[k] == P_NOT_VISIBLE && !show_only_directsamples_)){
                            continue;
                        }
                    }

                    s += slices[i]->sample(j).col_samples[k];
                }
            }

            
            float r, g, b;
            s.toSRGB(r, g, b);

            output_image[offset] = std::min(1.f, r) * 255 + 0.5f;
            output_image[offset + 1] = std::min(1.f, g) * 255 + 0.5f;
            output_image[offset + 2] = std::min(1.f, b) * 255 + 0.5f;
        }
    }

    if(show_error_){
        float error = calculateError(scene, vpls, min_dist_, output_image);
        std::cout << "Error: " << error << std::endl;
    }

    film->setBitmap(output_bitmap);

    if(!show_slices_){
        std::cout << fully_sampled << " " << total << std::endl;
        std::cout << "Fully sampled: " << (float)fully_sampled / total << std::endl;
    }
    
    return !cancel_;
}

MTS_NAMESPACE_END