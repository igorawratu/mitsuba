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

#include <mitsuba/core/plugin.h>

#include <eigen3/Eigen/Dense>

#include "definitions.h"
#include "nanoflann.hpp"

MTS_NAMESPACE_BEGIN

enum VISIBILITY{VISIBLE, NOT_VISIBLE, P_VISIBLE, P_NOT_VISIBLE, UNKNOWN};

struct RowSample{
    RowSample(const Point3f& p, const Normal& n, std::uint32_t x, std::uint32_t y, std::uint32_t cols, bool in_scene, Intersection intersection) : 
        position(p), normal(n), image_x(x), image_y(y), col_samples(cols, Spectrum(0.f)), visibility(cols, UNKNOWN), intersected_scene(in_scene), its(intersection){
    }

    Point3f position;
    Normal normal;
    std::uint32_t image_x, image_y;
    std::vector<Spectrum> col_samples;
    std::vector<VISIBILITY> visibility;
    bool intersected_scene;
    Intersection its;
};

struct KDTNode{
    KDTNode() : samples(), left(nullptr), right(nullptr){
    }

    std::vector<RowSample> samples;
    std::unique_ptr<KDTNode> left;
    std::unique_ptr<KDTNode> right;

    void Split(float norm_scale, std::uint32_t size_threshold){
        float maxf = std::numeric_limits<float>::max();
        float minf = std::numeric_limits<float>::min();
        Vector3f min_pos(maxf, maxf, maxf), max_pos(minf, minf, minf);
        Vector3f min_normal(maxf, maxf, maxf), max_normal(minf, minf, minf);

        for(size_t i = 0; i < samples.size(); ++i){
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            min_pos.y = std::min(samples[i].position.y, min_pos.y);
            min_pos.z = std::min(samples[i].position.z, min_pos.z);
            max_pos.x = std::max(samples[i].position.x, max_pos.x);
            max_pos.y = std::max(samples[i].position.y, max_pos.y);
            max_pos.z = std::max(samples[i].position.z, max_pos.z);

            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            min_normal.y = std::min(samples[i].normal.y, min_normal.y);
            min_normal.z = std::min(samples[i].normal.z, min_normal.z);
            max_normal.x = std::max(samples[i].normal.x, max_normal.x);
            max_normal.y = std::max(samples[i].normal.y, max_normal.y);
            max_normal.z = std::max(samples[i].normal.z, max_normal.z);
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

        std::array<std::function<bool(const RowSample& lhs, const RowSample& rhs)>, 6> sorters;
        sorters[0] =    [](const RowSample& lhs, const RowSample& rhs){
                            return lhs.position.x < rhs.position.x;
                        };
        sorters[1] =    [](const RowSample& lhs, const RowSample& rhs){
                            return lhs.position.y < rhs.position.y;
                        };
        sorters[2] =    [](const RowSample& lhs, const RowSample& rhs){
                            return lhs.position.z < rhs.position.z;
                        };
        sorters[3] =    [](const RowSample& lhs, const RowSample& rhs){
                            return lhs.normal.x < rhs.normal.x;
                        };
        sorters[4] =    [](const RowSample& lhs, const RowSample& rhs){
                            return lhs.normal.y < rhs.normal.y;
                        };
        sorters[5] =    [](const RowSample& lhs, const RowSample& rhs){
                            return lhs.normal.z < rhs.normal.z;
                        };

        std::array<std::function<bool(float, const RowSample&)>, 6> searchers;
        searchers[0] =  [](float value, const RowSample& entry){
                            return value < entry.position.x;
                        };
        searchers[1] =  [](float value, const RowSample& entry){
                            return value < entry.position.y;
                        };
        searchers[2] =  [](float value, const RowSample& entry){
                            return value < entry.position.z;
                        };
        searchers[3] =  [](float value, const RowSample& entry){
                            return value < entry.normal.x;
                        };
        searchers[4] =  [](float value, const RowSample& entry){
                            return value < entry.normal.y;
                        };
        searchers[5] =  [](float value, const RowSample& entry){
                            return value < entry.normal.z;
                        };

        std::sort(ranges.begin(), ranges.end(), 
            [](const std::pair<std::uint8_t, float>& lhs, const std::pair<std::uint8_t, float>& rhs){
                return lhs.second > rhs.second;
            });

        std::sort(samples.begin(), samples.end(), sorters[ranges[0].first]);
        
        auto split_iter = std::upper_bound(samples.begin(), samples.end(), midpoints[ranges[0].first], 
            searchers[ranges[0].first]);

        std::uint32_t split_index = split_iter - samples.begin();
        std::uint32_t smallest_slice_size = std::min(split_index, std::uint32_t(samples.size() - split_index));

        //for degenerate cases where one child is too small
        split_index = smallest_slice_size < (size_threshold / 10) ? samples.size() / 2 : split_index;

        left = std::unique_ptr<KDTNode>(new KDTNode());
        right = std::unique_ptr<KDTNode>(new KDTNode());

        left->samples.insert(left->samples.end(), samples.begin(), samples.begin() + split_index);
        right->samples.insert(right->samples.end(), samples.begin() + split_index, samples.end());

        samples.clear();
    }
};

void splitKDTree(KDTNode* node, std::uint32_t size_threshold, float min_dist){
    if(node == nullptr || node->samples.size() < size_threshold){
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

Spectrum sample(Scene* scene, Sampler* sampler, const Intersection& its, const VPL& vpl, const Point& position, 
    const Normal& normal, float min_dist, bool check_occlusion){

    //only dealing with emitter and surface VPLs curently.
    if (vpl.type != EPointEmitterVPL && vpl.type != ESurfaceVPL){
        return Spectrum(0.f);
    }

    if(check_occlusion){
        Ray shadow_ray(position, normalize(vpl.its.p - position), 0.f);

        Float t;
        ConstShapePtr shape;
        Normal norm;
        Point2 uv;

        if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
            if(abs((position - vpl.its.p).length() - t) > 0.0001f ){
                return Spectrum(0.f);
            }
        }
    }

    BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
        bsdf_sample_record.wi = its.toLocal(normalize(vpl.its.p - its.p));
        bsdf_sample_record.wo = its.toLocal(normal);

    Spectrum albedo = its.getBSDF()->eval(bsdf_sample_record);

    float d = std::max((its.p - vpl.its.p).length(), min_dist);
    float attenuation = 1.f / (d * d);

    float n_dot_ldir = std::max(0.f, dot(normalize(normal), normalize(vpl.its.p - its.p)));
    float ln_dot_ldir = std::max(0.f, dot(normalize(vpl.its.shFrame.n), normalize(its.p - vpl.its.p)));

    return (vpl.P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
}

std::unique_ptr<KDTNode> constructKDTree(Scene* scene, const std::vector<VPL>& vpls, 
    std::uint32_t size_threshold, float min_dist){

    auto kdt_root = std::unique_ptr<KDTNode>(new KDTNode());

    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for (std::int32_t y = 0; y < film->getSize().y; ++y) {
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

            kdt_root->samples.emplace_back(its.p, its.geoFrame.n, x, y, vpls.size(), intersected, its);

            auto& curr_sample = kdt_root->samples.back();

            if(its.isEmitter()){
                std::fill(curr_sample.col_samples.begin(), 
                    curr_sample.col_samples.end(), its.Le(-ray.d));
                continue;
            }

            for (std::uint32_t i = 0; i < vpls.size(); ++i) {
                curr_sample.col_samples[i] = sample(scene, sampler, its, vpls[i], its.p, 
                    its.geoFrame.n, min_dist, false);
            }
        }
    }

    splitKDTree(kdt_root.get(), size_threshold, min_dist);

    return kdt_root;
}

struct FLANNPoint
{
    Point position;
    std::uint32_t idx;

    inline float get_pt(const size_t dim) const
	{
		if (dim == 0) return position.x;
		else if (dim == 1) return position.y;
		else return position.z;
	}
};

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, FLANNPointCloud<FLANNPoint>>, 
    FLANNPointCloud<FLANNPoint>, 3> KDTree;

std::vector<VISIBILITY> knnPredictor(KDTNode* slice, std::uint32_t neighbours, std::uint32_t col, std::mt19937& rng){
    assert(slice != nullptr && neighbours > 0 && slice->samples.size() > 0
        && col < slice->samples[0].visibility.size());

    FLANNPointCloud<FLANNPoint> point_set;
    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        if(slice->samples[i].visibility[col] == VISIBLE || slice->samples[i].visibility[col] == NOT_VISIBLE){
            FLANNPoint p;
            p.position = slice->samples[i].position;
            p.idx = i;

            point_set.pts.push_back(p);
        }
    }

    KDTree kdt(3, point_set, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    kdt.buildIndex();

    std::uint32_t num_neighbours = std::min(neighbours, (std::uint32_t)point_set.pts.size());
    std::vector<VISIBILITY> output(slice->samples.size());

    if(num_neighbours == 0){
        for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
            output[i] = P_VISIBLE;
        }

        return output;
    }

    std::uniform_real_distribution<float> gen(0.f, 1.f);
    nanoflann::KNNResultSet<float> result_set(num_neighbours);
    std::unique_ptr<size_t[]> indices(new size_t[num_neighbours]);
    std::unique_ptr<float[]> distances(new float[num_neighbours]);

    result_set.init(indices.get(), distances.get());

    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        if(slice->samples[i].visibility[col] == VISIBLE || slice->samples[i].visibility[col] == NOT_VISIBLE){
            output[i] = slice->samples[i].visibility[col];
        }
        else{
            float query_pt[3] = 
                {slice->samples[i].position.x, slice->samples[i].position.y, slice->samples[i].position.z};
            kdt.findNeighbors(result_set, query_pt, nanoflann::SearchParams(10));

            std::uint32_t num_visible = 0;
            for(std::uint32_t j = 0; j < num_neighbours; ++j){
                std::uint32_t sample_index = point_set.pts[indices[j]].idx;
                if(slice->samples[sample_index].visibility[col] == VISIBLE){
                    num_visible++;
                }
            }

            float p_v = (float)num_visible / num_neighbours;
            output[i] = gen(rng) < p_v ? P_VISIBLE : P_NOT_VISIBLE;
        }
    }

    return output;
}
std::vector<VISIBILITY> linearPredictor(KDTNode* slice, std::uint32_t col, float min_dist, std::mt19937& rng){
    assert(slice != nullptr && slice->samples.size() > 0);

    Eigen::MatrixXf q(slice->samples.size(), 4);

    std::vector<std::uint32_t> sampled;
    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        if(slice->samples[i].visibility[col] == VISIBLE || slice->samples[i].visibility[col] == NOT_VISIBLE){
            sampled.push_back(i);
        }
        q(i, 0) = slice->samples[i].position.x / min_dist;
        q(i, 1) = slice->samples[i].position.y / min_dist;
        q(i, 2) = slice->samples[i].position.z / min_dist;
        q(i, 3) = 1;
    }

    Eigen::MatrixXf x(sampled.size(), 4);
    Eigen::MatrixXf y(sampled.size(), 1);

    for(std::uint32_t i = 0; i < sampled.size(); ++i){
        x(i, 0) = slice->samples[sampled[i]].position.x / min_dist;
        x(i, 1) = slice->samples[sampled[i]].position.y / min_dist;
        x(i, 2) = slice->samples[sampled[i]].position.z / min_dist;
        x(i, 3) = 1;

        y(i, 0) = slice->samples[sampled[i]].visibility[col] == VISIBLE ? 1 : -1;
    }

    Eigen::MatrixXf xxt = x.transpose() * x;
    Eigen::MatrixXf w = (x.transpose() * x).inverse() * x.transpose() * y;

    Eigen::MatrixXf fq = q * w;

    std::vector<VISIBILITY> output(slice->samples.size());

    std::uniform_real_distribution<float> gen(0.f, 1.f);

    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        if(slice->samples[i].visibility[col] == VISIBLE || slice->samples[i].visibility[col] == NOT_VISIBLE){
            output[i] = slice->samples[i].visibility[col];
        }
        else{
            float pv0 = std::min(1.f, std::max(0.f, 0.5f - fq(i, 0)));
            output[i] = gen(rng) < pv0 ? P_NOT_VISIBLE : P_VISIBLE;
            output[i] = fq(i, 0) > 0 ? P_VISIBLE : P_NOT_VISIBLE;
        }
    }

    return output;
}

std::vector<VISIBILITY> naiveBayes(KDTNode* slice, std::uint32_t col, std::uint32_t neighbours, 
    const std::vector<VPL>& vpls, std::mt19937& rng){

    assert(slice != nullptr && slice->samples.size() > 0 && slice->samples[0].visibility.size() > 0 && 
        vpls.size() > col);
    
    std::uint32_t num_neighbours = std::min(neighbours, (std::uint32_t)vpls.size() - 1);
    std::vector<VISIBILITY> output(slice->samples.size());
    if(num_neighbours == 0){
        for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
            output[i] = P_VISIBLE;
        }

        return output;
    }

    FLANNPointCloud<FLANNPoint> point_set;
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        if(i == col){
            continue;
        }

        FLANNPoint p;
        p.position = vpls[i].its.p;
        p.idx = i;

        point_set.pts.push_back(p);
    }

    KDTree kdt(3, point_set, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    kdt.buildIndex();

    std::uniform_real_distribution<float> gen(0.f, 1.f);
    nanoflann::KNNResultSet<float> result_set(num_neighbours);
    std::unique_ptr<size_t[]> indices(new size_t[num_neighbours]);
    std::unique_ptr<float[]> distances(new float[num_neighbours]);

    result_set.init(indices.get(), distances.get());

    float query_pt[3] = {vpls[col].its.p.x, vpls[col].its.p.y, vpls[col].its.p.z};
    kdt.findNeighbors(result_set, query_pt, nanoflann::SearchParams(10));

    float tot_vis = 0.f;
    float n_vis = 0.f;
    
    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        for(std::uint32_t j = 0; j < num_neighbours; ++j){
            VISIBILITY v = slice->samples[i].visibility[point_set.pts[indices[j]].idx];
            if(v == P_VISIBLE || v == VISIBLE){
                tot_vis += 1.f;
                n_vis += 1.f;
            }
        }

        VISIBILITY v = slice->samples[i].visibility[col];
        if(v == P_VISIBLE || v == VISIBLE){
            tot_vis += 1.f;
        }
    }

    float p_v = tot_vis / (float)(slice->samples.size() * num_neighbours + 1);
    float pnj_v = n_vis / tot_vis;
    float p_i = 1.f / slice->samples.size();
    float p_nj = (float)num_neighbours / (num_neighbours + 1);

    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        float tot_iv = 0.f;
        for(std::uint32_t j = 0; j < num_neighbours; ++j){
            VISIBILITY v = slice->samples[i].visibility[point_set.pts[indices[j]].idx];
            if(v == P_VISIBLE || v == VISIBLE){
                tot_iv += 1.f;
            }
        }

        VISIBILITY v = slice->samples[i].visibility[col];
        if(v == P_VISIBLE || v == VISIBLE){
            tot_iv += 1.f;
        }

        float pi_v = tot_iv / tot_vis;
        float pv_ij = (pi_v * pnj_v * p_v) / (p_i * p_nj);

        output[i] = gen(rng) < pv_ij ? P_VISIBLE : P_NOT_VISIBLE;
    }

    return output;
}

void performInitialVisibilitySamples(KDTNode* slice, float sample_percentage, std::mt19937& rng, 
    Scene* scene, const std::vector<VPL>& vpls){

    assert(slice != nullptr);

    std::uint32_t num_samples = sample_percentage * slice->samples.size();
    
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
        std::vector<std::uint32_t> sample_set(slice->samples.size());
        for(std::uint32_t j = 0; j < slice->samples.size(); ++j){
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
            Point ray_origin = slice->samples[to_sample[j]].position;
            Ray shadow_ray(ray_origin, normalize(vpls[i].its.p - ray_origin), 0.f);

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            slice->samples[to_sample[j]].visibility[i] = VISIBLE;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if(abs((ray_origin - vpls[i].its.p).length() - t) > 0.0001f ){
                    slice->samples[to_sample[j]].visibility[i] = NOT_VISIBLE;
                }
            }
        }
    }
}

std::set<std::uint32_t> sampleAndPredictVisibility(KDTNode* slice, float sample_percentage, std::mt19937& rng, Scene* scene,
    const std::vector<VPL>& vpls, const std::set<std::uint32_t>& active_columns, float error_threshold, float min_dist){

    std::set<std::uint32_t> new_active_cols;

    std::uint32_t num_samples = sample_percentage * slice->samples.size();
    std::vector<std::vector<VISIBILITY>> new_predictions(active_columns.size());

    for(auto iter = active_columns.begin(); iter != active_columns.end(); iter++){
        std::vector<std::uint32_t> unsampled;
        for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
            auto& curr_samples = slice->samples[i];
            if(curr_samples.visibility[*iter] != VISIBLE && curr_samples.visibility[*iter] != NOT_VISIBLE){
                unsampled.push_back(i);
            }
        }

        std::vector<std::vector<VISIBILITY>> predictions;
        
        predictions.push_back(linearPredictor(slice, *iter, min_dist, rng));
        predictions.push_back(naiveBayes(slice, *iter, 3, vpls, rng));
        predictions.push_back(knnPredictor(slice, 3, *iter, rng));

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
            Point ray_origin = slice->samples[to_sample[i]].position;
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

        for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
            slice->samples[i].visibility[*iter] = predictions[smallest_idx][i];
        }

        float normalized_error = (float)smallest_error / to_sample.size();
        if(normalized_error > error_threshold){
            new_active_cols.insert(*iter);
        }
    }

    return new_active_cols;
}

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> sliceToMatrices(KDTNode* node){
    assert(node->samples.size() > 0 && node->samples[0].col_samples.size() > 0);

    Eigen::MatrixXf slice_mat(node->samples.size(), node->samples[0].col_samples.size() * 3);
    Eigen::MatrixXf error_mat(node->samples.size(), node->samples[0].col_samples.size() * 3);

    for(std::uint32_t row = 0; row < node->samples.size(); ++row){
        for(std::uint32_t col = 0; col < node->samples[row].col_samples.size(); ++col){
            float r, g , b;
            node->samples[row].col_samples[col].toLinearRGB(r, g, b);
            if(node->samples[row].visibility[col] == VISIBLE || node->samples[row].visibility[col] == P_VISIBLE){
                slice_mat(row, col * 3) = r;
                slice_mat(row, col * 3 + 1) = g;
                slice_mat(row, col * 3 + 2) = b;
            }
            else{
                slice_mat(row, col * 3) = slice_mat(row, col * 3 + 1) = slice_mat(row, col * 3 + 2) = 0.f;
            }
            
            error_mat(row, col * 3) = 2.f / r;
            error_mat(row, col * 3 + 1) = 2.f / g;
            error_mat(row, col * 3 + 2) = 2.f / b;
        }
    }

    return std::make_tuple(slice_mat, error_mat);
}

void matrixToSlice(KDTNode* node, const Eigen::MatrixXf& mat){
    assert(node->samples.size() == (std::uint32_t)mat.rows() && 
        (node->samples[0].col_samples.size() * 3) == (std::uint32_t)mat.cols());

    for(std::uint32_t row = 0; row < node->samples.size(); ++row){
        for(std::uint32_t col = 0; col < node->samples[row].col_samples.size(); ++col){
            float r = mat(row, col * 3);
            float g = mat(row, col * 3 + 1);
            float b = mat(row, col * 3 + 2);

            node->samples[row].col_samples[col].fromLinearRGB(r, g, b);
        }
    }
}

void shrink(Eigen::MatrixXf& mat, float amount){
    for(std::uint32_t i = 0; i < mat.rows(); ++i){
        for(std::uint32_t j = 0; j < mat.cols(); ++j){
            int sign = mat(i, j) < 0 ? -1 : 1;
            mat(i, j) = sign * std::max(0.f, abs(mat(i, j)) - amount);
        }
    }
}

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> separate(const Eigen::MatrixXf& mat, const Eigen::MatrixXf& two_q,
    const std::uint32_t max_iterations, float beta, const std::set<std::pair<std::uint32_t, std::uint32_t>>& sampled){
    float c = two_q.rows() * two_q.cols();

    std::uint32_t rank_estimate = mat.bdcSvd().nonzeroSingularValues();;
    Eigen::MatrixXf identity = Eigen::MatrixXf::Identity(two_q.rows(), two_q.cols());
    Eigen::MatrixXf inv_two_q = two_q.cwiseInverse();
    Eigen::MatrixXf x;
    Eigen::MatrixXf y = Eigen::MatrixXf::Zero(rank_estimate, mat.cols());
    Eigen::MatrixXf z = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf h = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf lambda = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf pi = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());

    Eigen::MatrixXf b;

    for(std::uint32_t iteration = 0; iteration < max_iterations; ++iteration){
        b = z - lambda/beta;
        Eigen::MatrixXf ru = b * y.transpose();
        
        x = ru.householderQr().householderQ();
        x = x.block(0, 0, x.rows(), ru.cols());
        y = x.transpose() * b;

        z = 0.5 * (mat + h - x * y - (lambda + pi) / beta);
        shrink(z, 1.f / beta);

        h = z + pi / beta;
        
        for(auto iter = sampled.begin(); iter != sampled.end(); ++iter){
            h(iter->first, iter->second) = 0.f;
        }

        float gamma = sqrt(c / (two_q.cwiseProduct(h) - identity).squaredNorm());

        h = gamma * (h - inv_two_q) + inv_two_q;

        lambda = lambda + beta * (z + x * y - mat);
        pi = pi + beta * (z - h);

        //check for convergence???
    }

    return std::make_tuple(x * y, z);
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
                Spectrum col = sample(scene, sampler, slice->samples[row].its, vpls[light], 
                slice->samples[row].position, slice->samples[row].normal, min_dist, true);

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
        std::uint32_t show_slices, std::uint32_t only_directsamples) : 
        clusterer_(std::move(clusterer)), min_dist_(min_dist), sample_percentage_(sample_percentage),
        error_threshold_(error_threshold), reincorporation_density_threshold_(reincorporation_density_threshold),
        slice_size_(slice_size), max_prediction_iterations_(max_prediction_iterations), 
        max_separation_iterations_(max_separation_iterations), show_slices_(show_slices > 0),
        show_only_directsamples_(only_directsamples > 0), cancel_(false){

}

MatrixSeparationRenderer::MatrixSeparationRenderer(MatrixSeparationRenderer&& other) : clusterer_(std::move(other.clusterer_)),
        min_dist_(other.min_dist_), sample_percentage_(other.sample_percentage_), error_threshold_(other.error_threshold_),
        reincorporation_density_threshold_(other.reincorporation_density_threshold_), slice_size_(other.slice_size_),
        max_prediction_iterations_(other.max_prediction_iterations_),
        max_separation_iterations_(other.max_separation_iterations_), show_slices_(other.show_slices_), 
        show_only_directsamples_(other.show_only_directsamples_), cancel_(false){
    
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
    }
    
    return *this;
}

MatrixSeparationRenderer::~MatrixSeparationRenderer(){

}

bool MatrixSeparationRenderer::render(Scene* scene){
    Intersection its;
    std::vector<VPL> vpls = clusterer_->getClusteringForPoint(its);
    
    std::cout << "constructing kd tree" << std::endl;
    auto root = constructKDTree(scene, vpls, slice_size_, min_dist_);
    std::vector<KDTNode*> slices;

    std::cout << "getting slices" << std::endl;
    getSlices(root.get(), slices);

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    for(std::uint32_t i = 0; i < slices.size() && !show_slices_; ++i){
        std::cout << "slice " << i << " of " << slices.size() << std::endl;

        std::cout << "initial visibility samples" << std::endl;
        performInitialVisibilitySamples(slices[i], sample_percentage_, rng, scene, vpls);

        std::set<std::uint32_t> active_cols;
        for(std::uint32_t j = 0; j < vpls.size(); ++j){
            active_cols.insert(j);
        }

        std::uint32_t iter = 0;

        std::cout << "predicting visibility samples" << std::endl;
        while(active_cols.size() > 0 && iter < max_prediction_iterations_){
            active_cols = sampleAndPredictVisibility(slices[i], sample_percentage_, rng, scene, vpls,  
                active_cols, error_threshold_, min_dist_);
            iter++;
        }

        std::set<std::pair<std::uint32_t, std::uint32_t>> sampled;
        for(std::uint32_t j = 0; j < slices[i]->samples.size(); ++j){
            for(std::uint32_t k = 0; k < slices[i]->samples[j].visibility.size(); ++k){
                if(slices[i]->samples[j].visibility[k] == VISIBLE || slices[i]->samples[j].visibility[k] == NOT_VISIBLE){
                    sampled.insert(std::make_pair(j, k));
                }
            }
        }


        /*std::cout << "separating" << std::endl;
        Eigen::MatrixXf d, two_q;
        std::tie(d, two_q) = sliceToMatrices(slices[i]);
        
        float beta = 500.f / d.norm();
        
        Eigen::MatrixXf l, s;
        std::tie(l, s) = separate(d, two_q, max_separation_iterations_, beta, sampled);

        std::cout << "reincorporating high rank info" << std::endl;
        auto gaussian_kernel = createGaussianKernel(7);

        reincorporateDenseHighRank(l, s, slices[i], 0.00001f, gaussian_kernel, reincorporation_density_threshold_, 
            scene, vpls, min_dist_);

        matrixToSlice(slices[i], l);*/
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

    std::uint32_t fully_sampled = 0;
    std::uint32_t total = 0;

    for(std::uint32_t i = 0; i < slices.size(); ++i){
        int y = rand() % 255;
        int u = rand() % 255;
        int v = ((float)i / slices.size()) * 255.f;

        int sb = 1.164 * (y - 16.f) + 2.018 * (u - 128.f);
        int sg = 1.164 * (y - 16) - 0.813 * (v - 128) - 0.391 * (u - 128);
        int sr = 1.164 * (y - 16) + 1.596 * (v - 128);

        for(std::uint32_t j = 0; j < slices[i]->samples.size(); ++j){
            std::uint32_t offset = (slices[i]->samples[j].image_x + slices[i]->samples[j].image_y * 
                output_bitmap->getSize().x) * output_bitmap->getBytesPerPixel();

            if(show_slices_){
                output_image[offset] = sr;
                output_image[offset + 1] = sg;
                output_image[offset + 2] = sb;

                continue;
            }

            Spectrum s(0.f);

            for(std::uint32_t k = 0; k < slices[i]->samples[j].col_samples.size(); ++k){
                total++;
                if(slices[i]->samples[j].visibility[k] == NOT_VISIBLE || slices[i]->samples[j].visibility[k] == VISIBLE){
                    fully_sampled++;
                }

                if(slices[i]->samples[j].visibility[k] == NOT_VISIBLE ||
                    (slices[i]->samples[j].visibility[k] == P_NOT_VISIBLE && !show_only_directsamples_)){
                    continue;
                }

                s += slices[i]->samples[j].col_samples[k];

            }

            float r, g, b;
            s.toSRGB(r, g, b);

            output_image[offset] = std::min(1.f, r) * 255 + 0.5f;
            output_image[offset + 1] = std::min(1.f, g) * 255 + 0.5f;
            output_image[offset + 2] = std::min(1.f, b) * 255 + 0.5f;
        }
    }

    film->setBitmap(output_bitmap);

    std::cout << "Fully sampled: " << (float)fully_sampled / total << std::endl;

    return false;
}

MTS_NAMESPACE_END