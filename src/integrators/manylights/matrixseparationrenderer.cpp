#include "matrixseparationrenderer.h"

#include <vector>
#include <limits>
#include <utility>
#include <memory>
#include <random>
#include <tuple>
#include <cmath>

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

    void Split(float norm_scale){
        float maxf = std::numeric_limits<float>::max();
        float minf = std::numeric_limits<float>::min();
        Vector3f min_pos(maxf, maxf, maxf), max_pos(minf, minf, minf);
        Vector3f min_normal(maxf, maxf, maxf), max_normal(minf, minf, minf);

        for(size_t i = 0; i < samples.size(); ++i){
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            max_pos.x = std::min(samples[i].position.x, max_pos.x);
            max_pos.x = std::min(samples[i].position.x, max_pos.x);
            max_pos.x = std::min(samples[i].position.x, max_pos.x);

            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            max_normal.x = std::min(samples[i].normal.x, max_normal.x);
            max_normal.x = std::min(samples[i].normal.x, max_normal.x);
            max_normal.x = std::min(samples[i].normal.x, max_normal.x);
        }

        std::array<std::pair<std::uint8_t, float>, 6> ranges;
        ranges[0] = std::make_pair(0, max_pos.x - min_pos.x);
        ranges[1] = std::make_pair(1, max_pos.y - min_pos.y);
        ranges[2] = std::make_pair(2, max_pos.z - min_pos.z);
        ranges[3] = std::make_pair(3, (max_normal.x - min_normal.x) * norm_scale);
        ranges[4] = std::make_pair(4, (max_normal.x - min_normal.x) * norm_scale);
        ranges[5] = std::make_pair(5, (max_normal.x - min_normal.x) * norm_scale);

        std::sort(ranges.begin(), ranges.end(), 
            [](const std::pair<std::uint8_t, float>& lhs, const std::pair<std::uint8_t, float>& rhs){
                return lhs.second > rhs.second;
            });

        switch(ranges[0].first){
            case 0:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.x > rhs.position.x;
                });
                break;
            case 1:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.y > rhs.position.y;
                });
                break;
            case 2:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.z > rhs.position.z;
                });
                break;
            case 3:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.normal.x > rhs.normal.x;
                });
                break;
            case 4:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.normal.y > rhs.normal.y;
                });
                break;
            case 5:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.normal.z > rhs.normal.z;
                });
                break;
            default:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.x > rhs.position.x;
                });
                break;
        }

        left = std::unique_ptr<KDTNode>(new KDTNode());
        right = std::unique_ptr<KDTNode>(new KDTNode());

        std::copy(samples.begin(), samples.begin() + samples.size() / 2, left->samples.begin());
        std::copy(samples.begin() + samples.size() / 2, samples.end(), right->samples.begin());
    }
};

void splitKDTree(KDTNode* node, std::uint32_t size_threshold, float min_dist){
    if(node == nullptr || node->samples.size() < size_threshold){
        return;
    }

    if(node->left == nullptr && node->right == nullptr){
        node->Split(min_dist);
    }
    
    splitKDTree(node->left.get(), size_threshold, min_dist);
    splitKDTree(node->right.get(), size_threshold, min_dist);
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
            kdt_root->samples.emplace_back(its.p, its.geoFrame.n, x, y, vpls.size(), intersected, its);

            auto& curr_sample = kdt_root->samples.back();

            if(!intersected){
                continue;
            }

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
std::vector<VISIBILITY> linearPredictor(KDTNode* slice, std::uint32_t col, std::mt19937& rng){
    assert(slice != nullptr && slice->samples.size() > 0);

    Eigen::MatrixXf q(slice->samples.size(), 4);

    std::uint32_t num_sampled_values = 0;
    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        if(slice->samples[i].visibility[col] == VISIBLE || slice->samples[i].visibility[col] == NOT_VISIBLE){
            num_sampled_values++;
        }
        q(i, 0) = slice->samples[i].position.x;
        q(i, 1) = slice->samples[i].position.y;
        q(i, 2) = slice->samples[i].position.z;
        q(i, 3) = 1;
    }

    Eigen::MatrixXf x(num_sampled_values, 4);
    Eigen::MatrixXf y(num_sampled_values, 1);

    for(std::uint32_t i = 0; i < num_sampled_values; ++i){
        x(i, 0) = slice->samples[i].position.x;
        x(i, 1) = slice->samples[i].position.y;
        x(i, 2) = slice->samples[i].position.z;
        x(i, 3) = 1;

        y(i, 0) = slice->samples[i].visibility[col] == VISIBLE ? 1 : -1;
    }

    Eigen::MatrixXf w = (x.transpose() * x).inverse() * x * y;
    Eigen::MatrixXf fq = q * w;

    std::vector<VISIBILITY> output(slice->samples[0].visibility.size());

    std::uniform_real_distribution<float> gen(0.f, 1.f);

    for(std::uint32_t i = 0; i < slice->samples.size(); ++i){
        if(slice->samples[i].visibility[col] == VISIBLE || slice->samples[i].visibility[col] == NOT_VISIBLE){
            output[i] = slice->samples[i].visibility[col];
        }
        else{
            float pv0 = std::min(1.f, std::max(0.f, 0.5f - fq(i, 0)));
            output[i] = gen(rng) < pv0 ? NOT_VISIBLE : VISIBLE;
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
            else if(v == UNKNOWN){
                tot_vis += 0.5f;
                n_vis += 0.5f;
            }
        }

        VISIBILITY v = slice->samples[i].visibility[col];
        if(v == P_VISIBLE || v == VISIBLE){
            tot_vis += 1.f;
        }
        else if(v == UNKNOWN){
            tot_vis += 0.5f;
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
            else if(v == UNKNOWN){
                tot_iv += 0.5f;
            }
        }

        VISIBILITY v = slice->samples[i].visibility[col];
        if(v == P_VISIBLE || v == VISIBLE){
            tot_iv += 1.f;
        }
        else if(v == UNKNOWN){
            tot_iv += 0.5f;
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
    const std::vector<VPL>& vpls, const std::set<std::uint32_t>& active_columns, float error_threshold){

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
        
        predictions.push_back(knnPredictor(slice, 3, *iter, rng));
        predictions.push_back(linearPredictor(slice, *iter, rng));
        predictions.push_back(naiveBayes(slice, *iter, 3, vpls, rng));

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

        float normalized_error = (float)smallest_error / slice->samples.size();
        if(normalized_error > error_threshold){
            new_active_cols.insert(*iter);
        }
    }

    return new_active_cols;
}

Eigen::MatrixXf sliceToMatrix(KDTNode* node){
    assert(node->samples.size() > 0 && node->samples[0].col_samples.size() > 0);

    Eigen::MatrixXf slice_mat(node->samples.size(), node->samples[0].col_samples.size() * 3);

    for(std::uint32_t row = 0; row < node->samples.size(); ++row){
        for(std::uint32_t col = 0; col < node->samples[row].col_samples.size(); ++col){
            float r = 0.f, g = 0.f, b = 0.f;

            if(node->samples[row].visibility[col] == VISIBLE || node->samples[row].visibility[col] == P_VISIBLE){
                node->samples[row].col_samples[col].toLinearRGB(r, g, b);
            }
            
            slice_mat(row, col * 3) = r;
            slice_mat(row, col * 3 + 1) = g;
            slice_mat(row, col * 3 + 2) = b;
        }
    }

    return slice_mat;
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
    Eigen::MatrixXf y(rank_estimate, mat.cols());
    Eigen::MatrixXf z(mat.rows(), mat.cols());
    Eigen::MatrixXf h(mat.rows(), mat.cols());
    Eigen::MatrixXf lambda(mat.rows(), mat.cols());
    Eigen::MatrixXf pi(mat.rows(), mat.cols());

    Eigen::MatrixXf b;

    for(std::uint32_t iteration = 0; iteration < max_iterations; ++iteration){
        b = z - lambda/beta;
        Eigen::MatrixXf ru = b * y.transpose();
        
        x = ru.householderQr().householderQ();
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

bool MatrixSeparationRenderer::Render(Scene* scene){
    return false;
}

MTS_NAMESPACE_END