#include "matrixseparationrenderer.h"

#include <vector>
#include <limits>
#include <utility>
#include <memory>
#include <random>

#include <mitsuba/core/plugin.h>

#include "definitions.h"

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

            Normal n = its.geoFrame.n;

            for (std::uint32_t i = 0; i < vpls.size(); ++i) {
                Spectrum albedo(0.f);

                BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
                    bsdf_sample_record.wi = its.toLocal(normalize(vpls[i].its.p - its.p));
                    bsdf_sample_record.wo = its.toLocal(n);

                albedo = its.getBSDF()->eval(bsdf_sample_record);

                //only dealing with emitter and surface VPLs curently.
                if (vpls[i].type != EPointEmitterVPL && vpls[i].type != ESurfaceVPL){
                    continue;
                }

                float d = std::max((its.p - vpls[i].its.p).length(), min_dist);
                float attenuation = 1.f / (d * d);

                float n_dot_ldir = std::max(0.f, dot(normalize(n), normalize(vpls[i].its.p - its.p)));
                float ln_dot_ldir = std::max(0.f, dot(normalize(vpls[i].its.shFrame.n), normalize(its.p - vpls[i].its.p)));

                curr_sample.col_samples[i] = (vpls[i].P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
            }
        }
    }

    splitKDTree(kdt_root.get(), size_threshold, min_dist);

    return kdt_root;
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
        
        //PERFORM PREDICTIONS HERE!!!!111!1!11!!one!eleven

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
            Ray shadow_ray(ray_origin, normalize(vpls[i].its.p - ray_origin), 0.f);

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            samples[i] = VISIBLE;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if(abs((ray_origin - vpls[i].its.p).length() - t) > 0.0001f ){
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

bool MatrixSeparationRenderer::Render(Scene* scene){
    return false;
}

MTS_NAMESPACE_END