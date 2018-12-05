#include "matrixseparationrenderer.h"

#include <vector>
#include <limits>
#include <utility>
#include <memory>

#include <mitsuba/core/plugin.h>

#include "definitions.h"

MTS_NAMESPACE_BEGIN

struct RowSample{
    RowSample(const Point3f& p, const Normal& n, std::uint32_t x, std::uint32_t y, std::uint32_t cols, bool in_scene, Intersection intersection) : 
        position(p), normal(n), image_x(x), image_y(y), col_samples(cols), intersected_scene(in_scene), its(intersection){
    }

    Point3f position;
    Normal normal;
    std::uint32_t image_x, image_y;
    std::vector<Spectrum> col_samples;
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
                std::fill(curr_sample.col_samples.begin(), 
                    curr_sample.col_samples.end(), Spectrum(0.f));
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



bool MatrixSeparationRenderer::Render(Scene* scene){
    return false;
}

MTS_NAMESPACE_END