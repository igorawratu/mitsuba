#include "matrixreconstruction.h"

#include <mitsuba/core/plugin.h>
#include <chrono>
#include <algorithm>
#include <set>

MTS_NAMESPACE_BEGIN

MatrixReconstructionRenderer::MatrixReconstructionRenderer(){

}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(const MatrixReconstructionRenderer& other){

}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other){

}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (const MatrixReconstructionRenderer& other){

}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (MatrixReconstructionRenderer&& other){

}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::vector<std::pair<std::uint32_t, std::uint32_t>> generateComputableIndices(const std::pair<std::uint32_t, std::uint32_t>& bucket_size, 
    const Vector2i& image_size, std::uint32_t num_lights, std::uint32_t light_samples){

    std::uint32_t x_buckets = image_size.x / bucket_size.first + image_size.x % bucket_size.first > 0 ? 1 : 0;
    std::uint32_t y_buckets = image_size.y / bucket_size.second + image_size.y % bucket_size.second > 0 ? 1 : 0;
    light_samples = std::min(light_samples, num_lights);

    std::vector<std::pair<std::uint32_t, std::uint32_t>> indices_to_compute;
    indices_to_compute.reserve(x_buckets * y_buckets * light_samples);

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    for(std::uint32_t i = 0; i < y_buckets; ++i){
        for(std::uint32_t j = 0; j < x_buckets; ++j){
            std::uint32_t x_start = j * bucket_size.first;
            std::uint32_t y_start = i * bucket_size.second;

            auto gen_y = std::bind(std::uniform_int_distribution<std::uint32_t>(0, 
                            std::min(bucket_size.second, image_size.y - y_start) - 1), rng);
            auto gen_x = std::bind(std::uniform_int_distribution<std::uint32_t>(0, 
                            std::min(bucket_size.first, image_size.x - x_start) - 1), rng);
            auto gen_light = std::bind(std::uniform_int_distribution<std::uint32_t>(0, num_lights - 1), rng);


            std::uint32_t max_tries = num_lights * 10;
            std::uint32_t current_tries = 0;
            std::set<std::uint32_t> light_indices;
            
            while(current_tries++ < max_tries && light_indices.size() < light_samples){
                light_indices.insert(gen_light());
            }

            for(auto iter = light_indices.begin(); iter != light_indices.end(); ++iter){
                indices_to_compute.emplace_back(gen_x() + gen_y() * image_size.x, *iter);
            }
        }
    }

    return indices_to_compute;
}

const double PI = 3.14159265359;

void calculateSparseSamples(Scene* scene, const std::vector<VPL>& vpls, Eigen::MatrixXf& matrix,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& indices, const Vector2i& size, float min_dist){
    
    ref<Sensor> sensor = scene->getSensor();
	ref<Film> film = sensor->getFilm();

    Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

    //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
    Point2 aperture_sample(0.5f, 0.5f);
    Float time_sample(0.5f);

    for(size_t i = 0; i < indices.size(); ++i){
        std::uint32_t x = indices[i].first % size.x;
        std::uint32_t y = indices[i].first / size.x;

        Ray ray;

        Point2 sample_position(x + 0.5f, y + 0.5f);

        sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

        Intersection its;

        if (!scene->rayIntersect(ray, its)) {
            continue;
        }

        Normal n = its.geoFrame.n;

        Point ray_origin = its.p;
        Ray shadow_ray(ray_origin, normalize(vpls[i].its.p - ray_origin), ray.time);

        Float t;
        ConstShapePtr shape;
        Normal norm;
        Point2 uv;

        if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
            if(abs((ray_origin - vpls[i].its.p).length() - t) > 0.0001f ){
                continue;
            }
        }

        BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
        bsdf_sample_record.wi = its.toLocal(normalize(vpls[i].its.p - its.p));
        bsdf_sample_record.wo = its.toLocal(n);

        Spectrum albedo = its.getBSDF()->eval(bsdf_sample_record);

        float d = std::max((its.p - vpls[i].its.p).length(), min_dist);
        float attenuation = 1.f / (d * d);

        float n_dot_ldir = std::max(0.f, dot(normalize(n), normalize(vpls[i].its.p - its.p)));
        float ln_dot_ldir = std::max(0.f, dot(normalize(vpls[i].its.shFrame.n), normalize(its.p - vpls[i].its.p)));

        Spectrum lightContribution = (vpls[i].P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
        Float r, g, b;
        lightContribution.toLinearRGB(r, g, b);
        matrix(indices[i].first * 3, indices[i].second) = r;
        matrix(indices[i].first * 3 + 1, indices[i].second) = r;
        matrix(indices[i].first * 3 + 2, indices[i].second) = r;
    }
}

void MatrixReconstructionRenderer::Render(Scene* scene, const std::vector<VPL>& vpls,
    const std::pair<std::uint32_t, std::uint32_t>& bucket_size, const std::uint32_t& light_samples, float min_dist){

    if(scene == nullptr || vpls.size() == 0){
        return;
    }

    ref<Sensor> sensor = scene->getSensor();
	ref<Film> film = sensor->getFilm();

    auto size = film->getSize();
    if(size.x == 0 || size.y == 0){
        return;
    }

    Eigen::MatrixXf lighting_matrix = Eigen::MatrixXf::Zero(size.x * size.y * 3, vpls.size());

    auto indices_to_compute = generateComputableIndices(bucket_size, size, vpls.size(), light_samples);
    calculateSparseSamples(scene, vpls, lighting_matrix, indices_to_compute, size, min_dist);
    
    
}

MTS_NAMESPACE_END