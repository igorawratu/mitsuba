#include "matrixreconstruction.h"

#include <mitsuba/core/plugin.h>
#include <chrono>
#include <algorithm>
#include <set>

MTS_NAMESPACE_BEGIN

MatrixReconstructionRenderer::MatrixReconstructionRenderer() : cancel_(false){

}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(const MatrixReconstructionRenderer& other) : cancel_(other.cancel_){

}

MatrixReconstructionRenderer::MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other) : cancel_(other.cancel_){

}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (const MatrixReconstructionRenderer& other){
    if(this != &other){
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer& MatrixReconstructionRenderer::operator = (MatrixReconstructionRenderer&& other){
    if(this != &other){
        cancel_ = other.cancel_;
    }
    return *this;
}

MatrixReconstructionRenderer::~MatrixReconstructionRenderer(){

}

std::vector<std::pair<std::uint32_t, std::uint32_t>> generateComputableIndices(const std::pair<std::uint32_t, std::uint32_t>& bucket_size, 
    const Vector2i& image_size, std::uint32_t num_lights, std::uint32_t light_samples){

    std::uint32_t x_buckets = image_size.x / bucket_size.first + image_size.x % (bucket_size.first > 0 ? 1 : 0);
    std::uint32_t y_buckets = image_size.y / bucket_size.second + image_size.y % (bucket_size.second > 0 ? 1 : 0);
    light_samples = std::min(light_samples, num_lights);

    std::vector<std::pair<std::uint32_t, std::uint32_t>> indices_to_compute;
    indices_to_compute.reserve(x_buckets * y_buckets * light_samples);

    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);

    for(std::uint32_t i = 0; i < y_buckets; ++i){
        for(std::uint32_t j = 0; j < x_buckets; ++j){
            std::uint32_t x_start = j * bucket_size.first;
            std::uint32_t y_start = i * bucket_size.second;

            std::uniform_int_distribution<std::uint32_t> gen_y(0, std::min(bucket_size.second, image_size.y - y_start) - 1);
            std::uniform_int_distribution<std::uint32_t> gen_x(0, std::min(bucket_size.first, image_size.x - x_start) - 1);
            std::uniform_int_distribution<std::uint32_t> gen_light(0, num_lights - 1);

            std::uint32_t max_tries = num_lights * 10;
            std::uint32_t current_tries = 0;
            std::set<std::pair<std::uint32_t, std::uint32_t>> light_indices;
            
            while(current_tries++ < max_tries && light_indices.size() < light_samples){
                light_indices.insert(std::make_pair(
                    (gen_x(rng) + x_start) + (gen_y(rng) + y_start) * image_size.x, gen_light(rng)));
            }

            for(auto iter = light_indices.begin(); iter != light_indices.end(); ++iter){
                indices_to_compute.emplace_back(*iter);
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
        const VPL& vpl = vpls[indices[i].second];

        Ray ray;

        Point2 sample_position(x + 0.5f, y + 0.5f);

        sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            continue;
        }

        Normal n = its.geoFrame.n;

        Point ray_origin = its.p;
        Ray shadow_ray(ray_origin, normalize(vpl.its.p - ray_origin), ray.time);

        Float t;
        ConstShapePtr shape;
        Normal norm;
        Point2 uv;
        if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
            if(abs((ray_origin - vpl.its.p).length() - t) > 0.0001f ){
                continue;
            }
        }

        BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
        bsdf_sample_record.wi = its.toLocal(normalize(vpl.its.p - its.p));
        bsdf_sample_record.wo = its.toLocal(n);

        Spectrum albedo = its.getBSDF()->eval(bsdf_sample_record);

        float d = std::max((its.p - vpl.its.p).length(), min_dist);
        float attenuation = 1.f / (d * d);

        float n_dot_ldir = std::max(0.f, dot(normalize(n), normalize(vpl.its.p - its.p)));
        float ln_dot_ldir = std::max(0.f, dot(normalize(vpl.its.shFrame.n), normalize(its.p - vpl.its.p)));

        Spectrum lightContribution = (vpl.P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
        Float r, g, b;
        lightContribution.toSRGB(r, g, b);
        matrix(indices[i].first * 3, indices[i].second) = r;
        matrix(indices[i].first * 3 + 1, indices[i].second) = g;
        matrix(indices[i].first * 3 + 2, indices[i].second) = b;
    }
}

void copyMatrixToBuffer(std::uint8_t* output_image, Eigen::MatrixXf& light_matrix, Vector2i image_size){
    for(std::uint32_t i = 0; i < light_matrix.rows() / 3; ++i){
        float r = 0, g = 0, b = 0;
        for(int j = 0; j < light_matrix.cols(); ++j){
            r += light_matrix(i*3, j);
            g += light_matrix(i*3 + 1, j);
            b += light_matrix(i*3 + 2, j);
        }

        output_image[i*3] = std::min(1.f, r) * 255 + 0.5f;
        output_image[i*3 + 1] = std::min(1.f, g) * 255 + 0.5f;
        output_image[i*3 + 2] = std::min(1.f, b) * 255 + 0.5f;
    }
}

bool MatrixReconstructionRenderer::Render(Scene* scene, const std::vector<VPL>& vpls,
    const std::pair<std::uint32_t, std::uint32_t>& bucket_size, const std::uint32_t& light_samples, float min_dist,
    std::uint8_t* output_image){

    if(scene == nullptr || vpls.size() == 0 || output_image == nullptr){
        return true;
    }

    cancel_ = false;

    ref<Sensor> sensor = scene->getSensor();
	ref<Film> film = sensor->getFilm();

    auto size = film->getSize();
    if(size.x == 0 || size.y == 0){
        return true;
    }

    Eigen::MatrixXf lighting_matrix = Eigen::MatrixXf::Zero(size.x * size.y * 3, vpls.size());

    auto indices_to_compute = generateComputableIndices(bucket_size, size, vpls.size(), light_samples);
    calculateSparseSamples(scene, vpls, lighting_matrix, indices_to_compute, size, min_dist);
    
    //nuclear norm minimization here

    copyMatrixToBuffer(output_image, lighting_matrix, size);

    return cancel_;
}

MTS_NAMESPACE_END