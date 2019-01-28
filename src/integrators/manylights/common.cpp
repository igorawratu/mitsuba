#include "common.h"

#include <tuple>

MTS_NAMESPACE_BEGIN

float calculateError(Scene* scene, const std::vector<VPL>& vpls, float min_dist, std::uint8_t* image_buffer){
    float tot_error = 0.f;
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

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

            Point2i curr_pixel(x, y);
            Spectrum accumulator(0.f);

            if (scene->rayIntersect(ray, its)) {
                Normal n = its.geoFrame.n;

                Spectrum albedo(0.f);

                if(its.isEmitter()){
                    accumulator = its.Le(-ray.d);
                }
                else{
                    for (std::uint32_t i = 0; i < vpls.size(); ++i) {
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

                        albedo = its.getBSDF()->eval(bsdf_sample_record);

                        //only dealing with emitter and surface VPLs curently.
                        if (vpls[i].type != EPointEmitterVPL && vpls[i].type != ESurfaceVPL){
                            continue;
                        }

                        float d = std::max((its.p - vpls[i].its.p).length(), min_dist);
                        float attenuation = 1.f / (d * d);

                        float n_dot_ldir = std::max(0.f, dot(normalize(n), normalize(vpls[i].its.p - its.p)));
                        float ln_dot_ldir = std::max(0.f, dot(normalize(vpls[i].its.shFrame.n), normalize(its.p - vpls[i].its.p)));

                        accumulator += (vpls[i].P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
                    }
                }
            }

            float r, g, b;
            accumulator.toSRGB(r, g, b);

            std::uint32_t offset = (x + y * film->getSize().x) * 3;

            float error = fabs((float)image_buffer[offset] / 255.f - r)
                + fabs((float)image_buffer[offset + 1] / 255.f - g)
                + fabs((float)image_buffer[offset + 2] / 255.f - b);

            tot_error += error;

            error = error / 3.f * 20.f;

            image_buffer[offset] = std::min(1.f, error) * 255 + 0.5f;
            image_buffer[offset + 1] = 0;
            image_buffer[offset + 2] = (1.f - std::min(1.f, error)) * 255 + 0.5f;
        }
    }

    return tot_error;
}

MTS_NAMESPACE_END