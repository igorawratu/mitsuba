#include "common.h"

#include <tuple>

MTS_NAMESPACE_BEGIN

std::tuple<float, float, float> floatToRGB(float v){
    float r = 0.f, g = 0.f, b = 0.f;

    if(v < 0.25f){
        r = 0.f;
        b = 1.f;
        g = v * 4.f;
    }
    else if(v >= 0.25f && v < 0.5f){
        r = 0.f;
        b = 1.f - (v - 0.25f) * 4.f;
        g = 1.f;
    }
    else if(v >= 0.5f && v < 0.75f){
        r = (v - 0.5f) * 4.f;
        b = 0.f;
        g = 1.f;
    }
    else if(v >= 0.75f){
        r = 1.f;
        b = 0.f;
        g = 1.f - (v - 0.75f) * 4.f;
    }

    return std::make_tuple(r, g, b);
}

bool raySphereIntersect(Point ray_origin, Vector3f ray_d, Point sphere_center, float sphere_radius){
    float d_sq = (ray_origin - sphere_center).lengthSquared();
    float sr2 = sphere_radius * sphere_radius;
    if(d_sq < sr2){
        return true;
    }

    Vector3f v(sphere_center - ray_origin);
    Vector3f projection = dot(v, ray_d) / ray_d.lengthSquared() * ray_d;

    return dot(projection, ray_d) >= 0.f && projection.lengthSquared() >= sr2;
}

std::tuple<float, float, float> calculateMisWeight(const Vector3f& wi, const VPL& vpl, const Intersection& its,
    float solid_angle){
    
    float light_prob;

    Vector3f light_wo = -wi;
    if(vpl.emitter != nullptr){
        DirectionSamplingRecord dsr(light_wo);
        light_prob = vpl.emitter->pdfDirection(dsr, vpl.psr);
    }
    else{
        BSDFSamplingRecord bsdfsr(vpl.its, light_wo);
        light_prob = vpl.its.getBSDF()->pdf(bsdfsr);
    }

    BSDFSamplingRecord bsdfsr(its, its.toLocal(wi));
    bsdfsr.typeMask = BSDF::EReflection;
    float bsdf_prob = its.getBSDF()->pdf(bsdfsr);

    return std::make_tuple(1.f / solid_angle, bsdf_prob, light_prob);
}

Spectrum sampleCone(const Vector3f& wi, const VPL& vpl, const Intersection& its, float d, Sampler* sampler,
    float solid_angle){
    float r = sampler->next1D() * vpl.radius;
    float theta = sampler->next1D() * 2.f * M_PI;
    Vector3f cone_sample(cos(theta) * sqrt(r), sin(theta) * sqrt(r), d);
    Frame cone_frame(wi);
    cone_sample = normalize(cone_frame.toWorld(cone_sample));
    Vector3f sample_wi = normalize(cone_sample);

    BSDFSamplingRecord bsdf_sample_record(its, its.toLocal(sample_wi));
    bsdf_sample_record.typeMask = BSDF::EReflection;
    Spectrum curr_col = vpl.P * its.getBSDF()->eval(bsdf_sample_record);

    Spectrum weight(0.f);

    if(vpl.emitter != nullptr){
        DirectionSamplingRecord dir(-sample_wi);
        curr_col *= vpl.emitter->evalDirection(dir, vpl.psr);
    }
    else if(vpl.its.getBSDF() != nullptr){
        BSDFSamplingRecord light_sample_record(vpl.its, vpl.its.toLocal(-sample_wi));
        curr_col *= vpl.its.getBSDF()->eval(light_sample_record);
    }
    else{
        curr_col *= Spectrum(std::max(0.f, dot(vpl.its.shFrame.n, -sample_wi))) / PI;
    }

    float cone_prob, brdf_prob, light_prob;
    std::tie(cone_prob, brdf_prob, light_prob) = calculateMisWeight(sample_wi, vpl, its, solid_angle);
    curr_col /= cone_prob + brdf_prob;

    return curr_col * dot(sample_wi, normalize(vpl.its.p - its.p));
}

Spectrum sampleBsdf(const VPL& vpl, const Intersection& its, Sampler* sampler, float solid_angle, float cos_theta){
    BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
    bsdf_sample_record.typeMask = BSDF::EReflection;
    Float pdf;
    Spectrum curr_col = its.getBSDF()->sample(bsdf_sample_record, pdf, sampler->next2D());
    curr_col *= pdf;
    Vector3f wo = its.toWorld(bsdf_sample_record.wo);
    
    float wi_dot_vp = std::max(0.f, dot(wo, normalize(vpl.its.p - its.p)));
    if(wi_dot_vp > cos_theta){
        if(vpl.emitter != nullptr){
            DirectionSamplingRecord dir(-wo);
            curr_col *= vpl.P * vpl.emitter->evalDirection(dir, vpl.psr);
        }
        else if(vpl.its.getBSDF() != nullptr){
            BSDFSamplingRecord light_sample_record(vpl.its, vpl.its.toLocal(-wo));
            curr_col *= vpl.P * vpl.its.getBSDF()->eval(light_sample_record);
        }
        else{
            curr_col *= vpl.P * Spectrum(std::max(0.f, dot(vpl.its.shFrame.n, -wo))) / PI;
        }

        float cone_prob, brdf_prob, light_prob;
        std::tie(cone_prob, brdf_prob, light_prob) = calculateMisWeight(wo, vpl, its, solid_angle);
        if(brdf_prob > std::numeric_limits<float>::epsilon()){
            curr_col /= brdf_prob + cone_prob;
        }
        else curr_col = Spectrum(0.f);

        return curr_col * wi_dot_vp;
    }

    return Spectrum(0.f);
}

Spectrum sample(Scene* scene, Sampler* sampler, Intersection& its, const Ray& initial_ray, const VPL& vpl, float min_dist, 
    bool check_occlusion, std::uint32_t max_specular_bounces, bool perform_ray_intersection, bool& intersected,
    bool show_emitter, bool vsl, std::uint32_t& samples_taken){

    Ray ray = initial_ray;

    samples_taken = 1;

    if(perform_ray_intersection){
        std::uint32_t num_bounces = 0;

        while(true){
            intersected = scene->rayIntersect(ray, its);
            if(!intersected){
                break;
            }

            if(its.getBSDF()->getType() & BSDF::ESmooth || its.isEmitter()){
                break;
            }

            if(++num_bounces > max_specular_bounces){
                break;
            }

            BSDFSamplingRecord bsdf_sample_record(its, sampler);
            bsdf_sample_record.typeMask = its.getBSDF()->isDielectric() ? BSDF::EDeltaTransmission | BSDF::ENull : 
                BSDF::EDeltaReflection;
            its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());

            ray = Ray(its.p, bsdf_sample_record.its.toWorld(bsdf_sample_record.wo), ray.time);
        }
    }

    if(!intersected){
        if(scene->hasEnvironmentEmitter()){
            return scene->evalEnvironment(RayDifferential(ray));
        }
        else return Spectrum(0.f);
    }

    if(show_emitter && its.isEmitter()){
        return its.Le(-ray.d);
    }
    
    Vector3f wi = vpl.type == EDirectionalEmitterVPL ? -Vector3f(vpl.its.shFrame.n) : 
        normalize(vpl.its.p - its.p);

    if(check_occlusion){
        Ray shadow_ray(its.p, wi, 0.f);
        Float t;
        ConstShapePtr shape;
        Normal norm;
        Point2 uv;

        if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
            if(vpl.type == EDirectionalEmitterVPL || 
                (its.p - vpl.its.p).length() - t > 1e-5f * min_dist){
                return Spectrum(0.f);
            }
        }
    }

    RayDifferential rd(ray);

    const BSDF *bsdf = its.getBSDF(rd);

    Spectrum c(0.f);

    if(vsl && vpl.type == ESurfaceVPL){
        float central_disc_area = M_PI * vpl.radius * vpl.radius;
        float d = (vpl.its.p - its.p).length();
        
        float hypot_length = sqrt(d * d + vpl.radius * vpl.radius);
        float cos_theta = d / hypot_length;
        //float cos_theta = cos(asin(std::min(vpl.radius / d, 1.f)));

        float solid_angle = 2.f * M_PI * (1.f - cos_theta);
        std::uint32_t num_samples = std::max(1u, std::uint32_t(sqrt(1.f - cos_theta) * 100.f));
        samples_taken = num_samples;
        Spectrum total(0.f);

        for(std::uint32_t i = 0; i < num_samples; ++i){
            total += sampleCone(wi, vpl, its, d, sampler, solid_angle);
            total += sampleBsdf(vpl, its, sampler, solid_angle, cos_theta);
        }

        c = total / (float(num_samples) * central_disc_area);
    }
    else{
        //only care about non-specular surfaces for now, just return specular reflectance
        if(!(bsdf->getType() & BSDF::ESmooth)){
            BSDFSamplingRecord bsdf_sample_record(its, sampler);
            c = bsdf->sample(bsdf_sample_record, sampler->next2D()) * dot(its.shFrame.n, wi) / PI;
        }
        else{
            BSDFSamplingRecord bsdf_sample_record(its, its.toLocal(wi));
            c = bsdf->eval(bsdf_sample_record);
        }

        c *= vpl.P;

        if(vpl.type != EDirectionalEmitterVPL){
            float d = std::max((its.p - vpl.its.p).length(), min_dist);
            float attenuation = 1.f / (d * d);
            c *= attenuation;
        }

        if(vpl.type == ESurfaceVPL){
            //Vector3f n = vpl.its.wi.z > 0.f ? Vector3f(vpl.its.shFrame.n) : Vector3f(-vpl.its.shFrame.n);

            if(vpl.emitter != nullptr){
                DirectionSamplingRecord dir(-wi);
                //c *= vpl.emitter->evalDirection(dir, vpl.psr);
                c *= Spectrum(std::min(1.f, std::max(0.f, dot(vpl.its.shFrame.n, -wi)))) / PI;
            }
            else if(vpl.its.getBSDF() != nullptr){
                BSDFSamplingRecord bsdf_sample_record(vpl.its, vpl.its.toLocal(-wi));
                //c *= vpl.its.getBSDF()->eval(bsdf_sample_record);
                c *= Spectrum(std::min(1.f, std::max(0.f, dot(vpl.its.shFrame.n, -wi)))) / PI;
            }
            //fallback to diffuse if no bsdf or emitter found
            else{
                c *= Spectrum(std::min(1.f, std::max(0.f, dot(vpl.its.shFrame.n, -wi)))) / PI;
            }
        }
    }
    return c;
}

bool sampleVisibility(Scene* scene, const Intersection& its, const VPL& vpl, float min_dist){
    Point ray_origin = its.p;
    Vector dir = vpl.type == EDirectionalEmitterVPL ? -Vector3f(vpl.its.shFrame.n) : 
        normalize(vpl.its.p - ray_origin);

    Ray shadow_ray(ray_origin, dir, 0.f);

    bool visible = true;

    Float t;
    ConstShapePtr shape;
    Normal norm;
    Point2 uv;

    if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
        if(std::abs((ray_origin - vpl.its.p).length() - t) > 1e-5f * min_dist){
            visible = false;
        }
    }

    return visible;
}

std::uint64_t upperPo2(std::uint64_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

MTS_NAMESPACE_END