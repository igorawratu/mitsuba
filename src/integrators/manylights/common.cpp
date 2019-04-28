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
                            if(vpls[i].type == EDirectionalEmitterVPL){
                                continue;
                            }

                            if(abs((ray_origin - vpls[i].its.p).length() - t) > 0.0001f ){
                                continue;
                            }
                        }

                        BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
                        Spectrum albedo(0.f);
                        for(std::uint8_t i = 0; i < 10; ++i){
                            albedo += its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());
                        }
                        albedo /= 10.f;

                        float n_dot_ldir = std::max(0.f, dot(n, normalize(vpls[i].its.p - its.p)));

                        Spectrum c = (vpls[i].P * n_dot_ldir * albedo) / PI;

                        if(vpls[i].type == ESurfaceVPL){
                            float ln_dot_ldir = std::max(0.f, dot(vpls[i].its.shFrame.n, normalize(its.p - vpls[i].its.p)));
                            c *= ln_dot_ldir;
                        }

                        if(vpls[i].type != EDirectionalEmitterVPL){
                            float d = std::max((its.p - vpls[i].its.p).length(), min_dist);
                            float attenuation = 1.f / (d * d);
                            c *= attenuation;
                        }

                        accumulator += c;
                    }
                }
            }
            else{
                if(scene->hasEnvironmentEmitter()){
                    accumulator = scene->evalEnvironment(RayDifferential(ray));
                }
            }

            float r, g, b;
            accumulator.toLinearRGB(r, g, b);
            r = std::min(1.f, r);
            g = std::min(1.f, g);
            b = std::min(1.f, b);

            //convert from srgb to linear for correct scaling
            std::uint32_t offset = (x + y * film->getSize().x) * 3;
            float ir, ig, ib;
            Spectrum converter;
            converter.fromSRGB((float)image_buffer[offset] / 255.f, (float)image_buffer[offset + 1] / 255.f,
                (float)image_buffer[offset + 2] / 255.f);
            converter.toLinearRGB(ir, ig, ib);

            float error = fabs(ir - r) + fabs(ig - g) + fabs(ib - b);

            tot_error += error;

            error *= 20.f;

            float er, eg, eb;
            std::tie(er, eg, eb) = floatToRGB(std::min(1.f, error));

            image_buffer[offset] = er * 255 + 0.5f;
            image_buffer[offset + 1] = eg * 255 + 0.5f;
            image_buffer[offset + 2] = eb * 255 + 0.5f;
        }
    }

    return tot_error;
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

const float CONE_PROB = 1.f / (2.f * M_PI);

std::tuple<float, float, float> calculateMisWeight(const Vector3f& wi, const VPL& vpl, const Intersection& its){
    float light_prob;

    Vector3f light_wo = -wi;
    if(vpl.emitter != nullptr){
        DirectionSamplingRecord dsr(light_wo);
        light_prob = pdfDirection(dsr, vpl.psr);
    }
    else{
        BSDFSamplingRecord bsdfsr(vpl.its, light_wo);
        light_prob = vpl.its.getBSDF()->pdf(bsdfsr);
    }

    BSDFSamplingRecord bsdfsr(its, wi);
    float bsdf_prob = its.getBSDF()->pdf(bsdfsr);

    return std::make_tuple(CONE_PROB, bsdf_prob, light_prob);
}

Spectrum sampleCone(const Vector3f& wi, const VPL& vpl, const Intersection& its, float d, Sampler* sampler){
    float r = sampler->next1D() * vpl.radius;
    float theta = sampler->next1D() * 2.f * M_PI;
    Vector3f cone_sample(cos(theta) * sqrt(r), sin(theta) * sqrt(r), d);
    Frame cone_frame(wi);
    cone_sample = normalize(cone_frame.toWorld(cone_sample));
    Vector3f sample_wi = normalize(cone_sample);

    BSDFSamplingRecord bsdf_sample_record(its, its.toLocal(sample_wi));
    Spectrum curr_col = vpl.P * bsdf->eval(bsdf_sample_record);

    Spectrum weight(0.f);

    if(vpl.emitter != nullptr){
        DirectionSamplingRecord dir;
        dir.d = -sample_wi;
        curr_col *= vpl.emitter->evalDirection(dir, vpl.psr);
    }
    else if(vpl.its.getBSDF() != nullptr){
        BSDFSamplingRecord light_sample_record(vpl.its, vpl.its.toLocal(-sample_wi));
        curr_col *= vpl.its.getBSDF()->eval(light_sample_record);
    }

    float cone_prob, brdf_prob, light_prop;
    std::tie(cone_prob, brdf_prob, light_prob) = calculateMisWeight(sample_wi, vpl, its);
    curr_col *= cone_prob / (cone_prob + brdf_prob + light_prob);

    return curr_col;
}

Spectrum sampleBsdf(const VPL& vpl, const Intersection& its, Sampler* sampler){
    BSDFSamplingRecord bsdf_sample_record(its, sampler);
    bsdf_sample_record.typeMask = BSDF::EReflection;
    Spectrum curr_col = its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());

    if(raySphereIntersect(its.p, bsdf_sample_record.wi, vpl.its.p, vpl.radius)){
        curr_col *= vpl.P;

        if(vpl.emitter != nullptr){
            DirectionSamplingRecord dir;
            dir.d = -bsdf_sample_record.wi;
            curr_col *= vpl.emitter->evalDirection(dir, vpl.psr);
        }
        else if(vpl.its.getBSDF() != nullptr){
            BSDFSamplingRecord light_sample_record(vpl.its, vpl.its.toLocal(-bsdf_sample_record.wi));
            curr_col *= vpl.its.getBSDF()->eval(light_sample_record);
        }

        float cone_prob, brdf_prob, light_prop;
        std::tie(cone_prob, brdf_prob, light_prob) = calculateMisWeight(bsdf_sample_record.wi, vpl, its);
        curr_col *= brdf_prob / (cone_prob + brdf_prob + light_prob);

        return curr_col
    }

    return Spectrum(0.f);
}

Spectrum sampleLight(const VPL& vpl, const Intersection& its, Sampler* sampler){
    Vector3f wi;
    Spectrum curr_col;
    if(vpl.emitter != nullptr){
        DirectionSamplingRecord dsr;
        curr_col = vpl.emitter->sampleDirection(dsr, vpl.psr, sampler->next2D());
        wi = -dsr.d;
    }
    else{
        BSDFSamplingRecord bsdf_sample_record(vpl.its, sampler);
        bsdf_sample_record.typeMask = BSDF::EReflection;
        curr_col = vpl.its.getBSDF()->sample(bsdf_sample_record, sampler->next2D()) * vpl.P;
        wi = -vpl.its.toWorld(bsdf_sample_record.wo);
    }

    if(raySphereIntersect(its.p, wi, vpl.its.p, vpl.radius)){
        BSDFSamplingRecord bsdf_sample_record(its, its.toLocal(wi));
        curr_col *= bsdf->eval(bsdf_sample_record);

        float cone_prob, brdf_prob, light_prop;
        std::tie(cone_prob, brdf_prob, light_prob) = calculateMisWeight(wi, vpl, its);
        curr_col *= brdf_prob / (cone_prob + brdf_prob + light_prob);

        return curr_col
    }

    return Spectrum(0.f);

}

//-------taken from http://miloshasan.net/VirtualSphericalLights/vsl.fx--------
float circleSegmentArea(float h, float r)
{
    h /= r;
    float  A = (h * sqrt(1-h*h) + asin(h) + M_PI/2);
    return A*r*r;
}

float distFromPlane(Point p, Vector3f np, Point q, Vector3f nq)
{
    float dpq = dot(np,nq);
    if(dpq*dpq>0.99f) return 1e20f;

    float l = dot(p-q, nq);

    return sqrt( l*l / (1.f-dpq*dpq) );
}

float effectiveLightArea(float lightR, float3 lightP, float3 lightN, float3 P, float3 N)
{
    float h = distFromPlane(lightP, lightN, P, N);
    return (h<lightR) ? circleSegmentArea(h,lightR) : M_PI*lightR*lightR;
}
//-----------------------------------------------------------------------------

Spectrum sample(Scene* scene, Sampler* sampler, Intersection& its, const Ray& initial_ray, const VPL& vpl, float min_dist, 
    bool check_occlusion, std::uint32_t max_specular_bounces, bool perform_ray_intersection, bool& intersected,
    bool show_emitter, bool vsl){

    Ray ray = initial_ray;

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
            bsdf_sample_record.typeMask = BSDF::EReflection;
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
                ((its.p - vpl.its.p).length() - t) > 1e-6f * min_dist){
                return Spectrum(0.f);
            }
        }
    }

    RayDifferential rd(ray);

    const BSDF *bsdf = its.getBSDF(rd);

    Spectrum c(0.f);

    if(vsl && vpl.type == ESurfaceVPL){
        float central_disc_area = effectiveLightArea(vpl.radius, vpl.its.p, vpl.its.shFrame.n, its.p, its.geoFrame.n);
        float d = (vpl.its.p - its.p).length();
        float solid_angle = 2 * M_PI * (1.f - cos(asin(std::min(vpl.radius / d, 1.f))));
        float solid_angle_ratio = solid_angle / (2.f * M_PI);
        std::uint32_t num_samples = std::max(1u, std::uint32_t(sqrt(solid_angle_ratio) * 1.f));
        Spectrum total(0.f);

        for(std::uint32_t i = 0; i < num_samples; ++i){
            total += sampleCone(wi, vpl, its, d, sampler);
            total += sampleBsdf(vpl, its, sampler);
            total += sampleLight(vpl, its, sampler);
        }

        c = solid_angle_ratio * total / (float(num_samples) * central_disc_area);
    }
    else{
        //only care about non-specular surfaces for now, just return specular reflectance
        if(!(bsdf->getType() & BSDF::ESmooth)){
            BSDFSamplingRecord bsdf_sample_record(its, sampler);
            c = vpl.P * bsdf->sample(bsdf_sample_record, sampler->next2D()) * dot(its.shFrame.n, wi) / PI;
        }
        else{
            BSDFSamplingRecord bsdf_sample_record(its, its.toLocal(wi));
            c = vpl.P * bsdf->eval(bsdf_sample_record);
        }

        if(vpl.type != EDirectionalEmitterVPL){
            float d = std::max((its.p - vpl.its.p).length(), min_dist);
            float attenuation = 1.f / (d * d);
            c *= attenuation;
        }

        if(vpl.type == ESurfaceVPL){
            if(vpl.emitter != nullptr){
                DirectionSamplingRecord dir;
                dir.d = -wi;
                c *= vpl.emitter->evalDirection(dir, vpl.psr);
            }
            else if(vpl.its.getBSDF() != nullptr){
                BSDFSamplingRecord bsdf_sample_record(vpl.its, vpl.its.toLocal(-wi));
                c *= vpl.its.getBSDF()->eval(bsdf_sample_record);
            }
            //fallback to diffuse if no bsdf or emitter found
            else{
                c *= Spectrum(std::max(0.f, dot(vpl.its.shFrame.n, -wi)));
            }
        }
    }

    return c;
}

MTS_NAMESPACE_END