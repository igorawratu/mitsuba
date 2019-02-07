#include "common.h"

#include "arpaca.hpp"

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

                /*if(its.isEmitter()){
                    accumulator = its.Le(-ray.d);
                }
                else*/{
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

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf> partialSvd(const Eigen::MatrixXf& mat, 
    const std::uint32_t num_singular_values){
    assert(num_singular_values < std::min(mat.rows(), mat.cols()));

    std::uint32_t dim = std::max(mat.rows(), mat.cols());
    Eigen::MatrixXf jordan_wielandt = Eigen::MatrixXf::Zero(dim * 2, dim * 2);
    jordan_wielandt.block(0, dim, mat.rows(), mat.cols()) = mat;
    jordan_wielandt.block(dim, 0, mat.cols(), mat.rows()) = mat.transpose();

    const arpaca::EigenvalueType type = arpaca::ALGEBRAIC_LARGEST;
    arpaca::SymmetricEigenSolver<float> solver = arpaca::Solve(jordan_wielandt, num_singular_values, type);
    const Eigen::MatrixXf& eigenvectors = solver.eigenvectors();
    const Eigen::VectorXf& eigenvalues = solver.eigenvalues();

    Eigen::MatrixXf singular_values = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf u = Eigen::MatrixXf::Zero(mat.rows(), mat.rows());
    Eigen::MatrixXf v = Eigen::MatrixXf::Zero(mat.cols(), mat.cols());

    float scale = sqrt(2.f);

    for(int i = 0; i < eigenvalues.size(); ++i){
        singular_values(i, i) = eigenvalues(i);
        u.col(i) = scale * eigenvectors.block(0, i, mat.rows(), 1);
        v.col(i) = scale * eigenvectors.block(dim, i, mat.cols(), 1);
    }

    return std::make_tuple(u, singular_values, v);
}

Eigen::MatrixXf softThreshRank(const Eigen::MatrixXf& mat, float theta, const std::uint32_t step_size){
    std::uint32_t max_rank = std::min(mat.rows(), mat.cols());
    std::uint32_t curr_step_size = step_size;

    std::uint32_t dim = std::max(mat.rows(), mat.cols());
    Eigen::MatrixXf jordan_wielandt = Eigen::MatrixXf::Zero(dim * 2, dim * 2);
    jordan_wielandt.block(0, dim, mat.rows(), mat.cols()) = mat;
    jordan_wielandt.block(dim, 0, mat.cols(), mat.rows()) = mat.transpose();

    arpaca::SymmetricEigenSolver<float> solver;
    while(true){
        solver = arpaca::Solve(jordan_wielandt, curr_step_size, arpaca::ALGEBRAIC_LARGEST);
        if(curr_step_size == max_rank || solver.eigenvalues()(0) < theta){
            break;
        }
        curr_step_size = std::min(max_rank, curr_step_size + step_size);
    }

    const Eigen::MatrixXf& eigenvectors = solver.eigenvectors();
    const Eigen::VectorXf& eigenvalues = solver.eigenvalues();

    Eigen::MatrixXf singular_values = Eigen::MatrixXf::Zero(mat.rows(), mat.cols());
    Eigen::MatrixXf u = Eigen::MatrixXf::Zero(mat.rows(), mat.rows());
    Eigen::MatrixXf v = Eigen::MatrixXf::Zero(mat.cols(), mat.cols());

    float scale = sqrt(2.f);
    for(int i = 0; i < eigenvalues.size(); ++i){
        singular_values(i, i) = std::max(0.f, eigenvalues(i) - theta);
        u.col(i) = scale * eigenvectors.block(0, i, mat.rows(), 1);
        v.col(i) = scale * eigenvectors.block(dim, i, mat.cols(), 1);
    }
    
    return u * singular_values * v.transpose();
}

MTS_NAMESPACE_END