#include "lightclustererrenderer.h"

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include "definitions.h"
#include "common.h"
#include "hwshader.h"
#include "blockingqueue.hpp"
#include <unordered_map>
#include "filewriter.h"
#include <chrono>

MTS_NAMESPACE_BEGIN

LightClustererRenderer::LightClustererRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, float min_dist, bool vsl, bool hw) : 
    min_dist_(min_dist), clusterer_(std::move(clusterer)), vsl_(vsl), hw_(hw), cancel_(false){
}

LightClustererRenderer::LightClustererRenderer(LightClustererRenderer&& other) : 
    min_dist_(other.min_dist_), clusterer_(std::move(other.clusterer_)),
    vsl_(other.vsl_), hw_(other.hw_), cancel_(other.cancel_){
}

LightClustererRenderer& LightClustererRenderer::operator = (LightClustererRenderer&& other){
    if(this != &other){
        this->clusterer_ = std::move(other.clusterer_);
        this->vsl_ = other.vsl_;
        this->cancel_ = other.cancel_;
        this->hw_ = other.hw_;
        this->min_dist_ = other.min_dist_;
    }

    return *this;
}

LightClustererRenderer::~LightClustererRenderer(){
}

std::vector<HWBFPix> generateReceivers(Scene* scene, std::uint32_t spp, Sampler* sampler, Vector2i size){
    ref<Sensor> sensor = scene->getSensor();

    std::vector<HWBFPix> receivers;
    receivers.resize(size.y * size.x * spp);

    #pragma omp parallel for
    for (std::int32_t y = 0; y < size.y; ++y) {
        for (std::int32_t x = 0; x < size.x; ++x) {
            std::uint32_t cell_dim = sqrt(spp) + 0.5f;
            float cell_side_len = 1.f / cell_dim;

            for(std::uint32_t j = 0; j < spp; ++j){
                float x_jitter = 0.5f * cell_side_len;
                float y_jitter = 0.5f * cell_side_len;
                float x_off = (j % cell_dim) * cell_side_len;
                float y_off = (j / cell_dim) * cell_side_len;

                Point2 sample_position(x + x_off + x_jitter, y + y_off + y_jitter);
                
                //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
                Point2 aperture_sample(0.5f, 0.5f);
                Float time_sample(0.5f);

                HWBFPix curr;

                sensor->sampleRay(curr.ray, sample_position, aperture_sample, time_sample);
                curr.x = x;
                curr.y = y;
                curr.col = Spectrum(0.f);

                std::uint32_t num_bounces = 0;

                bool intersected_scene = false;

                while(true){
                    intersected_scene = scene->rayIntersect(curr.ray, curr.its);
                    if(!intersected_scene){
                        break;
                    }

                    if((curr.its.getBSDF()->getType() & BSDF::ESmooth) || curr.its.isEmitter()
                        || ++num_bounces > 10){
                        break;
                    }

                    BSDFSamplingRecord bsdf_sample_record(curr.its, sampler);
                    bsdf_sample_record.typeMask = curr.its.getBSDF()->isDielectric() ? 
                        BSDF::EDeltaTransmission | BSDF::ENull : BSDF::EDeltaReflection;
                    curr.its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());

                    curr.ray = Ray(curr.its.p, bsdf_sample_record.its.toWorld(bsdf_sample_record.wo), curr.ray.time);
                }
                curr.intersected = intersected_scene;

                receivers[y * size.x * spp + x * spp + j] = curr;
            }
        }
    }

    return receivers;
}

void computeShadows(BlockingQueue<std::pair<std::uint32_t, std::uint32_t>>& input, 
    BlockingQueue<std::pair<std::uint32_t, std::uint32_t>>& output, 
    std::vector<HWBFPix>& receivers, const std::vector<VPL>& vpls, Scene* scene, float min_dist){
    
    std::pair<std::uint32_t, std::uint32_t> work_unit;
    while(input.pop(work_unit)){
        #pragma omp parallel for
        for(std::uint32_t i = work_unit.first; i < work_unit.second; ++i){
            
            receivers[i].visibility.resize(vpls.size(), 0);
            if(receivers[i].intersected){
                for(std::uint32_t j = 0; j < vpls.size(); ++j){
                    receivers[i].visibility[j] = sampleVisibility(scene, receivers[i].its, vpls[j], min_dist) ? 1 : 0;
                }
            }
        }

        output.push(work_unit);
    }
    
    output.close();
}

void shade(BlockingQueue<std::pair<std::uint32_t, std::uint32_t>>& input, 
    BlockingQueue<std::pair<std::uint32_t, std::uint32_t>>& output,
    std::vector<HWBFPix>& receivers, const std::vector<VPL>& vpls, float min_dist, bool vsl, Scene* scene){

    HWShader hw_shader;

    std::pair<std::uint32_t, std::uint32_t> work_unit;

    while(input.pop(work_unit)){
        hw_shader.renderHWBF(receivers, vpls, work_unit.first, work_unit.second, min_dist, vsl, scene);
        output.push(work_unit);
    }

    output.close();
}

//NB: only suitable when no clustering is needed, ie. with the passthrough clusterer
void LightClustererRenderer::renderHW(Scene* scene, std::uint32_t spp, const RenderJob* job){
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();
    ref<Bitmap> output_image = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, film->getSize());
    Vector2i size = output_image->getSize();
    std::uint8_t *image_buffer = output_image->getUInt8Data();
    memset(image_buffer, 0, output_image->getBytesPerPixel() * size.x * size.y);

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

   

    Intersection its;
    std::vector<VPL> vpls = clusterer_->getClusteringForPoint(its);
    std::uint32_t div_factor = std::max(std::uint32_t(1), std::uint32_t(vpls.size()) / 10000u);
    std::uint32_t batch_size = 250000 / div_factor;
    std::cout << "generating receivers..." << std::endl;
    std::vector<HWBFPix> receivers = generateReceivers(scene, spp, sampler, size);
    std::cout << "rendering..." << std::endl;

    BlockingQueue<std::pair<std::uint32_t, std::uint32_t>> toshadowtest;
    BlockingQueue<std::pair<std::uint32_t, std::uint32_t>> toshade(2);
    BlockingQueue<std::pair<std::uint32_t, std::uint32_t>> tofinish;

    std::uint32_t offset = 0;
    std::uint32_t num_units = 0;
    while(offset < receivers.size()){
        num_units++;
        std::uint32_t end = std::min(std::uint32_t(receivers.size()), offset + batch_size);
        toshadowtest.push(std::make_pair(offset, end));
        offset = end;
    }
    toshadowtest.close();

    std::thread visibility_tester(computeShadows, std::ref(toshadowtest), std::ref(toshade), std::ref(receivers), 
        std::ref(vpls), scene, min_dist_);

    std::thread shader(shade, std::ref(toshade), std::ref(tofinish), std::ref(receivers), std::ref(vpls), 
        min_dist_, vsl_, scene);

    std::unordered_map<std::uint32_t, Spectrum> accumulator;

    ProgressReporter progress("Rendering", num_units, job);
    std::uint32_t finished_units = 0;

    std::pair<std::uint32_t, std::uint32_t> work_unit;
    while(tofinish.pop(work_unit)){
        for(std::uint32_t i = work_unit.first; i < work_unit.second; ++i){
            std::uint32_t buffer_pos = receivers[i].x + receivers[i].y * size.x;
            if(accumulator.find(buffer_pos) != accumulator.end()){
                accumulator[buffer_pos] += receivers[i].col;
            }
            else{
                accumulator[buffer_pos] = receivers[i].col;
            }
        }
        progress.update(finished_units++);
    }

    progress.finish();

    std::cout << std::endl << "Completing render..." << std::endl;

    visibility_tester.join();
    shader.join();

    for(auto iter = accumulator.begin(); iter != accumulator.end(); ++iter){
        Spectrum col = iter->second;
        col /= spp;

        float r, g, b;
        col.toSRGB(r, g, b);

        r = std::max(0.f, std::min(1.f, r));
        g = std::max(0.f, std::min(1.f, g));
        b = std::max(0.f, std::min(1.f, b));

        std::uint32_t buffer_pos = iter->first;

        image_buffer[buffer_pos*3] = r * 255 + 0.5f;
        image_buffer[buffer_pos*3 + 1] = g * 255 + 0.5f;
        image_buffer[buffer_pos*3 + 2] = b * 255 + 0.5f;
    }

    film->setBitmap(output_image);
}

void LightClustererRenderer::renderNHW(Scene* scene, std::uint32_t spp, const RenderJob *job){
    {
        std::lock_guard<std::mutex> lock(cancel_lock_);
        cancel_ = false;
    }
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();
    ref<Bitmap> output_image = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, film->getSize());
    std::uint8_t *image_buffer = output_image->getUInt8Data();
    Vector2i size = output_image->getSize();
    memset(image_buffer, 0, output_image->getBytesPerPixel() * size.x * size.y);

    Properties props("independent");
    Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    std::uint64_t total_lights = 0;
    std::uint64_t total_rendered_pixels = 0;
    std::mutex light_counter_mutex;
    for (std::int32_t y = 0; y < output_image->getSize().y; ++y) {
        {
            std::lock_guard<std::mutex> lock(cancel_lock_);
            if (cancel_) {
                break;
            }
        }
        #pragma omp parallel for
        for (std::int32_t x = 0; x < output_image->getSize().x; ++x) {
            Spectrum accumulator(0.f);

            std::uint32_t cell_dim = sqrt(spp) + 0.5f;
            float cell_side_len = 1.f / cell_dim;

            for(std::uint32_t j = 0; j < spp; ++j){
                Ray ray;

                float x_jitter = 0.5f * cell_side_len;
                float y_jitter = 0.5f * cell_side_len;
                float x_off = (j % cell_dim) * cell_side_len;
                float y_off = (j / cell_dim) * cell_side_len;

                Point2 sample_position(x + x_off + x_jitter, y + y_off + y_jitter);
                
                //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
                Point2 aperture_sample(0.5f, 0.5f);
                Float time_sample(0.5f);

                sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

                Intersection its;

                std::vector<VPL> vpls = clusterer_->getClusteringForPoint(its);
                {
                    std::lock_guard<std::mutex> counter_lock(light_counter_mutex);
                    total_lights += vpls.size();
                    total_rendered_pixels++;
                }
                
                bool intersected;
                for (std::uint32_t i = 0; i < vpls.size(); ++i) {
                    std::uint32_t num_samples;
                    accumulator += sample(scene, sampler, its, ray, vpls[i], min_dist_, true, 
                        10, i == 0, intersected, true, vsl_, num_samples);

                    if(!intersected || its.isEmitter()){
                        break;
                    }
                }

                
            }
            accumulator /= spp;

            float r, g, b;
            accumulator.toSRGB(r, g, b);
            
            //can set the buffer directly since we have direct control over the format of the image
            std::uint32_t offset = (x + y * output_image->getSize().x) * output_image->getBytesPerPixel();
            
            image_buffer[offset] = std::min(1.f, r) * 255 + 0.5f;
            image_buffer[offset + 1] = std::min(1.f, g) * 255 + 0.5f;
            image_buffer[offset + 2] = std::min(1.f, b) * 255 + 0.5f;
        }
    }

    std::cout << (float)total_lights / (float)total_rendered_pixels << std::endl;

    film->setBitmap(output_image);
}

bool LightClustererRenderer::render(Scene* scene, std::uint32_t spp, const RenderJob *job){
    auto start = std::chrono::high_resolution_clock::now();
    if(hw_){
        renderHW(scene, spp, job);
    }
    else{
        renderNHW(scene, spp, job);
    }
    auto end = std::chrono::high_resolution_clock::now();
    
    std::vector<float> timing;
    timing.push_back(std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count());
    writeOutputData(scene, "bftimings-", false, timing, ',');

    return true;
}

void LightClustererRenderer::setCancel(bool cancel){
    std::lock_guard<std::mutex> lock(cancel_lock_);

    cancel_ = cancel;
}

MTS_NAMESPACE_END