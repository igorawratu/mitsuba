#include "lightclustererrenderer.h"

#include <mitsuba/core/plugin.h>
#include "definitions.h"
#include "common.h"

MTS_NAMESPACE_BEGIN

LightClustererRenderer::LightClustererRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, float min_dist, bool vsl) : 
    min_dist_(min_dist), clusterer_(std::move(clusterer)), vsl_(vsl), cancel_(false){
}

LightClustererRenderer::LightClustererRenderer(LightClustererRenderer&& other) : clusterer_(std::move(other.clusterer_)),
    vsl_(other.vsl_), cancel_(other.cancel_){
}

LightClustererRenderer& LightClustererRenderer::operator = (LightClustererRenderer&& other){
    if(this != &other){
        this->clusterer_ = std::move(other.clusterer_);
        this->vsl_ = other.vsl_;
        this->cancel_ = other.cancel_;
    }

    return *this;
}

LightClustererRenderer::~LightClustererRenderer(){
}

bool LightClustererRenderer::render(Scene* scene){
    {
        std::lock_guard<std::mutex> lock(cancel_lock_);
        cancel_ = false;
    }

    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();
    ref<Bitmap> output_image = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, film->getSize());
    std::uint8_t *image_buffer = output_image->getUInt8Data();

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
            Ray ray;

            Point2 sample_position(x + 0.5f, y + 0.5f);
            
            //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
            Point2 aperture_sample(0.5f, 0.5f);
            Float time_sample(0.5f);

            sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

            Intersection its;

            Spectrum accumulator(0.f);

            std::vector<VPL> vpls = clusterer_->getClusteringForPoint(its);
            {
                std::lock_guard<std::mutex> counter_lock(light_counter_mutex);
                total_lights += vpls.size();
                total_rendered_pixels++;
            }
            
            bool intersected;
            for (std::uint32_t i = 0; i < vpls.size(); ++i) {
                accumulator += sample(scene, sampler, its, ray, vpls[i], min_dist_, true, 
                    10, i == 0, intersected, true, vsl_);

                if(!intersected || its.isEmitter()){
                    break;
                }
            }

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

    return !cancel_;
}

void LightClustererRenderer::setCancel(bool cancel){
    std::lock_guard<std::mutex> lock(cancel_lock_);

    cancel_ = cancel;
}

MTS_NAMESPACE_END