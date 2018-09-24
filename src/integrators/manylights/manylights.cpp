#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/core/statistics.h>

#include <mitsuba/core/scopeguard.h>

#include <vector>
#include <algorithm>
#include <mutex>

MTS_NAMESPACE_BEGIN

size_t generateVPLs(const Scene *scene, size_t offset, size_t count, int max_depth, bool prune, 
	std::vector<VPL> &vpls) {
	if (max_depth <= 1)
		return 0;

	Properties props("halton");
	props.setInteger("scramble", 0);
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

	const Sensor *sensor = scene->getSensor();
	Float time = sensor->getShutterOpen() + sensor->getShutterOpenTime() * sampler->next1D();

	const Frame standard_frame(Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1));
	int retries = 0;

	while (vpls.size() < count) {
		sampler->setSampleIndex(++offset);

		if (vpls.empty() && ++retries > 10000) {
			/* Unable to generate VPLs in this scene -- give up. */
			return 0;
		}

		PositionSamplingRecord point_sample(time);
		DirectionSamplingRecord direction_sample;
		Spectrum weight = scene->sampleEmitterPosition(point_sample, sampler->next2D());

		size_t start = vpls.size();

		const Emitter *emitter = static_cast<const Emitter *>(point_sample.object);

		//samples an emitter, adds it to the vpl vector, and then extracts a point and direction sample to be used later for additional
		//vpl generation
		if (!emitter->isEnvironmentEmitter() && emitter->needsDirectionSample()) {
			VPL vpl(EPointEmitterVPL, weight);
			vpl.its.p = point_sample.p;
			vpl.its.time = time;
			vpl.its.shFrame = point_sample.n.isZero() ? standard_frame : Frame(point_sample.n);
			vpl.emitter = emitter;
			vpls.push_back(vpl);

			weight *= emitter->sampleDirection(direction_sample, point_sample, sampler->next2D());
		}
		else {
			/* Hack to get the proper information for directional VPLs */
			DirectSamplingRecord direct_sample(scene->getKDTree()->getAABB().getCenter(), point_sample.time);

			Spectrum direct_sample_weight = emitter->sampleDirect(direct_sample, sampler->next2D())	/ scene->pdfEmitterDiscrete(emitter);

			if (direct_sample_weight.isZero())
				continue;

			VPL vpl(EDirectionalEmitterVPL, direct_sample_weight);
			vpl.its.p = Point(0.0);
			vpl.its.time = time;
			vpl.its.shFrame = Frame(-direct_sample.d);
			vpl.emitter = emitter;
			vpls.push_back(vpl);
			direction_sample.d = -direct_sample.d;

			Point2 offset = warp::squareToUniformDiskConcentric(sampler->next2D());
			Vector transformed_offset = Frame(direct_sample.d).toWorld(Vector(offset.x, offset.y, 0));
			BSphere bounding_sphere = scene->getKDTree()->getAABB().getBSphere();
			point_sample.p = bounding_sphere.center + (transformed_offset - direction_sample.d) * bounding_sphere.radius;
			weight = direct_sample_weight * M_PI * bounding_sphere.radius * bounding_sphere.radius;
		}

		int depth = 2;
		Ray ray(point_sample.p, direction_sample.d, time);
		Intersection its;

		//generates vpls from additional bounces
		while (!weight.isZero() && (depth < max_depth || max_depth == -1)) {
			if (!scene->rayIntersect(ray, its))
				break;

			const BSDF *bsdf = its.getBSDF();
			BSDFSamplingRecord bsdf_sample(its, sampler, EImportance);
			Spectrum bsdf_sample_weight = bsdf->sample(bsdf_sample, sampler->next2D());
			if (bsdf_sample_weight.isZero())
				break;

			float albedo = std::fmin(0.95f, bsdf_sample_weight.max());
			if (sampler->next1D() > albedo) {
				break;
			}
			else {
				weight /= albedo;
			}

			VPL vpl(ESurfaceVPL, weight);
			vpl.its = its;

			if (BSDF::getMeasure(bsdf_sample.sampledType) == ESolidAngle) {
				vpls.push_back(vpl);
			}
				
			weight *= bsdf_sample_weight;

			Vector wi = -ray.d, wo = its.toWorld(bsdf_sample.wo);
			ray = Ray(its.p, wo, 0.0f);

			float wi_dot_geometricnormal = dot(its.geoFrame.n, wi);
			float wo_dot_geometricnormal = dot(its.geoFrame.n, wo);
			if (wi_dot_geometricnormal * Frame::cosTheta(bsdf_sample.wi) <= 0 || wo_dot_geometricnormal * Frame::cosTheta(bsdf_sample.wo) <= 0) {
				break;
			}

			++depth;
		}

		size_t end = vpls.size();
		for (size_t i = start; i < end; ++i)
			vpls[i].emitterScale = 1.0f / (end - start);
	}

	return offset;
}

class ManyLightsIntegrator : public Integrator {
public:
	ManyLightsIntegrator(const Properties &props) : 
		Integrator(props),
		max_depth_(props.getInteger("maxDepth", 5)){
		
	}

	//generation of the point lights is done the same way as the vpl integrator for now, will need to extend this in the future
	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {

		Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);

		if (!(scene->getSensor()->getType() & Sensor::EProjectiveCamera))
			Log(EError, "The Many Lights integrator requires a projective camera "
				"(e.g. perspective/thinlens/orthographic/telecentric)!");

		vpls_.clear();
		size_t sample_count = scene->getSampler()->getSampleCount();
		vpls_.reserve(sample_count);
		float normalization = 1.f / generateVPLs(scene, 0, sample_count, max_depth_, true, vpls_);
		
		for (size_t i = 0; i < vpls_.size(); ++i) {
			vpls_[i].P *= normalization;
			vpls_[i].emitterScale *= normalization;
		}

		Log(EInfo, "Generated %i virtual point lights", vpls_.size());

		return true;
	}

	void cancel() {
		std::lock_guard<std::mutex> lock(cancel_lock_);
		cancel_ = true;
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job, int sceneResID, int sensorResID, int samplerResID) {
		
		{
			std::lock_guard<std::mutex> lock(cancel_lock_);
			cancel_ = false;
		}

		ref<Sensor> sensor = scene->getSensor();
		ref<Film> film = sensor->getFilm();
		if (output_image_ == nullptr || output_image_->getSize() != film->getSize()) {
			output_image_ = createBitmap(film->getSize());
		}

		std::uint8_t *image_buffer = output_image_->getUInt8Data();
		memset(image_buffer, 0, output_image_->getBytesPerPixel() * output_image_->getSize().y * output_image_->getSize().x);

		for (std::int32_t y = 0; y < output_image_->getSize().y; ++y) {
			for (std::int32_t x = 0; x < output_image_->getSize().x; ++x) {
				{
					std::lock_guard<std::mutex> lock(cancel_lock_);
					if (cancel_) {
						break;
					}
				}

				Ray ray;

				Point2 sample_position(x + 0.5f, y + 0.5f);
				
				//disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
				Point2 aperture_sample(0.5f, 0.5f);
				Float time_sample(0.5f);

				sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

				Float t;
				ConstShapePtr shape;
				Normal n;
				Point2 uv;

				size_t offset = (x + output_image_->getSize().x * y) * output_image_->getBytesPerPixel();

				if (scene->rayIntersect(ray, t, shape, n, uv)) {
					Point world_p = ray.o + ray.d * t;

					std::vector<VPL> vpls;
					getVPLs(vpls);

					float r = 0.f, g = 0.f, b = 0.f;
					for (std::uint32_t i = 0; i < vpls.size(); ++i) {
						float d = (world_p - vpls[i].its.p).length();
						d = std::min(0.01f, d);
						float power = 1. / (d * d);

						float lr, lg, lb;
						vpls[i].P.toLinearRGB(lr, lg, lb);

						float n_dot_ldir = std::max(0.f, dot(n, vpls[i].its.p - world_p));

						r += lr * n_dot_ldir * power;
						g += lg * n_dot_ldir * power;
						b += lb * n_dot_ldir * power;
					}

					r = std::min(r, 1.f);	
					g = std::min(g, 1.f);	
					b = std::min(b, 1.f);

					image_buffer[offset] = r * 255. + 0.5;
					image_buffer[offset + 1] = g * 255. + 0.5;
					image_buffer[offset + 2] = b * 255. + 0.5;
				}

				image_buffer[offset] = 0xff;
			}
		}

		film->setBitmap(output_image_);

		return !cancel_;
	}

	MTS_DECLARE_CLASS()

private:
	ref<Bitmap> createBitmap(const Vector2i& dimensions) {
		if(dimensions.x <= 0 || dimensions.y <= 0){
			return nullptr;
		}
		
		//make 8bit channels for now, can change it in the future
		return new Bitmap(Bitmap::ESpectrum, Bitmap::EUInt8, dimensions);
	}

	//brute force for now, implement in selection strategy selection later
	void getVPLs(std::vector<VPL>& vpls) {
		vpls = vpls_;
	}

private:
	ref<Bitmap> output_image_;
	std::vector<VPL> vpls_;
	std::uint32_t max_depth_;
	std::mutex cancel_lock_;
	bool cancel_;
};

MTS_IMPLEMENT_CLASS(ManyLightsIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(ManyLightsIntegrator, "Many-Lights-based integrator");
MTS_NAMESPACE_END
