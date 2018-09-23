#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/core/statistics.h>

#include <mitsuba/core/scopeguard.h>

#include <vector>
#include <algorithm>


MTS_NAMESPACE_BEGIN

size_t generateVPLs(const Scene *scene, size_t offset, size_t count, int max_depth, bool prune, 
	std::vector<VPL> &vpls) {
	if (max_depth <= 1)
		return 0;

	Properties props("halton");
	props.setInteger("scramble", 0);
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	auto sampler_scopeguard = makeScopeGuard([&sampler] {sampler->decRef(true); });
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

	}

	bool render(Scene *scene, RenderQueue *queue,
		const RenderJob *job, int sceneResID, int sensorResID, int samplerResID) {
		
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<VPL> vpls_;
	std::uint32_t max_depth_;
};

MTS_IMPLEMENT_CLASS(ManyLightsIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(ManyLightsIntegrator, "Many-Lights-based integrator");
MTS_NAMESPACE_END
