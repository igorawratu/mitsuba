#include <mitsuba/core/plugin.h>
#include <mitsuba/hw/vpl.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/qmc.h>
#include <mitsuba/core/scopeguard.h>

#include <vector>
#include <algorithm>
#include <mutex>
#include <iostream>

#include "lighttree.h"
#include "rowcolumnsampling.h"
#include "matrixreconstruction.h"
#include "lightclustererrenderer.h"
#include "manylightsbase.h"
#include "passthroughclusterer.h"
#include "matrixseparationrenderer.h"

#include "definitions.h"

MTS_NAMESPACE_BEGIN

enum CLUSTERING_STRATEGY{
		NONE = 0, LIGHTCUTS, ROWCOLSAMPLING, MATRIXRECONSTRUCTION, MATRIXSEPARATION,
		CLUSTER_STRATEGY_MIN = NONE, CLUSTER_STRATEGY_MAX = MATRIXSEPARATION 
	};

float calculateMinDistance(const Scene *scene, const std::vector<VPL>& vpls, float clamping){
	float acc_near = 0.f;
	float acc_far = 0.f;
	std::uint32_t lights_processed = 0;

	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		bool surface_emitter = vpls[i].emitter->getType() & Emitter::EOnSurface;
		float near = std::numeric_limits<float>::max();
		float far = std::numeric_limits<float>::min();

		for(int k = 0; k < 128; ++k){
			if(vpls[i].type == EDirectionalEmitterVPL){
				continue;
			}

			Point2 sample(sample02(k));

			Ray ray(vpls[i].its.p,
				surface_emitter ? vpls[i].its.shFrame.toWorld(warp::squareToCosineHemisphere(sample))
				: warp::squareToUniformSphere(sample), 0);

			Float t;
			ConstShapePtr shape;
			Normal n;
			Point2 uv;

			float acc = 0.f;
			bool hit = false;

			for(std::uint8_t j = 0; j < 5; ++j){
				if(!scene->rayIntersect(ray, t, shape, n, uv)){
					break;
				}

				acc += t;
				if (!(shape->getBSDF()->getType() & BSDF::ETransmission)) {
					hit = true;
					break;
				}

				ray.o = ray(t);
			}
			
			if (hit && acc > 0.f) {
				near = std::min(near, acc);
				far = std::max(far, acc);
			}

		}

		if(near >= far){
			BSphere bsphere(scene->getKDTree()->getAABB().getBSphere());
			Float min_dist = 0;

			if ((vpls[i].type == ESurfaceVPL || vpls[i].type == EPointEmitterVPL) &&
					!bsphere.contains(vpls[i].its.p))
				min_dist = (bsphere.center - vpls[i].its.p).length() - bsphere.radius;

			far = min_dist + bsphere.radius * 2.25f;
			near = std::max(min_dist - 0.25f * bsphere.radius, far * 1e-5f);
		}
		else{
			near = std::max(near / 1.5f, far * 1e-5f);
			far *= 1.5f;
		}

		acc_near += near;
		acc_far += far;
		lights_processed++;
	}

	acc_near /= (float)lights_processed;
	acc_far /= (float)lights_processed;

	return acc_near + (acc_far - acc_near) * clamping;
}

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
			VPL vpl(ESurfaceVPL, weight);
			vpl.its.p = point_sample.p;
			vpl.its.time = time;
			vpl.its.shFrame = point_sample.n.isZero() ? standard_frame : Frame(point_sample.n);
			vpl.emitter = emitter;
			vpls.push_back(vpl);

			weight *= emitter->sampleDirection(direction_sample, point_sample, sampler->next2D());
		}
		else {
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

			float approx_albedo = std::fmin(0.95f, bsdf_sample_weight.max());

            if (sampler->next1D() > approx_albedo){
				break;
			}
            else{
				weight /= approx_albedo;
			}
			
			weight *= bsdf_sample_weight;
			
			VPL vpl(ESurfaceVPL, weight);
			vpl.its = its;
			
			if (BSDF::getMeasure(bsdf_sample.sampledType) == ESolidAngle) {
				vpls.push_back(vpl);
			}

			Vector wi = -ray.d, wo = its.toWorld(bsdf_sample.wo);

			ray = Ray(its.p, wo, 0.0f);

			float wi_dot_geometricnormal = dot(its.geoFrame.n, wi);
			float wo_dot_geometricnormal = dot(its.geoFrame.n, wo);
			if (wi_dot_geometricnormal * Frame::cosTheta(bsdf_sample.wi) <= 0 || 
				wo_dot_geometricnormal * Frame::cosTheta(bsdf_sample.wo) <= 0) {
				break;
			}

			++depth;
		}

		size_t end = vpls.size();
		float lights_added = end - start;
		for (size_t i = start; i < end; ++i){
			vpls[i].emitterScale = 1.0f / lights_added;
		}
	}

	return offset;
}

class ManyLightsIntegrator : public Integrator {
public:
	ManyLightsIntegrator(const Properties &props) : 
		Integrator(props),
		max_depth_(props.getInteger("maxDepth", 5)),
		properties_(props),
		clamping_(props.getFloat("clamping", 0.1f)),
		clustering_strategy_(NONE), 
		renderer_(nullptr){
		int strategy = props.getInteger("clusteringStrategy", 0);
		clustering_strategy_ = strategy < CLUSTER_STRATEGY_MIN || strategy > CLUSTER_STRATEGY_MAX ? 
			NONE : (CLUSTERING_STRATEGY)strategy;
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

		min_dist_ = calculateMinDistance(scene, vpls_, clamping_);

		renderer_ = initializeRendererByStrategy(scene, properties_, clustering_strategy_);

		Log(EInfo, "Generated %i virtual point lights", vpls_.size());

		return true;
	}

	void cancel() {
		if(renderer_ != nullptr){
			renderer_->setCancel(true);
		}
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job, int sceneResID, int sensorResID, int samplerResID) {
		if(renderer_ != nullptr){
			return renderer_->render(scene);
		}

		return false;
	}

	MTS_DECLARE_CLASS()

private:
	std::unique_ptr<ManyLightsRenderer> initializeRendererByStrategy(const Scene* scene, const Properties& props, CLUSTERING_STRATEGY strategy){
		switch(strategy){
			case NONE:
			{
				std::unique_ptr<ManyLightsClusterer> clusterer(new PassthroughClusterer(vpls_));
				return std::unique_ptr<ManyLightsRenderer>(new LightClustererRenderer(std::move(clusterer), min_dist_));
			}
			case LIGHTCUTS:
			{
				int max_lights = props.getInteger("numClusters", 64);
				float error_threshold = props.getFloat("error_threshold", 0.02);
				std::unique_ptr<ManyLightsClusterer> clusterer(new LightTree(vpls_, min_dist_, max_lights, error_threshold));
				return std::unique_ptr<ManyLightsRenderer>(new LightClustererRenderer(std::move(clusterer), min_dist_));
			}
			case ROWCOLSAMPLING:
			{
				int rows = props.getInteger("rows", 300);
				int cols = props.getInteger("numClusters", 64);

				std::unique_ptr<ManyLightsClusterer> clusterer(new RowColumnSampling(vpls_, rows, cols, std::make_tuple(scene->getFilm()->getSize().x, scene->getFilm()->getSize().y),
					scene, min_dist_));

				return std::unique_ptr<ManyLightsRenderer>(new LightClustererRenderer(std::move(clusterer), min_dist_));
			}
			case MATRIXRECONSTRUCTION:
			{
				std::uint32_t bucket_width = props.getInteger("bucketWidth", 2);
				std::uint32_t bucket_height = props.getInteger("bucketHeight", 1);
				float samples_per_bucket = props.getFloat("samplesPerBucket", 0.5);
				int light_samples = samples_per_bucket * vpls_.size();
				float step_size_factor = props.getFloat("stepSizeFactor", 1.5);
				float tolerance = props.getFloat("tolerance", 0.0001);
				float tau = props.getFloat("tau", 5);
				int max_iterations = props.getInteger("maxIterations", 1000);

				std::unique_ptr<ManyLightsClusterer> clusterer(new PassthroughClusterer(vpls_));
				return std::unique_ptr<ManyLightsRenderer>(new MatrixReconstructionRenderer(std::move(clusterer), 
					std::make_pair(bucket_width, bucket_height), light_samples, min_dist_, step_size_factor, 
					tolerance, tau, max_iterations));
			}
			case MATRIXSEPARATION:
			{
				float sample_percentage = props.getFloat("sample_perc", 0.05);
				float error_threshold = props.getFloat("prediction_error_thresh", 0.1);
				float reincorporation_density_threshold = props.getFloat("density_thresh", 0.2);
        		std::uint32_t slice_size = props.getInteger("slice_size", 1024);
				std::uint32_t max_prediction_iterations = props.getInteger("prediction_iter", 10);
				std::uint32_t max_separation_iterations = props.getInteger("separation_iter", 20);
				std::uint32_t show_slices = props.getInteger("show_slices", 0);
				std::uint32_t only_directsamples = props.getInteger("only_direct_samples", 0);
				bool separate = props.getInteger("separate", 1) > 0;
				bool show_error = props.getInteger("show_error", 0) > 0;

				std::unique_ptr<ManyLightsClusterer> clusterer(new PassthroughClusterer(vpls_));
				return std::unique_ptr<MatrixSeparationRenderer>(new MatrixSeparationRenderer(std::move(clusterer), 
					min_dist_, sample_percentage, error_threshold, reincorporation_density_threshold, slice_size, 
					max_prediction_iterations, max_separation_iterations, show_slices, only_directsamples, separate,
					show_error));
			}
			default:
				return nullptr;
		}
	}

private:
	ref<Bitmap> output_image_;
	std::vector<VPL> vpls_;
	std::uint32_t max_depth_;

	Properties properties_;
	float clamping_;
	CLUSTERING_STRATEGY clustering_strategy_;
	std::unique_ptr<ManyLightsRenderer> renderer_;
	float min_dist_;
};

MTS_IMPLEMENT_CLASS(ManyLightsIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(ManyLightsIntegrator, "Many-Lights-based integrator");
MTS_NAMESPACE_END
