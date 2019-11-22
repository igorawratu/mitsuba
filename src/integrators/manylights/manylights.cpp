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
#include "common.h"

#include <flann/flann.hpp>

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
		bool surface_emitter = vpls[i].type == ESurfaceVPL;
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

void setVPLRadii(std::vector<VPL>& vpls, float min_dist){
	//radius calculation, 11 to account 10 + 1 for adding a node's self in nearest neighbours
	std::uint32_t num_neighbours = std::min((std::uint32_t)vpls.size(), 11u);

	std::uint32_t num_sl = 0;
	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(vpls[i].type == ESurfaceVPL){
			num_sl++;
		}
	} 

    flann::Matrix<float> dataset(new float[num_sl * 3], num_sl, 3);
	std::uint32_t curr_light = 0;
	std::vector<std::uint32_t> pointlight_indices;
    for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(vpls[i].type == ESurfaceVPL){
			float* curr = (float*)dataset[curr_light++];
			curr[0] = vpls[i].its.p.x;
			curr[1] = vpls[i].its.p.y;
			curr[2] = vpls[i].its.p.z;
			pointlight_indices.emplace_back(i);
		}
		else{
			vpls[i].radius = 0.f;
		}
    }

    flann::Index<flann::L2<float>> index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();

	std::vector<std::vector<int>> indices;
	std::vector<std::vector<float>> distances;

    index.knnSearch(dataset, indices, distances, num_neighbours, flann::SearchParams(128));

	for(std::uint32_t i = 0; i < distances.size(); ++i){
		float max = std::numeric_limits<float>::min();
		float min = std::numeric_limits<float>::max();
		for(std::uint32_t j = 0; j < distances[i].size(); ++j){
			if(distances[i][j] > std::numeric_limits<float>::epsilon()){
				max = std::max(distances[i][j], max);
				min = std::min(distances[i][j], min);
			}
		}

		vpls[pointlight_indices[i]].radius = sqrt(max) * 2.f;
	}
}

size_t generateVPLs(const Scene *scene, size_t offset, size_t count, int max_depth, bool prune, 
	std::vector<VPL> &vpls) {
	if (max_depth <= 1)
		return 0;

	Properties props("halton");
	props.setInteger("scramble", 0);
	ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
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
		if (!emitter->isEnvironmentEmitter()) {
			EVPLType type = emitter->needsDirectionSample() ? ESurfaceVPL : EPointEmitterVPL;
			VPL vpl(type, weight);
			vpl.its.p = point_sample.p;
			vpl.its.time = time;
			vpl.its.shFrame = point_sample.n.isZero() ? standard_frame : Frame(point_sample.n);
			vpl.emitter = emitter;
			vpl.psr = point_sample;
			vpls.push_back(vpl);

			if(type == ESurfaceVPL){
				weight *= emitter->sampleDirection(direction_sample, point_sample, sampler->next2D());
			}
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
		while (!weight.isZero() && (depth < max_depth || max_depth == -1) && vpls.size() < count) {
			if (!scene->rayIntersect(ray, its))
				break;

			const BSDF *bsdf = its.getBSDF();
			BSDFSamplingRecord bsdf_sample(its, sampler, EImportance);
			Spectrum bsdf_sample_weight = bsdf->sample(bsdf_sample, sampler->next2D());
			if (bsdf_sample_weight.isZero())
				break;

			/*float approx_albedo = std::fmin(0.95f, bsdf_sample_weight.max());

            if (sampler->next1D() > approx_albedo){
				break;
			}
            else{
				weight /= approx_albedo;
			}*/
			
			VPL vpl(ESurfaceVPL, weight);
			vpl.its = its;

			weight *= bsdf_sample_weight;

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
		max_depth_(props.getInteger("max_depth", 5)),
		properties_(props),
		clamping_(props.getFloat("clamping", 0.1f)),
		clustering_strategy_(NONE), 
		renderer_(nullptr){
		int strategy = props.getInteger("clustering_strat", 0);
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
		std::uint32_t upper_sqrt = ceil(sqrt(scene->getSampler()->getSampleCount())) + 0.5f;
		spp_ = upper_sqrt * upper_sqrt;//upperPo2(scene->getSampler()->getSampleCount());
		size_t sample_count = properties_.getInteger("vpls", 1024);
		vpls_.reserve(sample_count);
		float normalization = 1.f / generateVPLs(scene, 0, sample_count, max_depth_, true, vpls_);
		
		for (size_t i = 0; i < vpls_.size(); ++i) {
			vpls_[i].P *= normalization;
			vpls_[i].emitterScale *= normalization;
		}

		min_dist_ = calculateMinDistance(scene, vpls_, clamping_);

		setVPLRadii(vpls_, min_dist_);

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
			return renderer_->render(scene, spp_, job);
		}

		return false;
	}

	MTS_DECLARE_CLASS()

private:
	std::unique_ptr<ManyLightsRenderer> initializeRendererByStrategy(const Scene* scene, const Properties& props, CLUSTERING_STRATEGY strategy){
		bool vsl = props.getInteger("vsl", 0) > 0;
		switch(strategy){
			case NONE:
			{
				std::unique_ptr<ManyLightsClusterer> clusterer(new PassthroughClusterer(vpls_));
				return std::unique_ptr<ManyLightsRenderer>(new LightClustererRenderer(std::move(clusterer), min_dist_, vsl));
			}
			case LIGHTCUTS:
			{
				int max_lights = props.getInteger("lightcuts-num_clusters", 64);
				float error_threshold = props.getFloat("lightcuts-error_threshold", 0.02);
				std::unique_ptr<ManyLightsClusterer> clusterer(new LightTree(vpls_, min_dist_, max_lights, error_threshold));
				return std::unique_ptr<ManyLightsRenderer>(new LightClustererRenderer(std::move(clusterer), min_dist_, vsl));
			}
			case ROWCOLSAMPLING:
			{
				int rows = props.getInteger("rowcol-rows", 300);
				int cols = props.getInteger("rowcol-num_clusters", 64);

				std::unique_ptr<ManyLightsClusterer> clusterer(new RowColumnSampling(vpls_, rows, cols, std::make_tuple(scene->getFilm()->getSize().x, scene->getFilm()->getSize().y),
					scene, min_dist_, vsl));

				return std::unique_ptr<ManyLightsRenderer>(new LightClustererRenderer(std::move(clusterer), min_dist_, vsl));
			}
			case MATRIXRECONSTRUCTION:
			{
				float sample_percentage = props.getFloat("completion-sample_perc", 0.1);
				float step_size_factor = props.getFloat("completion-step_factor", 1.5);
				float tolerance = props.getFloat("completion-tolerance", 0.01);
				float tau = props.getFloat("completion-tau", 5);
				int max_iterations = props.getInteger("completion-reconstruction_iterations", 20);
				std::uint32_t slice_size = props.getInteger("completion-slice_size", 1024);
				bool visibility_only = props.getInteger("completion-visibility_only", 0) > 0;
				bool adaptive_col_sampling = props.getInteger("completion-adaptive_col_sampling", 0) > 0;
				bool adaptive_importance_sampling = props.getInteger("completion-adaptive_importance_sampling", 0) > 0;
				bool adaptive_force_resample = props.getInteger("completion-adaptive_force_resample", 0) > 0;
				bool adaptive_recover_transpose = props.getInteger("completion-adaptive_recover_transpose", 0) > 0;
				bool truncated = props.getInteger("completion-use_truncated", 0) > 0;
				bool show_slices = props.getInteger("completion-show_slices", 0) > 0;
				bool show_stats = props.getInteger("completion-show_stats", 0) > 0;
				bool show_svd = props.getInteger("completion-show_svd", 0) > 0;
				float error_scale = props.getFloat("completion-error_scale", 1.f);
				ClusteringStrategy cs = props.getString("completion-cluster_strat", "ls") == "ls" ? 
					ClusteringStrategy::LS : ClusteringStrategy::MDLC;
				bool hw = props.getInteger("completion-hw", 0) > 0;
				bool bin_vis = props.getInteger("completion-bvis", 0) > 0;
				int num_clusters = props.getInteger("completion-clusters_per_slice", 1000);

				//std::unique_ptr<ManyLightsClusterer> clusterer(new PassthroughClusterer(vpls_));
				return std::unique_ptr<ManyLightsRenderer>(new MatrixReconstructionRenderer(vpls_, sample_percentage, min_dist_, 
					step_size_factor, tolerance, tau, max_iterations, slice_size, visibility_only, adaptive_col_sampling, 
					adaptive_importance_sampling, adaptive_force_resample, adaptive_recover_transpose, truncated, show_slices, vsl,
					show_stats, show_svd, cs, error_scale, hw, bin_vis, num_clusters));
			}
			case MATRIXSEPARATION:
			{
				float sample_percentage = props.getFloat("separation-sample_perc", 0.05);
				float error_threshold = props.getFloat("separation-prediction_error_thresh", 0.1);
				float reincorporation_density_threshold = props.getFloat("separation-density_thresh", 0.2);
        		std::uint32_t slice_size = props.getInteger("separation-slice_size", 1024);
				std::uint32_t max_prediction_iterations = props.getInteger("separation-prediction_iter", 10);
				std::uint32_t max_separation_iterations = props.getInteger("separation-separation_iter", 20);
				std::uint32_t show_slices = props.getInteger("separation-show_slices", 0);
				std::uint32_t only_directsamples = props.getInteger("separation-only_direct_samples", 0);
				bool separate = props.getInteger("separation-separate", 1) > 0;
				bool show_error = props.getInteger("separation-show_error", 0) > 0;
				bool show_sparse = props.getInteger("separation-show_sparse", 0) > 0;
				std::uint32_t predictor_mask = (std::uint32_t)props.getInteger("separation-predictor_mask", 7);
				bool show_rank = props.getInteger("separation-show_rank", 0) > 0;
				bool show_predictors = props.getInteger("separation-show_predictors", 0) > 0;
				float rank_increase_threshold = props.getFloat("separation-rank_increase_threshold", 0.1);
				float theta = props.getFloat("separation-theta", 0.01);

				std::unique_ptr<ManyLightsClusterer> clusterer(new PassthroughClusterer(vpls_));
				return std::unique_ptr<MatrixSeparationRenderer>(new MatrixSeparationRenderer(std::move(clusterer), 
					min_dist_, sample_percentage, error_threshold, reincorporation_density_threshold, slice_size, 
					max_prediction_iterations, max_separation_iterations, show_slices, only_directsamples, separate,
					show_error, show_sparse, predictor_mask, show_rank, show_predictors, rank_increase_threshold,
					theta, vsl));
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
	std::uint32_t spp_;
};

MTS_IMPLEMENT_CLASS(ManyLightsIntegrator, false, Integrator)
MTS_EXPORT_PLUGIN(ManyLightsIntegrator, "Many-Lights-based integrator");
MTS_NAMESPACE_END
