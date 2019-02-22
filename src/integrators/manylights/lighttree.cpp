#include "lighttree.h"

#include <mitsuba/core/plugin.h>

#include <tuple>
#include <queue>
#include <functional>
#include <array>
#include <map>

MTS_NAMESPACE_BEGIN

std::unique_ptr<LightTreeNode> createLightTree(const std::vector<VPL>& vpls, EVPLType vpl_type, float min_dist) {
	if(vpls.size() == 0 || vpls[0].type != vpl_type){
		return nullptr;
	}

	std::array<std::vector<std::unique_ptr<LightTreeNode>>, 2> nodes;
	std::uint8_t current_level = 0;
	
	for(size_t i = 0; i < vpls.size(); ++i){
		nodes[current_level].emplace_back(new LightTreeNode());
		
		if(vpl_type != EDirectionalEmitterVPL){
			nodes[current_level][i]->min_bounds = vpls[i].its.p;
			nodes[current_level][i]->max_bounds = vpls[i].its.p;
		}
		
		if(vpl_type != EPointEmitterVPL){
			nodes[current_level][i]->cone_ray1 = vpls[i].its.shFrame.n;
			nodes[current_level][i]->cone_ray2 = vpls[i].its.shFrame.n;
		}
		
		nodes[current_level][i]->vpl = vpls[i];
		nodes[current_level][i]->emission_scale = 1.f;
	}

	typedef std::tuple<float, std::uint32_t, std::uint32_t> SimEntry;

	auto comparator = [](SimEntry l, SimEntry r){
		return std::get<0>(l) > std::get<0>(r);
	};

	Properties props("independent");
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

	while(true){
		if(nodes[current_level].size() == 1){
			break;
		}

		std::priority_queue<SimEntry, std::vector<SimEntry>, decltype(comparator)> similarity_matrix(comparator);

		for(size_t i = 0; i < nodes[current_level].size(); ++i){
			for(size_t j = i + 1; j < nodes[current_level].size(); ++j){
				float d = 0.f;

				if(vpl_type != EDirectionalEmitterVPL){
					Point min(std::min(nodes[current_level][i]->min_bounds[0], nodes[current_level][j]->min_bounds[0]), 
						std::min(nodes[current_level][i]->min_bounds[1], nodes[current_level][j]->min_bounds[1]), 
						std::min(nodes[current_level][i]->min_bounds[2], nodes[current_level][j]->min_bounds[2]));
					Point max(std::max(nodes[current_level][i]->max_bounds[0], nodes[current_level][j]->max_bounds[0]), 
						std::max(nodes[current_level][i]->max_bounds[1], nodes[current_level][j]->max_bounds[1]), 
						std::max(nodes[current_level][i]->max_bounds[2], nodes[current_level][j]->max_bounds[2]));
					
					d += (max - min).length();
				}
				
				if(vpl_type != EPointEmitterVPL){
					Vector3 ray1, ray2;

					ray1 = nodes[current_level][i]->cone_ray1;
					ray2 = nodes[current_level][i]->cone_ray2;
					
					float dp1 = dot(ray1, ray2);
					float dp2 = dot(ray1, nodes[current_level][j]->cone_ray1);
					float dp3 = dot(ray1, nodes[current_level][j]->cone_ray2);

					if(dp2 > dp1 || dp3 > dp1){
						ray2 = dp2 > dp3 ? nodes[current_level][j]->cone_ray1 : nodes[current_level][j]->cone_ray2;
					}

					dp1 = dot(ray1, ray2);
					dp2 = dot(ray2, nodes[current_level][j]->cone_ray1);
					dp3 = dot(ray2, nodes[current_level][j]->cone_ray2);

					if (dp2 > dp1 || dp3 > dp1) {
						ray1 = dp2 > dp3 ? nodes[current_level][j]->cone_ray1 : nodes[current_level][j]->cone_ray2;
					}

					d += (1.f - dot(ray1, ray2)) * min_dist * 2.f;
				}

				similarity_matrix.push(std::make_tuple(d, i, j));
			}
		}

		size_t nodes_left = nodes[current_level].size();
		std::uint8_t next_level = (current_level + 1) % 2;
	
		while(nodes_left > 0){
			if (nodes_left == 1) {
				for (std::uint32_t i = 0; i < nodes[current_level].size(); ++i) {
					if (nodes[current_level][i] != nullptr) {
						nodes[next_level].push_back(std::move(nodes[current_level][i]));
						nodes[current_level][i] = nullptr;
						nodes_left--;
						break;
					}
				}

				continue;
			}

			if(similarity_matrix.empty()){
				std::cerr << "Error: something went terribly wrong" << std::endl;
				exit(0);
 			}

			SimEntry entry = similarity_matrix.top();
			similarity_matrix.pop();

			std::uint32_t id1 = std::get<1>(entry), id2 = std::get<2>(entry);

			if(nodes[current_level][id1] != nullptr && 
				nodes[current_level][id2] != nullptr){
				nodes[next_level].emplace_back(new LightTreeNode());
				LightTreeNode *c1 = nodes[current_level][id1].get(), *c2 = nodes[current_level][id2].get();
				
				if(vpl_type != EDirectionalEmitterVPL){
					nodes[next_level].back()->min_bounds = Point(std::min(c1->min_bounds[0], c2->min_bounds[0]), 
						std::min(c1->min_bounds[1], c2->min_bounds[1]), 
						std::min(c1->min_bounds[2], c2->min_bounds[2]));
					nodes[next_level].back()->max_bounds = Point(std::max(c1->max_bounds[0], c2->max_bounds[0]), 
						std::max(c1->max_bounds[1], c2->max_bounds[1]), 
						std::max(c1->max_bounds[2], c2->max_bounds[2]));
				}
				
				//need to ideally find a better method to union 2 cones
				if(vpl_type != EPointEmitterVPL){
					nodes[next_level].back()->cone_ray1 = c1->cone_ray1;
					nodes[next_level].back()->cone_ray2 = c1->cone_ray2;
					
					float dp1 = dot(nodes[next_level].back()->cone_ray1, nodes[next_level].back()->cone_ray2);
					float dp2 = dot(nodes[next_level].back()->cone_ray1, c2->cone_ray1);
					float dp3 = dot(nodes[next_level].back()->cone_ray1, c2->cone_ray2);

					if(dp2 > dp1 || dp3 > dp1){
						nodes[next_level].back()->cone_ray2 = dp2 > dp3 ? c2->cone_ray1 : c2->cone_ray2;
					}

					dp1 = dot(nodes[next_level].back()->cone_ray1, nodes[next_level].back()->cone_ray2);
					dp2 = dot(nodes[next_level].back()->cone_ray2, c2->cone_ray1);
					dp3 = dot(nodes[next_level].back()->cone_ray2, c2->cone_ray2);

					if (dp2 > dp1 || dp3 > dp1) {
						nodes[next_level].back()->cone_ray1 = dp2 > dp3 ? c2->cone_ray1 : c2->cone_ray2;
					}
				}
				
				float c1_intensity = c1->vpl.P.getLuminance() * c1->emission_scale;
				float c2_intensity = c2->vpl.P.getLuminance() * c2->emission_scale;

				float sample = sampler->next1D();

				if (sample < c1_intensity / (c1_intensity + c2_intensity)) {
					nodes[next_level].back()->vpl = c1->vpl;
					nodes[next_level].back()->emission_scale = c1->emission_scale + c2_intensity / c1->vpl.P.getLuminance();
				}
				else {
					nodes[next_level].back()->vpl = c2->vpl;
					nodes[next_level].back()->emission_scale = c2->emission_scale + c1_intensity / c2->vpl.P.getLuminance();
				}

				nodes[next_level].back()->left = std::move(nodes[current_level][id1]);
				nodes[next_level].back()->right = std::move(nodes[current_level][id2]);

				nodes[current_level][id1] = nullptr;
				nodes[current_level][id2] = nullptr;

				nodes_left -= 2;
			}
		}

		nodes[current_level].clear();
		current_level = next_level;
	}

	return std::move(nodes[current_level][0]);
}

//http://iquilezles.org/www/articles/distfunctions/distfunctions.htm
float boxUDF(Point p, Point box_min, Point box_max){
	assert(p.dim == box_min.dim && p.dim == box_max.dim && p.dim == 3);

	Point box_center = (box_min + box_max) / 2;
	Vector3f box_half_dims = (box_max - box_min) / 2;
	
	Point result;
	result[0] = std::max(0., (double)fabs(p[0] - box_center[0]) - (double)box_half_dims.x);
	result[1] = std::max(0., (double)fabs(p[1] - box_center[1]) - (double)box_half_dims.y);
	result[2] = std::max(0., (double)fabs(p[2] - box_center[2]) - (double)box_half_dims.z);

	return distance(result, Point(0.f));
};

float calculateClusterContribution(Point shading_point_position, Normal shading_point_normal,
	LightTreeNode* light_tree_node, EVPLType vpl_type, float min_dist){

	float geometric = 0.f;

	float d;
	switch(vpl_type){
		case EPointEmitterVPL:
			d = std::max(min_dist, boxUDF(shading_point_position, 
										light_tree_node->min_bounds, light_tree_node->max_bounds));
			geometric = 1.f / (d * d);
			break;
		case ESurfaceVPL:
			{
				d = std::max(min_dist, boxUDF(shading_point_position,
					light_tree_node->min_bounds, light_tree_node->max_bounds));

				//rotation required to transform space so that normal is pointing torwards unit z
				Vector3f unit_z = Vector3f(0.f, 0.f, 1.f);
				float dp = dot(unit_z, shading_point_normal);
				float angle = acos(dp);
				Vector3f axis = abs(1.f - abs(dp)) < 0.0001f ? Vector3f(1.f, 0.f, 0.f) : 
					cross(unit_z, shading_point_normal);

				Transform rotation = Transform::rotate(axis, angle);

				//since we are rotating, the min and max points no longer specify the actual min and max, thus we need
				//the full box
				Point bounding_points[8];
				bounding_points[0] = light_tree_node->min_bounds;
				bounding_points[1] = Point(light_tree_node->min_bounds[0], light_tree_node->min_bounds[1], light_tree_node->max_bounds[2]);
				bounding_points[2] = Point(light_tree_node->min_bounds[0], light_tree_node->max_bounds[1], light_tree_node->min_bounds[2]);
				bounding_points[3] = Point(light_tree_node->max_bounds[0], light_tree_node->min_bounds[1], light_tree_node->min_bounds[2]);
				bounding_points[4] = Point(light_tree_node->min_bounds[0], light_tree_node->max_bounds[1], light_tree_node->max_bounds[2]);
				bounding_points[5] = Point(light_tree_node->max_bounds[0], light_tree_node->min_bounds[1], light_tree_node->max_bounds[2]);
				bounding_points[6] = Point(light_tree_node->max_bounds[0], light_tree_node->max_bounds[1], light_tree_node->min_bounds[2]);
				bounding_points[7] = light_tree_node->max_bounds;

				float max_x2 = std::numeric_limits<float>::min(), max_y2 = std::numeric_limits<float>::min(),
					max_z = std::numeric_limits<float>::min();

				float min_x2 = std::numeric_limits<float>::max(), min_y2 = std::numeric_limits<float>::max();

				for (int i = 0; i < 8; ++i) {
					bounding_points[i] = rotation.transformAffine(bounding_points[i]);

					float x2 = bounding_points[i][0] * bounding_points[i][0];
					float y2 = bounding_points[i][1] * bounding_points[i][1];

					max_x2 = std::max(x2, max_x2);
					min_x2 = std::min(x2, min_x2);
					max_y2 = std::max(y2, max_y2);
					min_y2 = std::min(y2, min_y2);
					max_z = std::max(bounding_points[i][2], max_z);
				}

				float cosine_term = max_z / sqrt(max_z >= 0 ? (min_x2 + min_y2 + max_z * max_z) :
					(max_x2 + max_y2 + max_z * max_z));

				float cosine = dot(normalize(shading_point_position - light_tree_node->vpl.its.p), 
					light_tree_node->vpl.its.shFrame.n);
				geometric = cosine_term / (d * d);
			}

			break;
		case EDirectionalEmitterVPL:
			geometric = 1.f;
			break;
		default:
			break;
	}

	//only dealing with lambertian diffuse right now
	float material = 1.f;

	return geometric * material * light_tree_node->vpl.P.getLuminance() * light_tree_node->emission_scale;
}

LightTree::LightTree() : point_tree_root_(nullptr), directional_tree_root_(nullptr), oriented_tree_root_(nullptr), min_dist_(0.f) {

}

LightTree::LightTree(const std::vector<VPL>& vpls, float min_dist, std::uint32_t max_lights, float error_threshold) : 
	min_dist_(min_dist), max_lights_(max_lights), error_threshold_(error_threshold) {
	for (std::uint32_t i = 0; i < vpls.size(); ++i) {
		switch (vpls[i].type) {
			case EPointEmitterVPL:
				point_vpls_.push_back(vpls[i]);
				break;
			case ESurfaceVPL:
				oriented_vpls_.push_back(vpls[i]);
				break;
			case EDirectionalEmitterVPL:
				directional_vpls_.push_back(vpls[i]);
				break;
			default:
				continue;
		}
	}

	point_tree_root_ = createLightTree(point_vpls_, EPointEmitterVPL, min_dist_);
	oriented_tree_root_ = createLightTree(oriented_vpls_, ESurfaceVPL, min_dist_);
	directional_tree_root_ = createLightTree(directional_vpls_, EDirectionalEmitterVPL, min_dist_);
}

LightTree::LightTree(const LightTree& other) : point_vpls_(other.point_vpls_), directional_vpls_(other.directional_vpls_),
	oriented_vpls_(other.oriented_vpls_), point_tree_root_(new LightTreeNode(*other.point_tree_root_)), 
	directional_tree_root_(new LightTreeNode(*other.directional_tree_root_)),
	oriented_tree_root_(new LightTreeNode(*other.oriented_tree_root_)), min_dist_(other.min_dist_){
}

LightTree::LightTree(LightTree&& other) : point_vpls_(other.point_vpls_), directional_vpls_(other.directional_vpls_),
	oriented_vpls_(other.oriented_vpls_), point_tree_root_(std::move(other.point_tree_root_)), 
	directional_tree_root_(std::move(other.directional_tree_root_)),
	oriented_tree_root_(std::move(other.oriented_tree_root_)), min_dist_(other.min_dist_) {

}

LightTree& LightTree::operator = (const LightTree& other) {
	if(&other != this){
		point_vpls_ = other.point_vpls_;
		directional_vpls_ = other.directional_vpls_;
		oriented_vpls_ = other.oriented_vpls_;
		point_tree_root_ = std::unique_ptr<LightTreeNode>(new LightTreeNode(*other.point_tree_root_));
		directional_tree_root_ = std::unique_ptr<LightTreeNode>(new LightTreeNode(*other.directional_tree_root_));
		oriented_tree_root_ = std::unique_ptr<LightTreeNode>(new LightTreeNode(*other.oriented_tree_root_));
		min_dist_ = other.min_dist_;
	}

	return *this;
}

LightTree& LightTree::operator = (LightTree&& other) {
	if(&other != this){
		point_vpls_ = other.point_vpls_;
		directional_vpls_ = other.directional_vpls_;
		oriented_vpls_ = other.oriented_vpls_;
		point_tree_root_ = std::move(other.point_tree_root_);
		directional_tree_root_ = std::move(other.directional_tree_root_);
		oriented_tree_root_ = std::move(other.oriented_tree_root_);
		min_dist_ = other.min_dist_;
	}

	return *this;
}

LightTree::~LightTree() {
}

std::vector<VPL> LightTree::getClusteringForPoint(const Intersection& its) {
	//std::lock_guard<std::mutex> lock(mutex_);

	std::vector<VPL> lights;

	float point_tree_radiance = point_tree_root_ != nullptr ? 
		calculateClusterContribution(its.p, its.geoFrame.n, point_tree_root_.get(), EPointEmitterVPL, min_dist_) : 0.f;
	float oriented_tree_radiance = oriented_tree_root_ != nullptr ? 
		calculateClusterContribution(its.p, its.geoFrame.n, oriented_tree_root_.get(), EPointEmitterVPL, min_dist_) : 0.f;
	float directional_tree_radiance = directional_tree_root_ != nullptr ? 
		calculateClusterContribution(its.p, its.geoFrame.n, directional_tree_root_.get(), EPointEmitterVPL, min_dist_) : 0.f;

	float total_estimated_radiance = point_tree_radiance + oriented_tree_radiance + directional_tree_radiance;

	typedef std::tuple<LightTreeNode*, float> ClusterAndScore;
	//largest first, but front of the queue is last out
	auto comparator = [](ClusterAndScore l, ClusterAndScore r){
		return std::get<1>(l) < std::get<1>(r);
	};

	std::priority_queue<ClusterAndScore, std::vector<ClusterAndScore>, decltype(comparator)> pqueue(comparator);
	std::map<LightTreeNode*, int> traversed;
	int counter = 0;
	if(point_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(point_tree_root_.get(), point_tree_radiance));
		traversed[point_tree_root_.get()] = counter++;
	}

	if(oriented_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(oriented_tree_root_.get(), oriented_tree_radiance));
		traversed[oriented_tree_root_.get()] = counter++;
	}

	if(directional_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(directional_tree_root_.get(), directional_tree_radiance));
		traversed[directional_tree_root_.get()] = counter++;
	}
	
	while((pqueue.size() + lights.size()) < max_lights_ && pqueue.size() > 0){
		auto entry = pqueue.top();
		
		//if the worst node is below threshold, then all nodes must be
		if(std::get<1>(entry) / total_estimated_radiance < error_threshold_){
			break;
		}

		LightTreeNode* node = std::get<0>(entry);
		
		if(node->left == nullptr && node->right == nullptr){
			//emission scale is 1 for leaf nodes so no need to fix intensity
			lights.push_back(node->vpl);
		}
		else{
			if(node->left == nullptr || node->right == nullptr){
				std::cerr << "A node in the lighttree should always have 2 children" << std::endl;
				exit(0);
			}

			if (traversed.find(node->left.get()) == traversed.end()) {
				float rad = calculateClusterContribution(its.p, its.geoFrame.n, node->left.get(), node->left->vpl.type, min_dist_);
				pqueue.push(std::make_tuple(node->left.get(), rad));
				traversed[node->left.get()] = counter++;
			}

			if (traversed.find(node->right.get()) == traversed.end()) {
				float rad = calculateClusterContribution(its.p, its.geoFrame.n, node->right.get(), node->right->vpl.type, min_dist_);
				pqueue.push(std::make_tuple(node->right.get(), rad));
				traversed[node->right.get()] = counter++;
			}
		}

		pqueue.pop();
	}

	while(pqueue.size() > 0 && lights.size() < max_lights_){
		auto entry = pqueue.top();
		lights.push_back(std::get<0>(entry)->vpl);
		lights.back().P *= std::get<0>(entry)->emission_scale;
		pqueue.pop();
	}

	//should be moved so it should be ok
	return lights;
}

MTS_NAMESPACE_END