#include "lighttree.h"

#include <mitsuba/core/plugin.h>

#include <tuple>
#include <queue>
#include <functional>
#include <array>
#include <map>
#include <random>
#include <chrono>
#include <math.h>

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
			nodes[current_level][i]->cone_ray = vpls[i].its.shFrame.n;
			nodes[current_level][i]->cone_halfangle = 0.f;
		}
		
		nodes[current_level][i]->vpl = vpls[i];
		float r, g, b;
		vpls[i].P.toLinearRGB(r, g, b);
		nodes[current_level][i]->emission_scale = Vector3f(r, g, b).length();
		nodes[current_level][i]->vpl.P /= nodes[current_level][i]->emission_scale;
	}

	typedef std::tuple<float, std::uint32_t, std::uint32_t, Vector3f, float> SimEntry;

	auto comparator = [](SimEntry l, SimEntry r){
		return std::get<0>(l) < std::get<0>(r);
	};

	std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<float> gen(0, 1);

	while(true){
		if(nodes[current_level].size() == 1){
			break;
		}

		std::priority_queue<SimEntry, std::vector<SimEntry>, decltype(comparator)> similarity_matrix(comparator);

		for(size_t i = 0; i < nodes[current_level].size(); ++i){
			for(size_t j = i + 1; j < nodes[current_level].size(); ++j){
				float d = 0.f;
				float union_angle_span = 0.f;
				Vector3f bcone(0.f);

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
					Vector3 ray1 = nodes[current_level][i]->cone_ray;
					Vector3 ray2 = nodes[current_level][j]->cone_ray;
					float ha1 = nodes[current_level][i]->cone_halfangle;
					float ha2 = nodes[current_level][j]->cone_halfangle;

					float angle = acos(dot(ray1, ray2));
					float max_child_half = std::max(ha1, ha2);
					float min_child_half = std::min(ha1, ha2);
					
					if((angle < max_child_half && (max_child_half - angle) > min_child_half) || (ray1 - ray2).length() < std::numeric_limits<float>::epsilon()){
						union_angle_span = max_child_half;
						bcone = ha1 > ha2 ? ray1 : ray2;
					}
					else{
						union_angle_span = angle + ha1 + ha2;
						union_angle_span = union_angle_span > M_PI ? 2 * M_PI : union_angle_span;
						union_angle_span = std::min(union_angle_span / 2.f, M_PI);

						Vector3f axis = cross(ray1, ray2);
						if(axis.length() < std::numeric_limits<float>::epsilon()){
							axis = cross(ray1, Vector3f(gen(rng), gen(rng), gen(rng)));
						}

						Point new_coneray = Transform::rotate(axis, -ha1).transformAffine(Point(ray1));
						bcone = Vector(Transform::rotate(axis, union_angle_span).transformAffine(new_coneray));
					}

					

					d += union_angle_span * union_angle_span * min_dist * 2.f;
				}

				similarity_matrix.push(std::make_tuple(d, i, j, bcone, union_angle_span));
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
				
				float c1_intensity = c1->emission_scale;
				float c2_intensity = c2->emission_scale;
				float ratio = c1_intensity / (c1_intensity + c2_intensity);

				float sample = gen(rng);

				if(c1_intensity > 0.f || c2_intensity > 0.f){
					if (sample < c1_intensity / (c1_intensity + c2_intensity)) {
						nodes[next_level].back()->vpl = c1->vpl;
					}
					else {
						nodes[next_level].back()->vpl = c2->vpl;
					}

					nodes[next_level].back()->cone_ray = std::get<3>(entry);
					nodes[next_level].back()->cone_halfangle = std::get<4>(entry);
					nodes[next_level].back()->vpl.its.p = Point(ratio * Vector3f(c1->vpl.its.p) + (1.f - ratio) * Vector3f(c2->vpl.its.p));
					Vector3f new_n = ratio * Vector3f(c1->vpl.its.shFrame.n) + (1.f - ratio) * Vector3f(c2->vpl.its.shFrame.n);
					if(new_n.length() < 1e-10f){
						new_n = cross(Vector3f(c1->vpl.its.shFrame.n), Vector3f(gen(rng), gen(rng), gen(rng)));
					}
					nodes[next_level].back()->vpl.its.shFrame = Frame(normalize(new_n));
					nodes[next_level].back()->vpl.P = ratio * c1->vpl.P + (1.f - ratio) * c2->vpl.P;
				}

				nodes[next_level].back()->emission_scale = c1_intensity + c2_intensity;

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

float distToBox(Point p, Point box_min, Point box_max){
	assert(p.dim == box_min.dim && p.dim == box_max.dim && p.dim == 3);

	float dist = 0.f;

	for(std::uint8_t i = 0; i < 3; ++i){
		float curr_dim_d = 0.f;
		if(p[i] < box_min[i]){
			curr_dim_d = box_min[i] - p[i];
		}
		else if(p[i] > box_max[i]){
			curr_dim_d = p[i] - box_max[i];
		}
		dist += curr_dim_d * curr_dim_d;
	}

	return sqrt(dist);
};

float calculateClusterContribution(Point shading_point_position, Normal shading_point_normal,
	LightTreeNode* light_tree_node, EVPLType vpl_type, float min_dist){

	float geometric = 0.f;

	Point bounding_points[8];
	bounding_points[0] = light_tree_node->min_bounds;
	bounding_points[1] = Point(light_tree_node->min_bounds[0], light_tree_node->max_bounds[1], light_tree_node->min_bounds[2]);
	bounding_points[2] = Point(light_tree_node->max_bounds[0], light_tree_node->max_bounds[1], light_tree_node->min_bounds[2]);
	bounding_points[3] = Point(light_tree_node->max_bounds[0], light_tree_node->min_bounds[1], light_tree_node->min_bounds[2]);
	bounding_points[4] = Point(light_tree_node->min_bounds[0], light_tree_node->min_bounds[1], light_tree_node->max_bounds[2]);
	bounding_points[5] = Point(light_tree_node->min_bounds[0], light_tree_node->max_bounds[1], light_tree_node->max_bounds[2]);
	bounding_points[6] = light_tree_node->max_bounds;
	bounding_points[7] = Point(light_tree_node->max_bounds[0], light_tree_node->min_bounds[1], light_tree_node->max_bounds[2]);

	float d;
	switch(vpl_type){
		case EPointEmitterVPL:
			d = std::max(min_dist, distToBox(shading_point_position, 
				light_tree_node->min_bounds, light_tree_node->max_bounds));
			geometric = 1.f / (d * d);
			break;
		case ESurfaceVPL:
			{
				d = std::max(min_dist, distToBox(shading_point_position,
					light_tree_node->min_bounds, light_tree_node->max_bounds));
				
				if(fabs(light_tree_node->cone_halfangle - M_PI) < std::numeric_limits<float>::epsilon()){
					geometric = 1.f / d * d;
				}
				else if((light_tree_node->max_bounds - light_tree_node->min_bounds).length() < 
					std::numeric_limits<float>::epsilon()){
					float degree = dot(normalize(shading_point_position - light_tree_node->min_bounds), light_tree_node->cone_ray);
					geometric = std::max(0.f, degree) / (d * d);
				}
				else{
					Vector3f unit_z(0.f, 0.f, 1.f);

					float angle = acos(dot(normalize(light_tree_node->cone_ray), unit_z));
					Vector3f axis = cross(light_tree_node->cone_ray, unit_z);

					//zero cross product means cone is pointing either in positive or negative z. Can set axis arbitrarily
					//as long as it is orthogonal to unit z
					if(axis.length() < 1e-20f){
						axis = Vector(1, 0, 0);
					}
					Transform rot = Transform::rotate(axis, angle);
					
					Point transformed_points[8];

					for(std::uint8_t i = 0; i < 8; ++i){
						Point p = Point(shading_point_position - bounding_points[i]);
						transformed_points[i] = rot.transformAffine(p);
					}

					float min_x = std::numeric_limits<float>::max();
					float min_y = std::numeric_limits<float>::max();
					float max_x = std::numeric_limits<float>::min();
					float max_y = std::numeric_limits<float>::min();
					float max_z = std::numeric_limits<float>::min();

					for(std::uint8_t i = 0; i < 8; ++i){
						min_x = std::min((float)transformed_points[i][0], min_x);
						min_y = std::min((float)transformed_points[i][1], min_y);
						max_x = std::max((float)transformed_points[i][0], max_x);
						max_y = std::max((float)transformed_points[i][1], max_y);
						max_z = std::max((float)transformed_points[i][2], max_z);
					}

					//box is behind light and won't get affected
					if(max_z <= 0.f){
						geometric = 0.f;
					}
					else{
						//intersecting
						if(min_x * max_x <= 0.f && min_y * max_y <= 0.f){
							geometric = 1.f;
						}
						else{
							float min_d = std::numeric_limits<float>::max();
							Point closest_point(0.f);

							std::pair<Point, Vector3f> lines[8];

							for(std::uint8_t i = 0; i < 4; ++i){
								std::uint32_t next = (i + 1) % 4;
								lines[i] = std::make_pair(transformed_points[i], transformed_points[next] - transformed_points[i]);
								lines[i + 4] = std::make_pair(transformed_points[i + 4], transformed_points[next + 4] - transformed_points[i + 4]);
							}

							for(std::uint8_t i = 0; i < 8; ++i){
								if(lines[i].second.length() <std::numeric_limits<float>::epsilon()){
									continue;
								}

								bool parallel = fabs(dot(normalize(lines[i].second), unit_z) - 1.f) > std::numeric_limits<float>::epsilon();
								//in the case where the line is parallel to the z axis, the closest point is one of the corners
								//which we will deal with when dealing with the corners
								if(parallel){
									continue;
								}

								Vector3f n1 = cross(lines[i].second, cross(unit_z, lines[i].second));
								Vector3f n2 = cross(unit_z, cross(lines[i].second, unit_z));

								Point on_unit_z = Point((dot(Vector3f(lines[i].first), n2) / dot(unit_z, n2)) * unit_z);

								float t = dot(Vector3f(-lines[i].first), n1) / dot(lines[i].second, n1);
								Point on_cube_edge = Point(lines[i].first + t * lines[i].second);

								float dist = (on_unit_z - on_cube_edge).length();
								if(dist < min_d && t >= 0.f && t <= 1.f){
									min_d = dist;
									closest_point = on_cube_edge;
								}
							}
							
							for(std::uint8_t i = 0; i < 8; ++i){
								//can just calculate distance in x, y since the projection to unit z is just setting them
								//to zero
								float dist = sqrt(transformed_points[i][0] * transformed_points[i][0] + 
									transformed_points[i][1] * transformed_points[i][1]);

								if(dist < min_d){
									min_d = dist;
									closest_point = transformed_points[i];
								}
							}

							float degree = dot(normalize(Vector3f(closest_point)), unit_z);
							geometric = std::max(0.f, degree);
						}
						
						geometric /= (d * d);
					}
				}
			}

			break;
		case EDirectionalEmitterVPL:
			geometric = 1.f;
			break;
		default:
			break;
	}

	//only dealing with lambertian diffuse right now
	float largest = std::numeric_limits<float>::min();
	for(std::uint8_t i = 0; i < 8; ++i){
		float dp = dot(normalize(bounding_points[i] - shading_point_position), shading_point_normal);
		largest = std::max(dp, largest);
	}

	float material = std::max(0.f, largest);

	return geometric * material * light_tree_node->emission_scale;
}

float calculateExactClusterContribution(Point shading_point_position, Normal shading_point_normal,
	LightTreeNode* light_tree_node, EVPLType vpl_type, float min_dist){
	
	if(light_tree_node == nullptr){
		return 0.f;
	}

	if(light_tree_node->left == nullptr && light_tree_node->right == nullptr){
		return calculateClusterContribution(shading_point_position, shading_point_normal, light_tree_node, 
			vpl_type, min_dist);
	}
	else{
		return calculateExactClusterContribution(shading_point_position, shading_point_normal,
			light_tree_node->left.get(), vpl_type, min_dist) + 
		calculateExactClusterContribution(shading_point_position, shading_point_normal,
			light_tree_node->right.get(), vpl_type, min_dist);
	}
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

	typedef std::tuple<LightTreeNode*, float> ClusterAndScore;
	//largest first, but front of the queue is last out
	auto comparator = [](ClusterAndScore l, ClusterAndScore r){
		return std::get<1>(l) < std::get<1>(r);
	};

	std::priority_queue<ClusterAndScore, std::vector<ClusterAndScore>, decltype(comparator)> pqueue(comparator);

	if(point_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(point_tree_root_.get(), std::numeric_limits<float>::max()));
	}

	if(oriented_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(oriented_tree_root_.get(), std::numeric_limits<float>::max()));
	}

	if(directional_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(directional_tree_root_.get(), std::numeric_limits<float>::max()));
	}
	
	while((pqueue.size() + lights.size()) < max_lights_ && pqueue.size() > 0){
		//std::cout << lights.size() << " " << pqueue.size() << std::endl;
		ClusterAndScore entry = pqueue.top();
		pqueue.pop();
		
		//if the worst node is below threshold, then all nodes must be
		if(std::get<1>(entry) < error_threshold_){
			break;
		}

		LightTreeNode* node = std::get<0>(entry);
		
		if(node->left == nullptr && node->right == nullptr){
			lights.push_back(std::get<0>(entry)->vpl);
			lights.back().P *= std::get<0>(entry)->emission_scale;
			continue;
		}

		if(node->left == nullptr || node->right == nullptr){
			std::cerr << "A node in the lighttree should always have 2 children" << std::endl;
			exit(0);
		}

		if(node->left.get() != nullptr){
			float rad = calculateClusterContribution(its.p, its.geoFrame.n, node->left.get(), node->left->vpl.type, min_dist_);
			float actual_rad = calculateExactClusterContribution(its.p, its.geoFrame.n, node->left.get(), node->left->vpl.type, min_dist_);
			
			if(actual_rad > 0.f){
				float error = fabs(rad - actual_rad) / actual_rad;
				pqueue.push(std::make_tuple(node->left.get(), error));
			}
		}

		if(node->right.get() != nullptr){
			float rad = calculateClusterContribution(its.p, its.geoFrame.n, node->right.get(), node->right->vpl.type, min_dist_);
			float actual_rad = calculateExactClusterContribution(its.p, its.geoFrame.n, node->right.get(), node->right->vpl.type, min_dist_);
			if(actual_rad > 0.f){
				float error = fabs(rad - actual_rad) / actual_rad;
				pqueue.push(std::make_tuple(node->right.get(), error));
			}
		}
	}

	while(pqueue.size() > 0 && lights.size() < max_lights_){
		auto entry = pqueue.top();
		lights.push_back(std::get<0>(entry)->vpl);
		lights.back().P *= std::get<0>(entry)->emission_scale;
		//std::cout << std::get<1>(entry) << " ";
		pqueue.pop();
	}
	//std::cout << std::endl;
	//should be moved so it should be ok
	return lights;
}

std::vector<VPL> LightTree::getClusteringForPoints(const std::vector<Intersection>& points){
	std::vector<VPL> lights;

	typedef std::tuple<LightTreeNode*, float> ClusterAndScore;
	//largest first, but front of the queue is last out
	auto comparator = [](ClusterAndScore l, ClusterAndScore r){
		return std::get<1>(l) < std::get<1>(r);
	};

	std::priority_queue<ClusterAndScore, std::vector<ClusterAndScore>, decltype(comparator)> pqueue(comparator);

	if(point_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(point_tree_root_.get(), std::numeric_limits<float>::max()));
	}

	if(oriented_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(oriented_tree_root_.get(), std::numeric_limits<float>::max()));
	}

	if(directional_tree_root_ != nullptr){
		pqueue.push(std::make_tuple(directional_tree_root_.get(), std::numeric_limits<float>::max()));
	}
	
	while((pqueue.size() + lights.size()) < max_lights_ && pqueue.size() > 0){
		//std::cout << lights.size() << " " << pqueue.size() << std::endl;
		ClusterAndScore entry = pqueue.top();
		pqueue.pop();

		LightTreeNode* node = std::get<0>(entry);
		
		if(node->left == nullptr && node->right == nullptr){
			lights.push_back(std::get<0>(entry)->vpl);
			lights.back().P *= std::get<0>(entry)->emission_scale;
			continue;
		}

		if(node->left == nullptr || node->right == nullptr){
			std::cerr << "A node in the lighttree should always have 2 children" << std::endl;
			exit(0);
		}

		if(node->left.get() != nullptr){
			float rad = 0.f;
			float actual_rad = 0.f;
			for(std::uint32_t i = 0; i < points.size(); ++i){
				rad += calculateClusterContribution(points[i].p, points[i].geoFrame.n, node->left.get(), node->left->vpl.type, min_dist_);
				actual_rad = calculateExactClusterContribution(points[i].p, points[i].geoFrame.n, node->left.get(), node->left->vpl.type, min_dist_);
			}

			if(actual_rad > 0.f){
				float error = fabs(rad - actual_rad) / actual_rad;
				pqueue.push(std::make_tuple(node->left.get(), error));
			}
		}

		if(node->right.get() != nullptr){
			float rad = 0.f;
			float actual_rad = 0.f;
			for(std::uint32_t i = 0; i < points.size(); ++i){
				rad += calculateClusterContribution(points[i].p, points[i].geoFrame.n, node->right.get(), node->right->vpl.type, min_dist_);
				actual_rad += calculateExactClusterContribution(points[i].p, points[i].geoFrame.n, node->right.get(), node->right->vpl.type, min_dist_);
			}

			if(actual_rad > 0.f){
				float error = fabs(rad - actual_rad) / actual_rad;
				pqueue.push(std::make_tuple(node->right.get(), error));
			}
		}
	}

	while(pqueue.size() > 0 && lights.size() < max_lights_){
		auto entry = pqueue.top();
		lights.push_back(std::get<0>(entry)->vpl);
		lights.back().P *= std::get<0>(entry)->emission_scale;
		pqueue.pop();
	}

	return lights;
}

MTS_NAMESPACE_END