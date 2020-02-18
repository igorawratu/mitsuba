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
#include <unordered_map>
#include "common.h"
#include <eigen3/Eigen/Dense>

MTS_NAMESPACE_BEGIN

const std::array<std::function<bool(float, const Intersection&)>, 6> divider_comp {
    [](float value, const Intersection& its){
        return value < its.p.x;
    },
    [](float value, const Intersection& its){
        return value < its.p.y;
    },
    [](float value, const Intersection& its){
        return value < its.p.z;
    },
    [](float value, const Intersection& its){
        return value < its.shFrame.n.x;
    },
    [](float value, const Intersection& its){
        return value < its.shFrame.n.y;
    },
    [](float value, const Intersection& its){
        return value < its.shFrame.n.z;
    }
};

const std::array<std::function<bool(const VPL&, const VPL&)>, 6> divider_sorter{
    [](const VPL& lhs, const VPL& rhs){
        return lhs.its.p.x < rhs.its.p.x;
    },
    [](const VPL& lhs, const VPL& rhs){
        return lhs.its.p.y < rhs.its.p.y;
    },
    [](const VPL& rhs, const VPL& lhs){
        return lhs.its.p.z < rhs.its.p.z;
    },
    [](const VPL& lhs, const VPL& rhs){
        return lhs.its.shFrame.n.x < rhs.its.shFrame.n.x;
    },
    [](const VPL& lhs, const VPL& rhs){
        return lhs.its.shFrame.n.y < rhs.its.shFrame.n.y;
    },
    [](const VPL& lhs, const VPL& rhs){
        return lhs.its.shFrame.n.z < rhs.its.shFrame.n.z;
    }
};

void divideByGreatestDim(const std::vector<VPL>& vpls,
	std::vector<VPL>& left, std::vector<VPL>& right, std::uint32_t min_size,
	float norm_scale, EVPLType vpl_type){
	left.clear();
	right.clear();

	float maxf = std::numeric_limits<float>::max();
	Vector3f min_pos(maxf, maxf, maxf), max_pos(-maxf, -maxf, -maxf);

	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		const VPL& vpl = vpls[i];

		min_pos.x = std::min(vpl.its.p.x, min_pos.x);
		min_pos.y = std::min(vpl.its.p.y, min_pos.y);
		min_pos.z = std::min(vpl.its.p.z, min_pos.z);
		max_pos.x = std::max(vpl.its.p.x, max_pos.x);
		max_pos.y = std::max(vpl.its.p.y, max_pos.y);
		max_pos.z = std::max(vpl.its.p.z, max_pos.z);
	}

	bool has_directional = vpl_type != EPointEmitterVPL;
	bool has_positional = vpl_type != EDirectionalEmitterVPL;
	bool consider_norm = (max_pos - min_pos).length() < norm_scale || vpl_type == EDirectionalEmitterVPL;
	norm_scale = consider_norm ? norm_scale : 0.f;

	Eigen::MatrixXf mean = Eigen::MatrixXf::Zero(6, 1);
	Eigen::MatrixXf cov_mat = Eigen::MatrixXf::Zero(6, 6);
	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		const VPL& vpl = vpls[i];
		Eigen::MatrixXf curr_pt(6, 1);
		curr_pt(0, 0) = has_positional ? vpl.its.p.x : 0.f;
		curr_pt(1, 0) = has_positional ? vpl.its.p.y : 0.f;
		curr_pt(2, 0) = has_positional ? vpl.its.p.z : 0.f;
		curr_pt(3, 0) = has_directional ? vpl.its.shFrame.n.x * norm_scale : 0.f;
		curr_pt(4, 0) = has_directional ? vpl.its.shFrame.n.y * norm_scale : 0.f;
		curr_pt(5, 0) = has_directional ? vpl.its.shFrame.n.z * norm_scale : 0.f;

		mean += curr_pt;

		cov_mat += curr_pt * curr_pt.transpose();
	}

	mean /= vpls.size();
	cov_mat /= vpls.size() - 1;

	auto svd = cov_mat.jacobiSvd(Eigen::ComputeThinU);
	Eigen::VectorXf axis = svd.matrixU().col(0);

	float proj_max = -std::numeric_limits<float>::max();
	float proj_min = std::numeric_limits<float>::max();

	std::vector<std::pair<float, std::uint32_t>> projections(vpls.size());
	float proj_mean = 0.f;
	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		const VPL& vpl = vpls[i];
		Eigen::VectorXf curr_pt(6);
		curr_pt(0) = has_positional ? vpl.its.p.x : 0.f;
		curr_pt(1) = has_positional ? vpl.its.p.y : 0.f;
		curr_pt(2) = has_positional ? vpl.its.p.z : 0.f;
		curr_pt(3) = has_directional ? vpl.its.shFrame.n.x * norm_scale : 0.f;
		curr_pt(4) = has_directional ? vpl.its.shFrame.n.y * norm_scale : 0.f;
		curr_pt(5) = has_directional ? vpl.its.shFrame.n.z * norm_scale : 0.f;

		float proj_rat = curr_pt.dot(axis);
		proj_mean += proj_rat;
		projections[i] = std::make_pair(proj_rat, i);

		proj_max = std::max(proj_rat, proj_max);
		proj_min = std::min(proj_rat, proj_min);
	}
	proj_mean /= vpls.size();

	float midpoint = (proj_max + proj_min) / 2.f;
	for(std::uint32_t i = 0; i < projections.size(); ++i){
		if(projections[i].first < proj_mean){
			left.push_back(vpls[projections[i].second]);
		}
		else{
			right.push_back(vpls[projections[i].second]);
		}
	}

	//the same point or points really close, no real need to sort, in reality i don't believe this degenerate
	//should really occur unless sampling from point lights
	if(left.size() < min_size || right.size() < min_size){
		left.clear();
		right.clear();
		std::uint32_t half = projections.size() / 2;
		for(std::uint32_t i = 0; i < projections.size(); ++i){
			if(i < half){
				left.push_back(vpls[projections[i].second]);
			}
			else{
				right.push_back(vpls[projections[i].second]);
			}
		}
	}
}

std::unique_ptr<LightTreeNode> createLightTree(const std::vector<VPL>& vpls, EVPLType vpl_type, float min_dist, std::mt19937& rng) {
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
		nodes[current_level][i]->emission_scale = nodes[current_level][i]->vpl.P.getLuminance();
		nodes[current_level][i]->vpl.P /= nodes[current_level][i]->emission_scale;
	}

	typedef std::tuple<float, std::uint32_t, std::uint32_t, Vector3f, float> SimEntry;

	auto comparator = [](SimEntry l, SimEntry r){
		return std::get<0>(l) < std::get<0>(r);
	};

	std::uniform_real_distribution<float> gen(0, 1);

	while(true){
		if(nodes[current_level].size() == 1){
			break;
		}

		std::vector<SimEntry> similarity_matrix;

		std::mutex mat_insert_mut;
		#pragma omp parallel for
		for(size_t i = 0; i < nodes[current_level].size(); ++i){
			#pragma omp parallel for
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

					Vector3 axis = cross(ray1, ray2);
					float angle = atan2(axis.length(), dot(ray1, ray2));

					float max_child_half = std::max(ha1, ha2);
					float min_child_half = std::min(ha1, ha2);

					//fully overlapping cones
					if(max_child_half - angle > min_child_half || (ray1 - ray2).length() < std::numeric_limits<float>::epsilon()
						|| M_PI - max_child_half < 0.0001f){
						union_angle_span = max_child_half;
						bcone = ha1 > ha2 ? ray1 : ray2;
					}
					else{
						union_angle_span = angle + ha1 + ha2;
						union_angle_span = union_angle_span > M_PI ? 2 * M_PI : union_angle_span;
						union_angle_span = std::min(union_angle_span / 2.f, M_PI);

						if(axis.length() < std::numeric_limits<float>::epsilon()){
							axis = cross(ray1, Vector3f(gen(rng), gen(rng), gen(rng)));
						}

						Point new_coneray = Transform::rotate(axis, -ha1).transformAffine(Point(ray1));
						bcone = Vector(Transform::rotate(axis, union_angle_span).transformAffine(new_coneray));
					}

					d += union_angle_span * union_angle_span * min_dist * 2.f;
				}

				{
					std::lock_guard<std::mutex> lock(mat_insert_mut);
					similarity_matrix.push_back(std::make_tuple(d, i, j, bcone, union_angle_span));
				}
				
			}
		}

		std::make_heap(similarity_matrix.begin(), similarity_matrix.end(), comparator);

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

			std::pop_heap(similarity_matrix.begin(), similarity_matrix.end(), comparator);
			SimEntry entry = similarity_matrix.back();
			similarity_matrix.pop_back();

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
				//float ratio = c1_intensity / (c1_intensity + c2_intensity);

				float sample = gen(rng);

				if(c1_intensity > 0.f || c2_intensity > 0.f){
					float lc_rat = c1_intensity / (c1_intensity + c2_intensity);
					if (sample < lc_rat) {
						nodes[next_level].back()->vpl = c1->vpl;
					}
					else {
						nodes[next_level].back()->vpl = c2->vpl;
					}
					
					//nodes[next_level].back()->vpl.P = c1->vpl.P * lc_rat + c2->vpl.P * (1.f - lc_rat);
					//nodes[next_level].back()->vpl.P /= nodes[next_level].back()->vpl.P.getLuminance();

					nodes[next_level].back()->cone_ray = std::get<3>(entry);
					nodes[next_level].back()->cone_halfangle = std::get<4>(entry);
					/*nodes[next_level].back()->vpl.its.p = Point(ratio * Vector3f(c1->vpl.its.p) + (1.f - ratio) * Vector3f(c2->vpl.its.p));
					Vector3f new_n = ratio * Vector3f(c1->vpl.its.shFrame.n) + (1.f - ratio) * Vector3f(c2->vpl.its.shFrame.n);
					if(new_n.length() < 1e-10f){
						new_n = cross(Vector3f(c1->vpl.its.shFrame.n), Vector3f(gen(rng), gen(rng), gen(rng)));
					}
					nodes[next_level].back()->vpl.its.shFrame = Frame(normalize(new_n));
					nodes[next_level].back()->vpl.P = ratio * c1->vpl.P + (1.f - ratio) * c2->vpl.P;*/
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

std::unique_ptr<LightTreeNode> tdCreateLightTree(const std::vector<VPL>& vpls, EVPLType vpl_type, float min_dist, std::uint32_t bottom_up_thresh, std::mt19937& rng){
	if(vpls.size() <= bottom_up_thresh){
		return createLightTree(vpls, vpl_type, min_dist, rng);
	}

	std::vector<VPL> left_vpls;
	std::vector<VPL> right_vpls;
	divideByGreatestDim(vpls, left_vpls, right_vpls, bottom_up_thresh, min_dist / 10.f, vpl_type);

	auto lc = tdCreateLightTree(left_vpls, vpl_type, min_dist, bottom_up_thresh, rng);
	auto rc = tdCreateLightTree(right_vpls, vpl_type, min_dist, bottom_up_thresh, rng);

	std::uniform_real_distribution<float> gen(0, 1);

	std::unique_ptr<LightTreeNode> curr_node(new LightTreeNode());
	float union_angle_span = 0.f;
	Vector3f bcone(0.f);

	if(vpl_type != EDirectionalEmitterVPL){
		curr_node->min_bounds = Point(std::min(lc->min_bounds[0], rc->min_bounds[0]), 
			std::min(lc->min_bounds[1], rc->min_bounds[1]), 
			std::min(lc->min_bounds[2], rc->min_bounds[2]));
		curr_node->max_bounds = Point(std::max(lc->max_bounds[0], rc->max_bounds[0]), 
			std::max(lc->max_bounds[1], rc->max_bounds[1]), 
			std::max(lc->max_bounds[2], rc->max_bounds[2]));
	}

	if(vpl_type != EPointEmitterVPL){
		Vector3 ray1 = lc->cone_ray;
		Vector3 ray2 = rc->cone_ray;
		float ha1 = lc->cone_halfangle;
		float ha2 = rc->cone_halfangle;

		float angle = atan2(cross(ray1, ray2).length(), dot(ray1, ray2));
		float max_child_half = std::max(ha1, ha2);
		float min_child_half = std::min(ha1, ha2);

		if(max_child_half - angle > min_child_half || (ray1 - ray2).length() < std::numeric_limits<float>::epsilon()
			|| M_PI - max_child_half < 0.0001f){
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
	}

	curr_node->cone_ray = bcone;
	curr_node->cone_halfangle = union_angle_span;
	
	float lc_intensity = lc->emission_scale;
	float rc_intensity = rc->emission_scale;

	float sample = gen(rng);
	float lc_rat = lc_intensity / (lc_intensity + rc_intensity);

	if(lc_intensity > 0.f || rc_intensity > 0.f){
		if (sample < lc_rat) {
			curr_node->vpl = lc->vpl;
		}
		else {
			curr_node->vpl = rc->vpl;
		}
	}

	//curr_node->vpl.P = lc->vpl.P * lc_rat + rc->vpl.P * (1.f - lc_rat);
	//curr_node->vpl.P /= curr_node->vpl.P.getLuminance();

	curr_node->emission_scale = lc_intensity + rc_intensity;

	curr_node->left = std::move(lc);
	curr_node->right = std::move(rc);

	return curr_node;
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

float calculateClusterEstimate(Scene* scene, Point shading_point_position, Normal shading_point_normal, LightTreeNode* node,
	float min_dist){
	Intersection its;
	its.p = shading_point_position;
	/*if(!sampleVisibility(scene, its, node->vpl, min_dist)){
		return 0.f;
	}*/

	float d;
	float geometric = 0.f;

	switch(node->vpl.type){
		case EPointEmitterVPL:
			d = std::max(min_dist, (node->vpl.its.p - shading_point_position).length());
			geometric = 1.f / (d * d);
			break;
		case ESurfaceVPL:
			{
				Vector3 dir = shading_point_position - node->vpl.its.p;
				d = std::max(min_dist, dir.length());
				geometric = dot(normalize(dir), node->vpl.its.shFrame.n) / (d * d);
			}
			break;
		case EDirectionalEmitterVPL:
			geometric = 1.f;
			break;
		default:
			break;
	}

	float material = 0.f;

	if(node->vpl.type != EDirectionalEmitterVPL){
		material = dot(normalize(node->vpl.its.p - shading_point_position), shading_point_normal);
	}
	else{
		material = dot(-node->vpl.its.shFrame.n, shading_point_normal);
	}

	return material * geometric * node->emission_scale;
}

float calculateClusterBounds(Point shading_point_position, Normal shading_point_normal,
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
					Vector3f normalized_cone_ray = normalize(light_tree_node->cone_ray);
					float angle = atan2(cross(normalized_cone_ray, unit_z).length(), dot(normalized_cone_ray, unit_z));
					Vector3f axis = cross(light_tree_node->cone_ray, unit_z);

					//zero cross product means cone is pointing either in positive or negative z. Can set axis arbitrarily
					//as long as it is orthogonal to unit z
					if(axis.lengthSquared() < 1e-20f){
						axis = Vector(1, 0, 0);
					}
					Transform rot = Transform::rotate(axis, angle);
					
					Point transformed_points[8];

					for(std::uint8_t i = 0; i < 8; ++i){
						Point p = Point(shading_point_position - bounding_points[i]);
						transformed_points[i] = rot.transformAffine(p);
					}

					float min_x2 = std::numeric_limits<float>::max();
					float min_y2 = std::numeric_limits<float>::max();
					float max_x2 = std::numeric_limits<float>::min();
					float max_y2 = std::numeric_limits<float>::min();
					float max_z = std::numeric_limits<float>::min();

					for(std::uint8_t i = 0; i < 8; ++i){
						float x2 = transformed_points[i][0] * transformed_points[i][0];
						float y2 = transformed_points[i][1] * transformed_points[i][1];
						min_x2 = std::min(x2, min_x2);
						min_y2 = std::min(y2, min_y2);
						max_x2 = std::max(x2, max_x2);
						max_y2 = std::max(y2, max_y2);
						max_z = std::max((float)transformed_points[i][2], max_z);
					}

					float cos_theta = max_z > 0.f ? max_z / sqrt(min_x2 + min_y2 + max_z * max_z) :
						max_z / sqrt(max_x2 + max_y2 + max_z * max_z);
					float degree = acos(cos_theta);

					degree = std::max(0.f, degree - light_tree_node->cone_halfangle);

					geometric = std::max(0., cos(degree));
				}
				
				geometric /= (d * d);
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

	if(vpl_type != EDirectionalEmitterVPL){
		for(std::uint8_t i = 0; i < 8; ++i){
			float dp = dot(normalize(bounding_points[i] - shading_point_position), shading_point_normal);
			float angle = acos(dp);
			float dif = std::max(0.f, angle - light_tree_node->cone_halfangle);
			largest = std::max(float(cos(dif)), largest);
		}
		largest = std::max(0.f, largest);
	}
	else{
		float angle = acos(dot(-light_tree_node->cone_ray, shading_point_normal));
		float dif = std::max(0.f, angle - light_tree_node->cone_halfangle);
		largest = std::max(0., cos(dif));
	}

	float material = std::max(0.f, largest);

	return geometric * material * light_tree_node->emission_scale;
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

	std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	point_tree_root_ = tdCreateLightTree(point_vpls_, EPointEmitterVPL, min_dist_, 1000, rng);
	oriented_tree_root_ = tdCreateLightTree(oriented_vpls_, ESurfaceVPL, min_dist_, 1000, rng);
	directional_tree_root_ = tdCreateLightTree(directional_vpls_, EDirectionalEmitterVPL, min_dist_, 1000, rng);
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

	std::unordered_map<LightTreeNode*, float> contribution_cache;
	contribution_cache.reserve((point_vpls_.size() + directional_vpls_.size() + oriented_vpls_.size()) * 3);

	while((pqueue.size() + lights.size()) < max_lights_ && pqueue.size() > 0){
		ClusterAndScore entry = pqueue.top();
		pqueue.pop();
		
		//if the worst node is below threshold, then all nodes must be
		//if(std::get<1>(entry) < error_threshold_){
		//	break;
		//}

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
			float rad = calculateClusterBounds(its.p, its.shFrame.n, node->left.get(), node->left->vpl.type, min_dist_);
			pqueue.push(std::make_tuple(node->left.get(), rad));
		}

		if(node->right.get() != nullptr){
			float rad = calculateClusterBounds(its.p, its.shFrame.n, node->right.get(), node->right->vpl.type, min_dist_);
			pqueue.push(std::make_tuple(node->right.get(), rad));
		}
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

std::vector<VPL> LightTree::getClusteringForPoints(Scene* scene, const std::vector<Intersection>& points){
	std::vector<VPL> lightcut;

	typedef std::tuple<LightTreeNode*, float> ClusterAndScore;
	//largest first, but front of the queue is last out
	auto comparator = [](ClusterAndScore l, ClusterAndScore r){
		return std::get<1>(l) < std::get<1>(r);
	};

	std::priority_queue<ClusterAndScore, std::vector<ClusterAndScore>, decltype(comparator)> pqueue(comparator);
	std::unordered_map<LightTreeNode*, float> contribution_cache;
	std::uint32_t total_lights = point_vpls_.size() + directional_vpls_.size() + oriented_vpls_.size();

	Point slice_centroid_pos(0.f);
	Vector3f slice_centroid_normal(0.f);
	for(std::uint32_t i = 0; i < points.size(); ++i){
		slice_centroid_pos += points[i].p;
		slice_centroid_normal += Vector3f(points[i].shFrame.n);
	}
	slice_centroid_pos /= points.size();
	slice_centroid_normal = normalize(slice_centroid_normal);

	/*std::vector<VPL> lights;
	if(point_tree_root_ != nullptr){
		float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), point_tree_root_.get(), 
			point_tree_root_->vpl.type, min_dist_);
		pqueue.push(std::make_tuple(point_tree_root_.get(), rad));
	}

	if(oriented_tree_root_ != nullptr){
		float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), oriented_tree_root_.get(), 
			oriented_tree_root_->vpl.type, min_dist_);
		pqueue.push(std::make_tuple(oriented_tree_root_.get(), rad));
	}

	if(directional_tree_root_ != nullptr){
		float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), directional_tree_root_.get(), 
			directional_tree_root_->vpl.type, min_dist_);
		pqueue.push(std::make_tuple(directional_tree_root_.get(), rad));
	}

	while((pqueue.size() + lights.size()) < max_lights_ && pqueue.size() > 0){
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
			float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), node->left.get(), 
				node->left->vpl.type, min_dist_);

			pqueue.push(std::make_tuple(node->left.get(), rad));
		}

		if(node->right.get() != nullptr){
			float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), 
				node->right.get(), node->right->vpl.type, min_dist_);

			pqueue.push(std::make_tuple(node->right.get(), rad));
		}
	}

	while(pqueue.size() > 0 && lights.size() < max_lights_){
		auto entry = pqueue.top();
		lights.push_back(std::get<0>(entry)->vpl);
		lights.back().P *= std::get<0>(entry)->emission_scale;
		pqueue.pop();
	}

	while(pqueue.size() > 0){
		pqueue.pop();
	}

	lightcut.insert(lightcut.end(), lights.begin(), lights.end());*/

	if(point_tree_root_ != nullptr){
		std::vector<VPL> lights;
		pqueue.push(std::make_tuple(point_tree_root_.get(), std::numeric_limits<float>::max()));
		float ratio = float(point_vpls_.size()) / float(total_lights);
		std::uint32_t lights_to_extract = ratio * max_lights_ + 0.5f;

		while((pqueue.size() + lights.size()) < lights_to_extract && pqueue.size() > 0){
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
				float diag = (node->left->max_bounds - node->left->min_bounds).length();
				float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), node->left.get(), 
					node->left->vpl.type, min_dist_);// * diag;

				pqueue.push(std::make_tuple(node->left.get(), rad));
			}

			if(node->right.get() != nullptr){
				float diag = (node->right->max_bounds - node->right->min_bounds).length();
				float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), 
					node->right.get(), node->right->vpl.type, min_dist_);// * diag;

				pqueue.push(std::make_tuple(node->right.get(), rad));
			}
		}

		while(pqueue.size() > 0 && lights.size() < lights_to_extract){
			auto entry = pqueue.top();
			lights.push_back(std::get<0>(entry)->vpl);
			lights.back().P *= std::get<0>(entry)->emission_scale;
			pqueue.pop();
		}

		while(pqueue.size() > 0){
			pqueue.pop();
		}

		lightcut.insert(lightcut.end(), lights.begin(), lights.end());
	}

	if(oriented_tree_root_ != nullptr){
		std::vector<VPL> lights;
		std::vector<float> errs;
		pqueue.push(std::make_tuple(oriented_tree_root_.get(), std::numeric_limits<float>::max()));
		float ratio = float(oriented_vpls_.size()) / float(total_lights);
		std::uint32_t lights_to_extract = ratio * max_lights_ + 0.5f;

		while((pqueue.size() + lights.size()) < lights_to_extract && pqueue.size() > 0){
			ClusterAndScore entry = pqueue.top();
			pqueue.pop();

			LightTreeNode* node = std::get<0>(entry);
			
			if(node->left == nullptr && node->right == nullptr){
				lights.push_back(std::get<0>(entry)->vpl);
				lights.back().P *= std::get<0>(entry)->emission_scale;
				errs.push_back(std::get<1>(entry));
				continue;
			}

			if(node->left == nullptr || node->right == nullptr){
				std::cerr << "A node in the lighttree should always have 2 children" << std::endl;
				exit(0);
			}

			if(node->left.get() != nullptr){
				float diag = (node->left->max_bounds - node->left->min_bounds).length();
				float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), 
					node->left.get(), node->left->vpl.type, min_dist_);// * (node->left->cone_halfangle * min_dist_ + diag);

				pqueue.push(std::make_tuple(node->left.get(), rad));
			}

			if(node->right.get() != nullptr){
				float diag = (node->right->max_bounds - node->right->min_bounds).length();
				float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), 
					node->right.get(), node->right->vpl.type, min_dist_);// * (node->right->cone_halfangle * min_dist_ + diag);

				pqueue.push(std::make_tuple(node->right.get(), rad));
			}
		}

		while(pqueue.size() > 0 && lights.size() < lights_to_extract){
			auto entry = pqueue.top();
			lights.push_back(std::get<0>(entry)->vpl);
			lights.back().P *= std::get<0>(entry)->emission_scale;
			errs.push_back(std::get<1>(entry));
			pqueue.pop();
		}

		while(pqueue.size() > 0){
			pqueue.pop();
		}

		lightcut.insert(lightcut.end(), lights.begin(), lights.end());
	}

	if(directional_tree_root_ != nullptr){
		std::vector<VPL> lights;
		std::vector<float> errs;
		pqueue.push(std::make_tuple(directional_tree_root_.get(), std::numeric_limits<float>::max()));

		float ratio = float(directional_vpls_.size()) / float(total_lights);
		std::uint32_t lights_to_extract = ratio * max_lights_ + 0.5f;

		while((pqueue.size() + lights.size()) < lights_to_extract && pqueue.size() > 0){
			ClusterAndScore entry = pqueue.top();
			pqueue.pop();

			LightTreeNode* node = std::get<0>(entry);
			
			if(node->left == nullptr && node->right == nullptr){
				lights.push_back(std::get<0>(entry)->vpl);
				lights.back().P *= std::get<0>(entry)->emission_scale;
				errs.push_back(std::get<1>(entry));
				continue;
			}

			if(node->left == nullptr || node->right == nullptr){
				std::cerr << "A node in the lighttree should always have 2 children" << std::endl;
				exit(0);
			}

			if(node->left.get() != nullptr){
				float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), 
					node->left.get(), node->left->vpl.type, min_dist_);// * node->left->cone_halfangle;

				pqueue.push(std::make_tuple(node->left.get(), rad));
			}

			if(node->right.get() != nullptr){
				float rad = calculateClusterBounds(slice_centroid_pos, Normal(slice_centroid_normal), 
					node->right.get(), node->right->vpl.type, min_dist_);// * node->right->cone_halfangle;

				pqueue.push(std::make_tuple(node->right.get(), rad));
			}
		}

		while(pqueue.size() > 0 && lights.size() < lights_to_extract){
			auto entry = pqueue.top();
			lights.push_back(std::get<0>(entry)->vpl);
			lights.back().P *= std::get<0>(entry)->emission_scale;
			errs.push_back(std::get<1>(entry));
			pqueue.pop();
		}

		while(pqueue.size() > 0){
			pqueue.pop();
		}

		lightcut.insert(lightcut.end(), lights.begin(), lights.end());
	}

	while(lightcut.size() > max_lights_)
		lightcut.pop_back();

	return lightcut;
}

MTS_NAMESPACE_END