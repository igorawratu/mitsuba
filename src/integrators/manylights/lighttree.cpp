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

void divideByGreatestDimLargest(const std::vector<VPL>& vpls, std::vector<VPL>& left, std::vector<VPL>& right, std::uint32_t min_size,
	float norm_scale, EVPLType vpl_type){
	left.clear();
	right.clear();

	float maxf = std::numeric_limits<float>::max();
	Vector3f min_pos(maxf, maxf, maxf), max_pos(-maxf, -maxf, -maxf);
	Vector3f min_normal(maxf, maxf, maxf), max_normal(-maxf, -maxf, -maxf);

	float norm_factor = vpl_type != EPointEmitterVPL ? 1.f : 0.f;
	float pos_factor = vpl_type != EDirectionalEmitterVPL ? 1.f : 0.f;

	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		min_pos.x = std::min(vpls[i].its.p.x, min_pos.x);
		min_pos.y = std::min(vpls[i].its.p.y, min_pos.y);
		min_pos.z = std::min(vpls[i].its.p.z, min_pos.z);
		max_pos.x = std::max(vpls[i].its.p.x, max_pos.x);
		max_pos.y = std::max(vpls[i].its.p.y, max_pos.y);
		max_pos.z = std::max(vpls[i].its.p.z, max_pos.z);

		min_normal.x = std::min(vpls[i].its.shFrame.n.x, min_normal.x);
		min_normal.y = std::min(vpls[i].its.shFrame.n.y, min_normal.y);
		min_normal.z = std::min(vpls[i].its.shFrame.n.z, min_normal.z);
		max_normal.x = std::max(vpls[i].its.shFrame.n.x, max_normal.x);
		max_normal.y = std::max(vpls[i].its.shFrame.n.y, max_normal.y);
		max_normal.z = std::max(vpls[i].its.shFrame.n.z, max_normal.z);
	}

	std::array<std::pair<std::uint8_t, float>, 6> ranges;
	ranges[0] = std::make_pair(0, (max_pos.x - min_pos.x) * pos_factor);
	ranges[1] = std::make_pair(1, (max_pos.y - min_pos.y) * pos_factor);
	ranges[2] = std::make_pair(2, (max_pos.z - min_pos.z) * pos_factor);
	ranges[3] = std::make_pair(3, (max_normal.x - min_normal.x) * norm_scale * norm_factor);
	ranges[4] = std::make_pair(4, (max_normal.y - min_normal.y) * norm_scale * norm_factor);
	ranges[5] = std::make_pair(5, (max_normal.z - min_normal.z) * norm_scale * norm_factor);

	std::array<float, 6> midpoints;
	midpoints[0] = (max_pos.x + min_pos.x) / 2.f;
	midpoints[1] = (max_pos.y + min_pos.y) / 2.f;
	midpoints[2] = (max_pos.z + min_pos.z) / 2.f;
	midpoints[3] = (max_normal.x + min_normal.x) / 2.f;
	midpoints[4] = (max_normal.y + min_normal.y) / 2.f;
	midpoints[5] = (max_normal.z + min_normal.z) / 2.f;

	std::sort(ranges.begin(), ranges.end(), 
		[](const std::pair<std::uint8_t, float>& lhs, const std::pair<std::uint8_t, float>& rhs){
			return lhs.second > rhs.second;
		});

	for(std::uint32_t i = 0; i < vpls.size(); ++i){
		if(divider_comp[ranges[0].first](midpoints[ranges[0].first], vpls[i].its)){
			right.push_back(vpls[i]);
		}
		else{
			left.push_back(vpls[i]);
		}
	}

	if(left.size() < min_size || right.size() < min_size){
		auto vpl_sorted = vpls;
		std::sort(vpl_sorted.begin(), vpl_sorted.end(), divider_sorter[ranges[0].first]);
		left.clear();
		right.clear();
		std::uint32_t midpoint = vpl_sorted.size() / 2;
		left.insert(left.end(), vpl_sorted.begin(), vpl_sorted.begin() + midpoint);
		right.insert(right.end(), vpl_sorted.begin() + midpoint, vpl_sorted.end());
	}
}

void divideByGreatestDimPA(const std::vector<VPL>& vpls,
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
			Vector3f cone_axis(vpls[i].its.shFrame.n);
			nodes[current_level][i]->bcone = DirConef(cone_axis);
		}
		
		nodes[current_level][i]->vpl = vpls[i];
		nodes[current_level][i]->emission_scale = nodes[current_level][i]->vpl.P.getLuminance();
		nodes[current_level][i]->vpl.P /= nodes[current_level][i]->emission_scale;
		nodes[current_level][i]->num_children = 0;
	}

	typedef std::tuple<float, std::uint32_t, std::uint32_t, DirConef> SimEntry;

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
				DirConef union_cone;

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
					DirConef cone1 = nodes[current_level][i]->bcone;
					DirConef cone2 = nodes[current_level][j]->bcone;
					union_cone = DirConef::Union(cone1, cone2);
					float angle = union_cone.GetAngleCos();

					d += angle * angle * min_dist * 2.f;
				}

				{
					std::lock_guard<std::mutex> lock(mat_insert_mut);
					similarity_matrix.push_back(std::make_tuple(d, i, j, union_cone));
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


					nodes[next_level].back()->bcone = std::get<3>(entry);
				}

				nodes[next_level].back()->emission_scale = c1_intensity + c2_intensity;

				nodes[next_level].back()->num_children = 2 + nodes[current_level][id1]->num_children + nodes[current_level][id2]->num_children;

				nodes[next_level].back()->left = std::move(nodes[current_level][id1]);
				nodes[next_level].back()->right = std::move(nodes[current_level][id2]);

				float d = (nodes[next_level].back()->left->vpl.its.p - nodes[next_level].back()->right->vpl.its.p).length();

				nodes[next_level].back()->vpl.radius = std::min(
					nodes[next_level].back()->left->vpl.radius + nodes[next_level].back()->right->vpl.radius,
					d / 2 + std::max(nodes[next_level].back()->left->vpl.radius, nodes[next_level].back()->right->vpl.radius));

				nodes[next_level].back()->vpl.radius *= nodes[next_level].back()->vpl.radius;

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

std::unique_ptr<LightTreeNode> tdCreateLightTree(const std::vector<VPL>& vpls, EVPLType vpl_type, float min_dist, std::uint32_t bottom_up_thresh, 
	std::mt19937& rng, bool use_pa_div){
	if(vpls.size() <= bottom_up_thresh){
		return createLightTree(vpls, vpl_type, min_dist, rng);
	}

	std::vector<VPL> left_vpls;
	std::vector<VPL> right_vpls;
	if(use_pa_div){
		divideByGreatestDimPA(vpls, left_vpls, right_vpls, bottom_up_thresh, min_dist / 10.f, vpl_type);
	}
	else{
		divideByGreatestDimLargest(vpls, left_vpls, right_vpls, bottom_up_thresh, min_dist / 10.f, vpl_type);
	}

	auto lc = tdCreateLightTree(left_vpls, vpl_type, min_dist, bottom_up_thresh, rng, use_pa_div);
	auto rc = tdCreateLightTree(right_vpls, vpl_type, min_dist, bottom_up_thresh, rng, use_pa_div);

	std::uniform_real_distribution<float> gen(0, 1);

	std::unique_ptr<LightTreeNode> curr_node(new LightTreeNode());
	DirConef union_cone;

	if(vpl_type != EDirectionalEmitterVPL){
		curr_node->min_bounds = Point(std::min(lc->min_bounds[0], rc->min_bounds[0]), 
			std::min(lc->min_bounds[1], rc->min_bounds[1]), 
			std::min(lc->min_bounds[2], rc->min_bounds[2]));
		curr_node->max_bounds = Point(std::max(lc->max_bounds[0], rc->max_bounds[0]), 
			std::max(lc->max_bounds[1], rc->max_bounds[1]), 
			std::max(lc->max_bounds[2], rc->max_bounds[2]));
	}

	if(vpl_type != EPointEmitterVPL){
		DirConef lc_bcone = lc->bcone;
		DirConef rc_bcone = rc->bcone;
		union_cone = DirConef::Union(lc_bcone, rc_bcone);
	}

	curr_node->bcone = union_cone;
	
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
	curr_node->num_children = 2 + lc->num_children + rc->num_children;

	curr_node->left = std::move(lc);
	curr_node->right = std::move(rc);

	float d = (curr_node->left->vpl.its.p - curr_node->right->vpl.its.p).length();

	curr_node->vpl.radius = std::min(
		curr_node->left->vpl.radius + curr_node->right->vpl.radius,
		d / 2 + std::max(curr_node->left->vpl.radius, curr_node->right->vpl.radius));

	curr_node->vpl.radius *= curr_node->vpl.radius;

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

float cosBoundB2P(const Point& max, const Point& min, const DirConef& cone, const Point& p){
	float cosHalfAngle = cone.GetAngleCos();

	if(cosHalfAngle <= 0.f){
		return 1.f;
	}

	Vector axis = cone.GetAxis();
	float cos = dot(Vector(0.0f, 0.0f, 1.0f), axis);

	Transform mt;
	if (std::abs(cos - 1.0f) < 1e-6) {
		mt = Transform::rotate(Vector(1.0f, 0.0f, 0.0f), 0.0f);
	}
	else if (std::abs(cos + 1.0f) < 1e-6) {
		mt = Transform::rotate(Vector(1.0f, 0.0f, 0.0f), -180.0f);
	}
	else {
		Vector vv = normalize(cross(Vector(0.0f, 0.0f, 1.0f), axis));
		mt = Transform::rotate(vv, radToDeg(-acosf(cos)));
	}

	Point bounding_points[8];
	bounding_points[0] = min;
	bounding_points[1] = Point(min[0], max[1], min[2]);
	bounding_points[2] = Point(max[0], max[1], min[2]);
	bounding_points[3] = Point(max[0], min[1], min[2]);
	bounding_points[4] = Point(min[0], min[1], max[2]);
	bounding_points[5] = Point(min[0], max[1], max[2]);
	bounding_points[6] = max;
	bounding_points[7] = Point(max[0], min[1], max[2]);

	AABB xbox;
	for (int i = 0; i < 8; i++) {
		xbox.expandBy(mt.transformAffine(Point(p - bounding_points[i])));
	}

	Point &xm = xbox.min;
	Point &xM = xbox.max;

	float cosTheta;
	if (xM.z > 0)
	{
		float minx2, miny2;
		if (xm.x * xM.x <= 0)
		{
			minx2 = 0.0f;
		}
		else
		{
			minx2 = std::min(xm.x * xm.x, xM.x * xM.x);
		}
			

		if (xm.y * xM.y <= 0)
		{
			miny2 = 0.0f;
		}
		else
		{
			miny2 = std::min(xm.y * xm.y, xM.y * xM.y);
		}
			
		Float maxz2 = xM.z * xM.z;
		cosTheta = xM.z / sqrt(minx2 + miny2 + maxz2);
	}
	else
	{
		cosTheta = 0.f;
	}
		

	if (cosTheta > cosHalfAngle)
	{
		return 1.0f;
	}
	else
	{
		float sinHalfAngle = sqrt(1 - cosHalfAngle * cosHalfAngle);
		float sinTheta = sqrt(1 - cosTheta * cosTheta);
		return math::clamp(cosTheta * cosHalfAngle + sinTheta * sinHalfAngle, 0.0f, 1.0f);
	}
}

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

float LightTree::calculateClusterBounds(Point shading_point_position, Normal shading_point_normal,
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
			d = std::max(min_dist, distToBox(shading_point_position,
				light_tree_node->min_bounds, light_tree_node->max_bounds));
			geometric = cosBoundB2P(light_tree_node->max_bounds, light_tree_node->min_bounds, 
				light_tree_node->bcone, shading_point_position) / (d * d);

			break;
		case EDirectionalEmitterVPL:
			geometric = 1.f;
			break;
		default:
			break;
	}

	//only dealing with lambertian diffuse right now
	float largest = -std::numeric_limits<float>::max();

	if(vpl_type != EDirectionalEmitterVPL){
		for(std::uint8_t i = 0; i < 8; ++i){
			float dp = dot(normalize(bounding_points[i] - shading_point_position), shading_point_normal);
			float angle = acos(dp);
			float coneAngle = acos(light_tree_node->bcone.GetAngleCos());
			float dif = std::max(0.f, angle - coneAngle);
			largest = std::max(float(cos(dif)), largest);
		}
		largest = std::max(0.f, largest);
	}
	else{
		float angle = acos(dot(-light_tree_node->bcone.GetAxis(), shading_point_normal));
		float coneAngle = acos(light_tree_node->bcone.GetAngleCos());
		float dif = std::max(0.f, angle - coneAngle);
		largest = std::max(0.f, float(cos(dif)));
	}

	float material = std::max(0.f, largest);

	return geometric * material * light_tree_node->emission_scale;
}

LightTree::LightTree() : point_tree_root_(nullptr), directional_tree_root_(nullptr), oriented_tree_root_(nullptr), min_dist_(0.f) {

}

LightTree::LightTree(const std::vector<VPL>& vpls, float min_dist, std::uint32_t max_lights, float error_threshold, bool divbypa) : 
	min_dist_(min_dist), max_lights_(max_lights), error_threshold_(error_threshold), divbypa_(divbypa) {
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

	point_tree_root_ = tdCreateLightTree(point_vpls_, EPointEmitterVPL, min_dist_, 1, rng, divbypa_);
	oriented_tree_root_ = tdCreateLightTree(oriented_vpls_, ESurfaceVPL, min_dist_, 1, rng, divbypa_);
	directional_tree_root_ = tdCreateLightTree(directional_vpls_, EDirectionalEmitterVPL, min_dist_, 1, rng, divbypa_);
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

	std::vector<VPL> lights;
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

	lightcut.insert(lightcut.end(), lights.begin(), lights.end());

	return lightcut;
}

MTS_NAMESPACE_END