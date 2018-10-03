#include "lighttree.h"

LightTree::LightTree() : point_tree_root_(nullptr), directional_tree_root_(nullptr), oriented_tree_root_(nullptr){

}

void calculateSimilarityMatrices(const std::vector<VPL>& point_vpls, const std::vector<VPL>& directional_vpls, const std::vector<VPL>& oriented_vpls,
	std::vector<float>& point_similarity_matrix, std::vector<float>& direction_similarity_matrix, std::vector<float>& oriented_similarity_matrix) {
	
	point_similarity_matrix.clear();
	if (point_vpls.size() > 0) {
		point_similarity_matrix.insert(point_similarity_matrix.begin(), point_vpls.size() * point_vpls.size(), 0.f);
	}

	direction_similarity_matrix.clear();
	if (directional_vpls.size() > 0) {
		direction_similarity_matrix.insert(direction_similarity_matrix.begin(), directional_vpls.size() * directional_vpls.size(), 0.f);
	}

	oriented_similarity_matrix.clear();
	if (oriented_vpls.size() > 0) {
		oriented_similarity_matrix.insert(oriented_similarity_matrix.begin(), oriented_vpls.size() * oriented_vpls.size(), 0.f);
	}

	for (std::uint32_t i = 0; i < point_vpls.size(); ++i) {
		for (std::uint32_t j = 0; j < point_vpls.size(); ++j) {
			std::uint32_t index = j + i * point_vpls.size();
			point_similarity_matrix[index] = (point_vpls[i].its.p - point_vpls[j].its.p).length();
		}
	}

	for (std::uint32_t i = 0; i < directional_vpls.size(); ++i) {
		for (std::uint32_t j = 0; j < directional_vpls.size(); ++j) {
			std::uint32_t index = j + i * directional_vpls.size();
			direction_similarity_matrix[index] = dot(directional_vpls[i].its.shFrame.n, directional_vpls[j].its.shFrame.n);
		}
	}

	for (std::uint32_t i = 0; i < oriented_vpls.size(); ++i) {
		for (std::uint32_t j = 0; j < oriented_vpls.size(); ++j) {
			std::uint32_t index = j + i * oriented_vpls.size();
			oriented_similarity_matrix[index] = (oriented_vpls[i].its.p - oriented_vpls[j].its.p).length() + dot(oriented_vpls[i].its.shFrame.n, oriented_vpls[j].its.shFrame.n);
		}
	}
}

std::unique_ptr<LightTreeNode> createLightTree(const std::vector<VPL>& vpls, const std::vector<float>& similarity_matrix) {

}

LightTree::LightTree(const std::vector<VPL>& vpls) {
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

	calculateSimilarityMatrices(point_vpls_, directional_vpls_, oriented_vpls_, point_similarity_matrix_, directional_similarity_matrix_, oriented_similarity_matrix_);
}

LightTree::LightTree(const LightTree& other) {

}

LightTree::LightTree(const LightTree&& other) {

}

LightTree& LightTree::operator = (const LightTree& other) {

}

LightTree& LightTree::operator = (LightTree&& other) {

}

LightTree::~LightTree() {

}

void LightTree::setVPLs(const std::vector<VPL>& vpls) {

}

std::vector<VPL> LightTree::getClusteringForPoint(Intersection its) {

}