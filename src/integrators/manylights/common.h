#ifndef COMMON_H
#define COMMON_H

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/vpl.h>
#include <cstdint>
#include <tuple>
#include <eigen3/Eigen/Dense>

#include "definitions.h"

MTS_NAMESPACE_BEGIN

float calculateError(Scene* scene, const std::vector<VPL>& vpls, float min_dist, std::uint8_t* image_buffer);

std::tuple<float, float, float> floatToRGB(float v);

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf> partialSvd(const Eigen::MatrixXf& mat, std::uint32_t num_singular_values);

Eigen::MatrixXf softThreshRank(const Eigen::MatrixXf& mat, float theta, const std::uint32_t initial, const std::uint32_t step_size);

Spectrum sample(Scene* scene, Sampler* sampler, const Intersection& its, const VPL& vpl, float min_dist, bool check_occlusion);

Eigen::MatrixXf computeMoorePenroseInverse(const Eigen::MatrixXf& m);

MTS_NAMESPACE_END

#endif