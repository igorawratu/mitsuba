#ifndef COMMON_H
#define COMMON_H

#include <mitsuba/core/plugin.h>
#include <mitsuba/render/vpl.h>
#include <cstdint>

#include "definitions.h"

MTS_NAMESPACE_BEGIN

float calculateError(Scene* scene, const std::vector<VPL>& vpls, float min_dist, std::uint8_t* image_buffer);

MTS_NAMESPACE_END

#endif