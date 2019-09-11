#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include <mitsuba/core/plugin.h>

#include <string>
#include <vector>

MTS_NAMESPACE_BEGIN

void writeOutputImage(Scene* scene, std::string filename, std::uint32_t width, std::uint32_t height, bool hdr, const std::vector<Spectrum>& data);
float writeOutputErrorImage(Scene* scene, std::string filename, std::uint32_t width, std::uint32_t height, bool hdr, const std::vector<Spectrum>& data1,
    const std::vector<Spectrum>& data2, float max_error);
void writeOutputData(std::string filename, bool new_file, const std::vector<float>& data, char delimiter);

MTS_NAMESPACE_END

#endif