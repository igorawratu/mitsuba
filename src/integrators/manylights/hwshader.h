#ifndef HWSHADER_H_
#define HWSHADER_H_

#include "common.h"
#include "matrixreconstruction.h"
#include <vector>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>

MTS_NAMESPACE_BEGIN

class HWShader{
public:
    HWShader();
    HWShader(const HWShader& other) = delete;
    ~HWShader();

    void renderSlices(const std::vector<KDTNode<ReconstructionSample>*>& slices,
        const std::vector<std::vector<VPL>*>& vpls, std::uint32_t cluster_size,
        float min_dist, bool vsl);

    void renderHWBF(std::vector<HWBFPix>& receivers, const std::vector<VPL>& vpls, std::uint32_t start,
        std::uint32_t end, float min_dist, bool vsl, Scene* scene);
    
    bool initialized(){
        return initialized_;
    }

private:
    bool initializePixelBuffers(std::uint32_t pixels);
    bool initializeLightBuffer(std::uint32_t slices, std::uint32_t clusters_per_slice);

private:
    cl_context context_;
    cl_command_queue command_queue_;
    cl_platform_id platform_id_;
    cl_device_id device_id_;
    cl_program program_;
    cl_kernel kernel_, vsl_kernel_;
    bool initialized_;
    
    cl_mem pixel_buffer_, light_buffer_, output_buffer_, coeff_buffer_;
    std::uint32_t curr_buffer_elements_, curr_light_elements_;
};

MTS_NAMESPACE_END

#endif