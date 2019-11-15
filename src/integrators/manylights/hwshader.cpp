#include "hwshader.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <chrono>

#include "hwshader_ocl.h"

MTS_NAMESPACE_BEGIN

HWShader::HWShader() : 
    context_(nullptr),
    command_queue_(nullptr),
    platform_id_(nullptr),
    device_id_(nullptr),
    program_(nullptr),
    kernel_(nullptr),
    vsl_kernel_(nullptr),
    initialized_(false),    
    pixel_buffer_(nullptr),
    light_buffer_(nullptr),
    output_buffer_(nullptr),
    curr_buffer_elements_(0){
    
    cl_uint num_devices;
    cl_uint num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id_, &num_platforms);
    ret = clGetDeviceIDs(platform_id_, CL_DEVICE_TYPE_DEFAULT, 1, &device_id_, &num_devices);

    if(ret != 0){
        std::cerr << "Unable to find opencl device" << std::endl;
        return;
    }

    context_ = clCreateContext( NULL, 1, &device_id_, NULL, NULL, &ret);
    if(ret != 0){
        std::cerr << "Unable to create opencl context" << std::endl;
        return;
    }

    command_queue_ = clCreateCommandQueue(context_, device_id_, 0, &ret);
    if(ret != 0){
        std::cerr << "Error creating command queue" << std::endl;
        return;
    }

    size_t source_code_size = strlen(hwshader_ocl);

    program_ = clCreateProgramWithSource(context_, 1, 
            (const char **)&hwshader_ocl, (const size_t *)&source_code_size, &ret);
    if(ret != 0){
        std::cerr << "Unable to create CL program" << std::endl;
        return;
    }

    ret = clBuildProgram(program_, 1, &device_id_, NULL, NULL, NULL);
    if(ret != 0){
        std::cerr << "Unable to build CL program" << std::endl;

        if(ret == CL_BUILD_PROGRAM_FAILURE){
            size_t log_size;
            clGetProgramBuildInfo(program_, device_id_, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
            std::unique_ptr<char[]> log(new char[log_size]);
            clGetProgramBuildInfo(program_, device_id_, CL_PROGRAM_BUILD_LOG, log_size, log.get(), NULL);
            printf("%s\n", log.get());
        }

        return;
    }

    kernel_ = clCreateKernel(program_, "shade", &ret);
    if(ret != 0){
        std::cerr << "Unable to create shading kernel" << std::endl;
        return;
    }

    vsl_kernel_ = clCreateKernel(program_, "shadeVSL", &ret);
    if(ret != 0){
        std::cerr << "Unable to create vsl kernel" << std::endl;
        return;
    }

    std::cout << "OpenCL initialized" << std::endl;

    initialized_ = true;
}

HWShader::~HWShader(){
    if(command_queue_){
        clFlush(command_queue_);
        clFinish(command_queue_);
    }
        
    if(kernel_){
        clReleaseKernel(kernel_);
    }
    
    if(vsl_kernel_){
        clReleaseKernel(vsl_kernel_);
    }

    if(program_){
        clReleaseProgram(program_);
    }
    
    if(pixel_buffer_){
        clReleaseMemObject(pixel_buffer_);
    }
    
    if(light_buffer_){
        clReleaseMemObject(light_buffer_);
    }

    if(output_buffer_){
        clReleaseMemObject(output_buffer_);
    }
    
    if(command_queue_){
        clReleaseCommandQueue(command_queue_);
    }
    
    if(context_){
        clReleaseContext(context_);
    }
}

struct PixelElement{
    cl_float3 diff_col;
    cl_float3 spec_col;
    cl_float3 p;
    cl_float3 n;
    cl_float roughness;
    cl_float3 eta;
    cl_float3 k;
    cl_float3 wi;
    cl_int type;
};

struct LightElement{
    cl_float3 power;
    cl_float3 diff_col;
    cl_float3 spec_col;
    cl_float3 p;
    cl_float3 n;
    cl_float roughness;
    cl_float3 eta;
    cl_float3 k;
    cl_float3 wi;
    cl_float coeff;
    cl_float rad;
    cl_int type;
    cl_int light_surface_type;
};

struct OutputElement{
    cl_float3 col;
};

bool HWShader::initializeBuffers(std::uint32_t size){
    cl_int ret;

    if(pixel_buffer_){
        ret = clReleaseMemObject(pixel_buffer_);
    }
    
    if(output_buffer_){
        ret = clReleaseMemObject(output_buffer_);
    }

    if(light_buffer_){
        ret = clReleaseMemObject(light_buffer_);
    }

    pixel_buffer_ = clCreateBuffer(context_, CL_MEM_READ_ONLY, 
        size * sizeof(PixelElement), NULL, &ret);
    if(ret != 0){
        std::cerr << "Unable to create pixel opencl buffer" << std::endl;
        return false;
    }

    light_buffer_ = clCreateBuffer(context_, CL_MEM_READ_ONLY, 
        size * sizeof(LightElement), NULL, &ret);
    if(ret != 0){
        std::cerr << "Unable to create light opencl buffer" << std::endl;
        return false;
    }

    output_buffer_ = clCreateBuffer(context_, CL_MEM_READ_WRITE, 
        size * sizeof(OutputElement), NULL, &ret);
    if(ret != 0){
        std::cerr << "Unable to create output opencl buffer" << std::endl;
        return false;
    }

    ret = clSetKernelArg(kernel_, 0, sizeof(cl_mem), (void *)&pixel_buffer_);
    if(ret != 0){
        std::cerr << "Unable to copy to set kernel argument 0" << std::endl;
        return false;
    }

    ret = clSetKernelArg(kernel_, 1, sizeof(cl_mem), (void *)&light_buffer_);
    if(ret != 0){
        std::cerr << "Unable to copy to set kernel argument 1" << std::endl;
        return false;
    }

    ret = clSetKernelArg(kernel_, 2, sizeof(cl_mem), (void *)&output_buffer_);
    if(ret != 0){
        std::cerr << "Unable to copy to set kernel argument 2" << std::endl;
        return false;
    }

    ret = clSetKernelArg(vsl_kernel_, 0, sizeof(cl_mem), (void *)&pixel_buffer_);
    if(ret != 0){
        std::cerr << "Unable to copy to set kernel argument 0" << std::endl;
        return false;
    }

    ret = clSetKernelArg(vsl_kernel_, 1, sizeof(cl_mem), (void *)&light_buffer_);
    if(ret != 0){
        std::cerr << "Unable to copy to set kernel argument 1" << std::endl;
        return false;
    }

    ret = clSetKernelArg(vsl_kernel_, 2, sizeof(cl_mem), (void *)&output_buffer_);
    if(ret != 0){
        std::cerr << "Unable to copy to set kernel argument 2" << std::endl;
        return false;
    }

    return true;
}

void HWShader::renderSlices(const std::vector<KDTNode<ReconstructionSample>*>& slices,
    const std::vector<std::vector<VPL>*>& vpls, std::uint32_t cluster_size, float min_dist, bool vsl){

    assert(vpls.size() == slices.size());

    if(!initialized_){
        return;
    }

    std::uint32_t num_elements = 0;
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        num_elements += slices[i]->sample_indices.size();
    }

    if(num_elements > curr_buffer_elements_){
        if(!initializeBuffers(num_elements)){
            initialized_ = false;
            return;
        }

        curr_buffer_elements_ = num_elements;
    }

    cl_int num_el_cl = num_elements;

    cl_int ret;

    if(vsl){
        ret = clSetKernelArg(vsl_kernel_, 3, sizeof(cl_int), (void *)&num_el_cl);
        if(ret != 0){
            initialized_ = false;
            std::cerr << "Unable to copy to set vsl kernel argument 3" << std::endl;
            return;
        }

        cl_float cl_min_dist = min_dist;
        ret = clSetKernelArg(vsl_kernel_, 4, sizeof(cl_float), (void *)&cl_min_dist);
        if(ret != 0){
            initialized_ = false;
            std::cerr << "Unable to set vsl kernel argument 4" << std::endl;
            return;
        }

        cl_int max_samples = 100;
        ret = clSetKernelArg(vsl_kernel_, 5, sizeof(cl_int), (void *)&max_samples);
        if(ret != 0){
            initialized_ = false;
            std::cerr << "Unable to set vsl kernel argument 5" << std::endl;
            return;
        }
    }
    else{
        ret = clSetKernelArg(kernel_, 3, sizeof(cl_int), (void *)&num_el_cl);
        if(ret != 0){
            initialized_ = false;
            std::cerr << "Unable to copy to set kernel argument 3" << std::endl;
            return;
        }

        cl_float cl_min_dist = min_dist;
        ret = clSetKernelArg(kernel_, 4, sizeof(cl_float), (void *)&cl_min_dist);
        if(ret != 0){
            initialized_ = false;
            std::cerr << "Unable to set kernel argument 4" << std::endl;
            return;
        }
    }
    

    std::unique_ptr<PixelElement[]> host_pixel_buffer(new PixelElement[num_elements]);
    std::unique_ptr<LightElement[]> host_light_buffer(new LightElement[num_elements]);
    std::unique_ptr<OutputElement[]> host_output_buffer(new OutputElement[num_elements]);
    for(std::uint32_t i = 0; i < num_elements; ++i){
        host_output_buffer[i].col.s[0] = 0.f;
        host_output_buffer[i].col.s[1] = 0.f;
        host_output_buffer[i].col.s[2] = 0.f;
    }

    ret = clEnqueueWriteBuffer(command_queue_, output_buffer_, CL_TRUE, 0,
        num_elements * sizeof(OutputElement), host_output_buffer.get(), 0, NULL, NULL);
    if(ret != 0){
        std::cerr << "Unable to copy to output opencl buffer" << std::endl;
        initialized_ = false;
        return;
    }

    std::uint32_t elements_processed = 0;
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
            std::uint32_t curr_element = elements_processed + j;

            ReconstructionSample& sample = slices[i]->sample(j);
            const BSDF* bsdf = sample.its.getBSDF();

            Spectrum diffuse_col = bsdf->getDiffuseReflectance(sample.its);
            diffuse_col.toLinearRGB(host_pixel_buffer[curr_element].diff_col.s[0], 
                host_pixel_buffer[curr_element].diff_col.s[1], host_pixel_buffer[curr_element].diff_col.s[2]);
            
            Spectrum spec_col = bsdf->getSpecularReflectance(sample.its);
            spec_col.toLinearRGB(host_pixel_buffer[curr_element].spec_col.s[0], 
                host_pixel_buffer[curr_element].spec_col.s[1], host_pixel_buffer[curr_element].spec_col.s[2]);

            host_pixel_buffer[curr_element].p.s[0] = sample.its.p.x;
            host_pixel_buffer[curr_element].p.s[1] = sample.its.p.y;
            host_pixel_buffer[curr_element].p.s[2] = sample.its.p.z;
            Normal n = sample.its.wi.z > 0.f ? sample.its.shFrame.n : -sample.its.shFrame.n;
            host_pixel_buffer[curr_element].n.s[0] = n.x;
            host_pixel_buffer[curr_element].n.s[1] = n.y;
            host_pixel_buffer[curr_element].n.s[2] = n.z;

            host_pixel_buffer[curr_element].roughness = std::min(1.5f, bsdf->getRoughness(sample.its, 0));
            
            Spectrum eta = bsdf->getEtaSpec(sample.its.wi);
            eta.toLinearRGB(host_pixel_buffer[curr_element].eta.s[0], 
                host_pixel_buffer[curr_element].eta.s[1], host_pixel_buffer[curr_element].eta.s[2]);

            Spectrum k = bsdf->getK(sample.its.wi);
            k.toLinearRGB(host_pixel_buffer[curr_element].k.s[0], host_pixel_buffer[curr_element].k.s[1],
                host_pixel_buffer[curr_element].k.s[2]);

            Vector wi = sample.its.toWorld(sample.its.wi);

            host_pixel_buffer[curr_element].wi.s[0] = wi.x;
            host_pixel_buffer[curr_element].wi.s[1] = wi.y;
            host_pixel_buffer[curr_element].wi.s[2] = wi.z;

            host_pixel_buffer[curr_element].type = bsdf->isConductor() ? 0 : 1;
        }
        elements_processed += slices[i]->sample_indices.size();
    }

    ret = clEnqueueWriteBuffer(command_queue_, pixel_buffer_, CL_TRUE, 0,
        num_elements * sizeof(PixelElement), host_pixel_buffer.get(), 0, NULL, NULL);
    if(ret != 0){
        std::cerr << "Unable to copy to pixel opencl buffer" << std::endl;
        initialized_ = false;
        return;
    }

    for(std::uint32_t i = 0; i < cluster_size; ++i){
        elements_processed = 0;
        for(std::uint32_t j = 0; j < slices.size(); ++j){
            VPL& vpl = (*vpls[j])[i];

            LightElement light_for_slice;

            vpl.P.toLinearRGB(light_for_slice.power.s[0], light_for_slice.power.s[1],
                    light_for_slice.power.s[2]);
            light_for_slice.n.s[0] = vpl.its.shFrame.n.x;
            light_for_slice.n.s[1] = vpl.its.shFrame.n.y;
            light_for_slice.n.s[2] = vpl.its.shFrame.n.z;

            light_for_slice.p.s[0] = vpl.its.p.x;
            light_for_slice.p.s[1] = vpl.its.p.y;
            light_for_slice.p.s[2] = vpl.its.p.z;

            light_for_slice.rad = vpl.radius;

            if(vpl.type == EDirectionalEmitterVPL){
                light_for_slice.type = 0;
            }
            else if(vpl.type == EPointEmitterVPL){
                light_for_slice.type = 1;
            }
            else{
                light_for_slice.type = 2;
                Vector light_wi = vpl.its.toWorld(vpl.its.wi);

                light_for_slice.wi.s[0] = light_wi.x;
                light_for_slice.wi.s[1] = light_wi.y;
                light_for_slice.wi.s[2] = light_wi.z;
                
                if(vpl.emitter != nullptr || vpl.its.getBSDF() == nullptr){
                    light_for_slice.diff_col.s[0] = 1.f;
                    light_for_slice.diff_col.s[1] = 1.f;
                    light_for_slice.diff_col.s[2] = 1.f;

                    light_for_slice.spec_col.s[0] = 0.f;
                    light_for_slice.spec_col.s[1] = 0.f;
                    light_for_slice.spec_col.s[2] = 0.f;

                    light_for_slice.eta.s[0] = 1.f;
                    light_for_slice.eta.s[0] = 1.f;
                    light_for_slice.eta.s[0] = 1.f;

                    light_for_slice.k.s[0] = 0.f;
                    light_for_slice.k.s[0] = 0.f;
                    light_for_slice.k.s[0] = 0.f;

                    light_for_slice.roughness = 1.5f;

                    light_for_slice.light_surface_type = 1;
                }
                else{
                    const BSDF* bsdf = vpl.its.getBSDF();
                    Spectrum ldiffuse_col = bsdf->getDiffuseReflectance(vpl.its);
                    Spectrum lspecular_col = bsdf->getSpecularReflectance(vpl.its);
                    Spectrum eta = bsdf->getEtaSpec(vpl.its.wi);
                    Spectrum k = bsdf->getK(vpl.its.wi);

                    ldiffuse_col.toLinearRGB(light_for_slice.diff_col.s[0], light_for_slice.diff_col.s[1],
                        light_for_slice.diff_col.s[2]);
                    lspecular_col.toLinearRGB(light_for_slice.spec_col.s[0], light_for_slice.spec_col.s[1],
                        light_for_slice.spec_col.s[2]);
                    eta.toLinearRGB(light_for_slice.eta.s[0], light_for_slice.eta.s[1],
                        light_for_slice.eta.s[2]);
                    k.toLinearRGB(light_for_slice.k.s[0], light_for_slice.k.s[1],
                        light_for_slice.k.s[2]);

                    light_for_slice.roughness = std::min(1.5f, bsdf->getRoughness(vpl.its, 0));

                    light_for_slice.light_surface_type = bsdf->isConductor() ? 0 : 1;
                }
            }

            for(std::uint32_t sample_idx = 0; sample_idx < slices[j]->sample_indices.size(); ++sample_idx){
                std::uint32_t curr_element = elements_processed + sample_idx;
                host_light_buffer[curr_element] = light_for_slice;

                std::uint32_t idx = i * slices[j]->sample_indices.size() + sample_idx;
                host_light_buffer[curr_element].coeff = slices[j]->visibility_coefficients[idx];
            }

            elements_processed += slices[j]->sample_indices.size();
        }

        ret = clEnqueueWriteBuffer(command_queue_, light_buffer_, CL_TRUE, 0,
            num_elements * sizeof(LightElement), host_light_buffer.get(), 0, NULL, NULL);
        if(ret != 0){
            std::cerr << "Unable to copy to light opencl buffer" << std::endl;
            initialized_ = false;
            return;
        }
        
        size_t local_item_size = 64;
        size_t global_item_size = (num_elements / local_item_size + 1) * local_item_size;
        
        if(vsl){
            ret = clEnqueueNDRangeKernel(command_queue_, vsl_kernel_, 1, NULL, 
                &global_item_size, &local_item_size, 0, NULL, NULL);
        }
        else{
            ret = clEnqueueNDRangeKernel(command_queue_, kernel_, 1, NULL, 
                &global_item_size, &local_item_size, 0, NULL, NULL);
        }
        
        if(ret != 0){
            std::cerr << "Unable to run opencl kernel with code " << ret << std::endl;
            initialized_ = false;
            return;
        }
    }

    ret = clEnqueueReadBuffer(command_queue_, output_buffer_, CL_TRUE, 0, 
            num_elements * sizeof(OutputElement), host_output_buffer.get(), 0, NULL, NULL);
    if(ret != 0){
        std::cerr << "Unable to read from output opencl buffer" << std::endl;
        initialized_ = false;
        return;
    }

    std::uint32_t curr_element = 0;
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        slices[i]->visibility_coefficients.clear();
        slices[i]->visibility_coefficients.shrink_to_fit();
        
        for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
            auto& sample = slices[i]->sample(j);

            sample.color = sample.its.isEmitter() ? sample.its.Le(-sample.ray.d) : 
                Spectrum(host_output_buffer[curr_element].col.s);

            curr_element++;
        }
    }
}

MTS_NAMESPACE_END