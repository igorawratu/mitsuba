#include "hwshader.h"

#include <iostream>
#include <memory>
#include <fstream>

MTS_NAMESPACE_BEGIN

HWShader::HWShader() : 
    context_(nullptr),
    command_queue_(nullptr),
    platform_id_(nullptr),
    device_id_(nullptr),
    program_(nullptr),
    kernel_(nullptr),
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

    std::ifstream source_file("/disc/swang/mitsuba/src/integrators/manylights/hwshader.cl");
    std::string source_code;
    source_code.assign(std::istreambuf_iterator<char>(source_file), std::istreambuf_iterator<char>());

    const char* source_code_cstr = source_code.c_str();
    size_t source_code_size = source_code.length();

    program_ = clCreateProgramWithSource(context_, 1, 
            (const char **)&source_code_cstr, (const size_t *)&source_code_size, &ret);
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

    kernel_ = clCreateKernel(program_, "vector_add", &ret);
    if(ret != 0){
        std::cerr << "Unable to create shading kernel" << std::endl;
        return;
    }

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

};

struct LightElement{

};

struct OutputElement{
    cl_float r;
    cl_float g;
    cl_float b;
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

    output_buffer_ = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, 
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

    return true;
}

void HWShader::renderSlices(const std::vector<KDTNode<ReconstructionSample>*>& slices,
    const std::vector<std::vector<VPL>*>& vpls, std::uint32_t cluster_size){

    assert(vpls.size() == slices.size());

    if(!initialized_){
        return;
    }

    for(std::uint32_t i = 0; i < slices.size(); ++i){
        for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
            slices[i]->sample(j).unoccluded_samples.resize(cluster_size);
        }
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

    std::unique_ptr<PixelElement[]> host_pixel_buffer(new PixelElement[num_elements]);
    std::unique_ptr<LightElement[]> host_light_buffer(new LightElement[num_elements]);
    std::unique_ptr<OutputElement[]> host_output_buffer(new OutputElement[num_elements]);

    std::uint32_t curr_element = 0;
    for(std::uint32_t i = 0; i < slices.size(); ++i){
        for(std::uint32_t j = 0; j < slices[i]->sample_indices.size(); ++j){
            //assign input elements here

            curr_element++;
        }
    }

    cl_int ret = clEnqueueWriteBuffer(command_queue_, pixel_buffer_, CL_TRUE, 0,
        num_elements * sizeof(PixelElement), host_pixel_buffer.get(), 0, NULL, NULL);
    if(ret != 0){
        std::cerr << "Unable to copy to pixel opencl buffer" << std::endl;
        initialized_ = false;
        return;
    }

    for(std::uint32_t i = 0; i < cluster_size; ++i){
        curr_element = 0;
        for(std::uint32_t j = 0; j < slices.size(); ++j){
            for(std::uint32_t k = 0; k < slices[i]->sample_indices.size(); ++k){
                //assign ith light props here
                
                curr_element++;
            }
        }

        ret = clEnqueueWriteBuffer(command_queue_, light_buffer_, CL_TRUE, 0,
            num_elements * sizeof(LightElement), host_light_buffer.get(), 0, NULL, NULL);
        if(ret != 0){
            std::cerr << "Unable to copy to light opencl buffer" << std::endl;
            initialized_ = false;
            return;
        }

        size_t global_item_size = num_elements;
        size_t local_item_size = 64;
        ret = clEnqueueNDRangeKernel(command_queue_, kernel_, 1, NULL, 
                &global_item_size, &local_item_size, 0, NULL, NULL);
        if(ret != 0){
            std::cerr << "Unable to run opencl kernel" << std::endl;
            initialized_ = false;
            return;
        }

        ret = clEnqueueReadBuffer(command_queue_, output_buffer_, CL_TRUE, 0, 
                num_elements * sizeof(OutputElement), host_output_buffer.get(), 0, NULL, NULL);
        if(ret != 0){
            std::cerr << "Unable to read from output opencl buffer" << std::endl;
            initialized_ = false;
            return;
        }

        curr_element = 0;
        for(std::uint32_t j = 0; j < slices.size(); ++j){
            for(std::uint32_t k = 0; k < slices[i]->sample_indices.size(); ++k){
                OutputElement& curr_out = host_output_buffer[curr_element++];
                slices[j]->sample(k).unoccluded_samples[i].fromLinearRGB(curr_out.r, curr_out.g, curr_out.b);
            }
        }
    }
}

MTS_NAMESPACE_END