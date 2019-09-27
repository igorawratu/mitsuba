struct PixelElement{
    cl_float dr, dg, db;
    cl_float sr, sg, sb;
    cl_float x, y, z;
    cl_float nx, ny, nz;
    cl_float roughness;
    cl_float eta_r, eta_g, eta_b, k_r, k_g, k_b;
    cl_float wi_x, wi_y, wi_z;
};

struct LightElement{
    cl_float r, g, b;
    cl_float nx, ny, nz;
    cl_float x, y, z;
    cl_float rad;
};

struct OutputElement{
    cl_float r, g, b;
};

__kernel void shade(__global const PixelElement* A, __global const LightElement* B, __global OutputElement *C){
    int i = get_global_id(0);
    
    
}