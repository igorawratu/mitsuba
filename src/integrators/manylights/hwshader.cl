struct PixelElement{
    float dr, dg, db;
    float sr, sg, sb;
    float x, y, z;
    float nx, ny, nz;
    float roughness;
    float eta_r, eta_g, eta_b, k_r, k_g, k_b;
    float wi_x, wi_y, wi_z;
};

struct LightElement{
    float r, g, b;
    float nx, ny, nz;
    float x, y, z;
    float rad;
};

struct OutputElement{
    float r, g, b;
};

__kernel void shade(__global const struct PixelElement* A, __global const struct LightElement* B, 
    __global struct OutputElement *C, int num_pixels){
    int i = get_global_id(0);
    if(i >= num_pixels){
        return;
    }
    
    C[i].r = A[i].dr / 1000;
    C[i].g = A[i].dg / 1000;
    C[i].b = A[i].db / 1000;
}