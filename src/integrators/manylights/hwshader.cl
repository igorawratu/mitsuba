struct PixelElement{
    float3 diff_col;
    float3 spec_col;
    float3 p;
    float3 n;
    float roughness;
    float3 eta;
    float3 k;
    float3 wo;
};

struct LightElement{
    float3 power;
    float3 diff_col;
    float3 spec_col;
    float3 p;
    float3 n;
    float roughness;
    float3 eta;
    float3 k;
    float3 wi;
    float coeff;
    float rad;
    int type;
};

struct OutputElement{
    float3 col;
};

#define PI 3.1415926

float3 lambertian(float3 albedo, float3 p, float3 n, float3 lp){
    float3 dir = lp - p;
    float3 wi = normalize(dir);
    float wi_dot_n = fmax(0.0f, dot(wi, n));
    
    float3 col = albedo * (float3)(wi_dot_n) / (float3)(PI);

    return col;
}

float d_ggx(float3 h, float3 n, float alpha){
    float alpha2 = alpha * alpha;
    float noh = clamp(dot(n, h), 0.0f, 1.0f);
    float noh2 = noh * noh;
    float denom = noh2 * alpha2 + (1.0f - noh2);

    return alpha2 / max(0.00001f, (float)(PI * denom * denom));
}

float g_ggx_smith(float alpha, float3 n, float3 v, float3 l){
    float alpha2 = alpha * alpha;
    float nov = clamp(dot(n, v), 0.00001f, 1.0f);
    float nol = clamp(dot(n, l), 0.00001f, 1.0f);
    float tanv = (1.0f - nov) / nov;
    float tanl = (1.0f - nol) / nol;

    float gv = 2.0f / (1.0f + sqrt(1.0f + alpha2 * tanv * tanv));
    float gl = 2.0f / (1.0f + sqrt(1.0f + alpha2 * tanl * tanl));

    return gv * gl;
}

float3 fresnel(float3 v, float3 h, float3 eta, float3 k) {
    float costheta = clamp(dot(v, h), 0.0f, 1.0f);
    float cos2theta = costheta * costheta;
    float sin2theta = 1.0f - cos2theta;
    float sin4theta = sin2theta * sin2theta;

    float3 temp1 =  eta * eta - k * k - (float3)(sin2theta);
    float3 a2pb2 = sqrt(temp1 * temp1 + k * k * eta * eta * 4.0f);
    float3 a = sqrt((a2pb2 + temp1) * 0.5f);

    float3 term1 = a2pb2 + (float3)(cos2theta);
    float3 term2 = a * (2.0f * costheta);

    float3 rs2 = (term1 - term2) / (term1 + term2);

    float3 term3 = a2pb2 * cos2theta + (float3)(sin4theta);
    float3 term4 = term2 * sin2theta;

    float3 rp2 = rs2 * (term3 - term4) / (term3 + term4);

    return 0.5f * (rp2 + rs2);
}

float3 cookTorrance(float3 pixel_p, float3 light_p, float3 pixel_wo, float3 pixel_n, float3 pixel_rough, float3 pixel_eta, float3 pixel_k,
    float3 pixel_spec_col){

    float3 l = normalize(light_p - pixel_p);
    float3 h = normalize(pixel_wo + l);
    float3 n = normalize(pixel_n);
    
    float d = d_ggx(h, n, pixel_rough);
    float g = g_ggx_smith(pixel_rough, n, pixel_wo, l);
    float3 f = fresnel(pixel_wo, h, pixel_eta, pixel_k) * pixel_spec_col;

    float nol = clamp(dot(n, l), 0.f, 1.f);
    float nov = clamp(dot(pixel_wo, n), 0.f, 1.f);

    return f * d * g / clamp(4.0f * (nol * nov + 0.05f), 0.f, 1.f);
}

float3 getLightCol(int light_type, float3 light_power, float3 pixel_p, float3 light_p, float3 light_n, float3 light_wi, float3 light_rough,
    float3 light_eta, float3 light_k, float3 light_spec_col){
    if(light_type != 2){
        return light_power;
    }

    float3 dir = normalize(pixel_p - light_p);
    
    float nodir = clamp(dot(light_n, dir), 0.0f, 1.0f);
    float3 diffuse_col = light_power * nodir / (float)(PI);

    float3 h = normalize(dir + light_wi);
    
    float d = d_ggx(h, light_n, light_rough);
    float g = g_ggx_smith(light_rough, light_n, dir, light_wi);
    float3 f = fresnel(dir, h, light_eta, light_k) * light_power * light_spec_col;

    float nol = clamp(dot(light_n, light_wi), 0.f, 1.f);
    float nov = clamp(dot(dir, light_n), 0.f, 1.f);

    float3 specular_col = light_rough < 1.000001f ? f * d * g / clamp(4.0f * (nol * nov + 0.05f), 0.f, 1.f) : (float3)(0.0f);

    return diffuse_col + specular_col;
}

__kernel void shade(__global const struct PixelElement* pixels, __global const struct LightElement* lights, 
    __global struct OutputElement *output, int num_pixels, float min_dist){
    int i = get_global_id(0);
    if(i >= num_pixels){
        return;
    }

    float coeff = (lights[i].coeff + 1.0f) / 2.0f;

    float3 diffuse = lambertian(pixels[i].diff_col, pixels[i].p, pixels[i].n, lights[i].p);

    float3 specular = pixels[i].roughness < 1.0001f ? 
        cookTorrance(pixels[i].p, lights[i].p, pixels[i].wo, pixels[i].n, pixels[i].roughness, pixels[i].eta, pixels[i].k, pixels[i].spec_col) : 
        (float3)(0.0f);

    float3 bsdf_color = diffuse + specular;

    float3 light_col = getLightCol(lights[i].type, lights[i].power, pixels[i].p, light[i].p, lights[i].n, lights[i].wi, light[i].roughness,
        lights[i].eta, lights[i].k, lights[i].spec_col);

    if(lights[i].type != 0){
        float3 dir = pixels[i].p - lights[i].p;
        float dist = fmax(min_dist, length(dir));
        light_col /= (dist * dist);
    }

    output[i].col = output[i].col + bsdf_color * light_col * (float3)(coeff);
}

float rand(unsigned int seed){
    unsigned int c = seed>>32;
    unsigned int x = seed & 0xFFFFFFFF;
    seed = (unsigned long)(x * (unsigned long)(4294883355) + c);
    return (float)(x ^ c) / ;
}

float3 rotate_to_frame(float3 old_frame, float3 new_frame, float3 v){
    if(old_frame == -new_frame){
        old_frame += (0.001f, 0.001f, 0.001f);
    }

    float4 cp = cross((float4)(old_frame, 0.0f), (float4)(new_frame, 0.0f));
    float sin = cp.length();
    float cos = dot(old_frame, new_frame);

    float v2s = 1.0f / (1.0f + c);
    float3[3] mat = {(float3)(1.0f, -cp.z + cp.z * cp.z * v2s, cp.y + cp.y * cp.y * v2s),
                    (float3)(cp.z + cp.z * cp.z * v2s, 1.0f, -cp.x + cp.x * cp.x * v2s),
                    (float3)(-cp.y + cp.y * cp.y * vs2, cp.x + cp.x * cp.x * v2s, 1.0f)};

    //assuming v on the right
    return (float3)(mat[0].x * v.x + mat[0].y * v.y + mat[0].z * v.z,
                    mat[1].x * v.x + mat[1].y * v.y + mat[1].z * v.z,
                    mat[2].x * v.x + mat[2].y * v.y + mat[2].z * v.z);
}

float getGGXPDF(struct PixelElement pixel, struct LightElement light){

}

float3 sampleCone(struct PixelElement pixel, struct LightElement light, float solid_angle){
    float r = rand() * light.rad;
    float theta = rand() * 2.0f * PI;
    float sqrt_r = sqrt(r);
    float3 dir = normalize(light.p - pixel.p);
    float3 sample_ray = normalize(rotate_to_frame(dir, normalize((float3)(cos(theta) * sqrt_r, sin(theta) * sqrt_r, d))));

    float3 diffuse = lambertian(pixel.diff_col, pixel.p, pixel.n, sample_ray);
    float3 specular = pixel.roughness < 1.0001f ? 
        cookTorrance(pixel.p, light.p, pixel.wo, pixel.n, pixel.roughness, pixel.eta, pixel.k, pixel.spec_col) : 
        (float3)(0.0f);
    float3 bsdf_col = diffuse + specular;

    float3 light_col = getLightCol(lights[i].type, lights[i].power, pixels[i].p, light[i].p, lights[i].n, lights[i].wi, light[i].roughness,
        lights[i].eta, lights[i].k, lights[i].spec_col);

    return bsdf_col * light_col;
}

float3 sampleBSDF(struct PixelElement pixel, struct LightElement light, float solid_angle, float cos_theta){
    float2 v = (float2)(rand(), rand());
    
}

__kernel void shadeVSL(__global const struct PixelElement* pixels, __global const struct LightElement* lights, 
    __global struct OutputElement *output, int num_pixels, int max_samples_per_vsl){
    int i = get_global_id(0);
    if(i >= num_pixels){
        return;
    }

    float coeff = (lights[i].coeff + 1.0f) / 2.0f;

    float central_disc_area = PI * lights[i].rad * lights[i].rad;
    float d = length(pixels[i].p - lights[i].p);
    float cos_theta = cos(asin(min(lights[i].rad / d, 1.0f)));
    float solid_angle = 2.0f * PI * (1.0f - cos_theta);
    int num_samples = max(1, (int)(sqrt(1.0f - cos_theta) * max_samples_per_vsl));

    float3 final_color = 0.0f;

    for(int i = 0; i < num_samples; ++i){
        final_color += sampleCone(pixels[i], lights[i], solid_angle);
        final_color += sampleBSDF(pixels[i], lights[i], solid_angle, cos_theta);
    }

    return final_color * coeff / (num_samples * central_disc_area);
}