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

float3 cookTorrance(struct PixelElement pixel, struct LightElement light){
    float3 l = normalize(light.p - pixel.p);
    float3 h = normalize(pixel.wo + l);
    float3 n = normalize(pixel.n);
    
    float d = d_ggx(h, n, pixel.roughness);
    float g = g_ggx_smith(pixel.roughness, n, pixel.wo, l);
    float3 f = fresnel(pixel.wo, h, pixel.eta, pixel.k) * pixel.spec_col;

    float nol = clamp(dot(n, l), 0.f, 1.f);
    float nov = clamp(dot(pixel.wo, n), 0.f, 1.f);

    return f * d * g / clamp(4.0f * (nol * nov + 0.05f), 0.f, 1.f);
}

float3 getLightCol(struct LightElement light, struct PixelElement pixel){
    if(light.type != 2){
        return light.power;
    }

    float3 dir = normalize(pixel.p - light.p);
    
    float nodir = clamp(dot(light.n, dir), 0.0f, 1.0f);
    float3 diffuse_col = light.power * nodir / (float)(PI);

    float3 h = normalize(dir + light.wi);
    
    float d = d_ggx(h, light.n, light.roughness);
    float g = g_ggx_smith(light.roughness, light.n, dir, light.wi);
    float3 f = fresnel(dir, h, light.eta, light.k) * light.power * light.spec_col;

    float nol = clamp(dot(light.n, light.wi), 0.f, 1.f);
    float nov = clamp(dot(dir, light.n), 0.f, 1.f);

    float3 specular_col = light.roughness < 1.000001f ? f * d * g / clamp(4.0f * (nol * nov + 0.05f), 0.f, 1.f) : (float3)(0.0f);

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

    float3 specular = pixels[i].roughness < 1.0001f ? cookTorrance(pixels[i], lights[i]) : (float3)(0.0f);

    float3 bsdf_color = diffuse + specular;

    float3 light_col = getLightCol(lights[i], pixels[i]);

    if(lights[i].type != 0){
        float3 dir = pixels[i].p - lights[i].p;
        float dist = fmax(min_dist, length(dir));
        light_col /= (dist * dist);
    }

    output[i].col = output[i].col + bsdf_color * light_col * (float3)(coeff);
}