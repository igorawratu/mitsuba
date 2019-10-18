struct PixelElement{
    float3 diff_col;
    float3 spec_col;
    float3 p;
    float3 n;
    float roughness;
    float3 eta;
    float3 k;
    float3 wi;
};

struct LightElement{
    float3 col;
    float3 n;
    float3 p;
    float rad;
    float coeff;
    int type;
};

struct OutputElement{
    float3 col;
};

#define PI 3.1415926

float3 lambertian(float3 albedo, float3 p, float3 n, float3 lp, float3 lcol, float3 ln){
    float3 dir = lp - p;
    float3 wi = normalize(dir);
    float wi_dot_n = fmax(0.0f, dot(wi, n));
    
    float3 col = albedo * (float3)(wi_dot_n) * lcol / (float3)(PI);

    return col;
}

float ggx_dist(float3 h, float3 n, float alpha){
    float noh = clamp(dot(h, n), 0.0f, 1.0f);
    float alpha2 = alpha * alpha;

    float denom = (noh * noh) * (alpha2 - 1.0f) + 1.0f;

    return alpha2 / (denom * denom * PI);
}

float schlick_geom(float alpha, float3 n, float3 v, float3 l){
    float k = alpha + 1.0f;
    k = (k * k) / 8.0f;
    
    float nov = clamp(dot(n, v), 0.0f, 1.0f);
    float nol = clamp(dot(n, l), 0.0f, 1.0f);

    float g1 = nov / (nov * (1.0f - k) + k);
    float g2 = nol / (nol * (1.0f - k) + k);

    return g1 * g2;
}

float3 schlick_fres(float3 v, float3 h, float3 f0){
    float voh = 1.0 - clamp(dot(v, h), 0.0f, 1.0f);
    voh = powr(voh, 5.0f);
    return f0 + ((float3)(1.0f) - f0) * (float3)(voh);
}

float3 cookTorrance(struct PixelElement pixel, struct LightElement light){
    float3 l = normalize(light.p - pixel.p);
    float3 h = normalize(pixel.wi + l);
    float3 n = normalize(pixel.n);
    
    float d = ggx_dist(h, n, pixel.roughness);
    float g = schlick_geom(pixel.roughness, n, pixel.wi, l);
    float3 f = schlick_fres(pixel.wi, h, pixel.spec_col);

    float nol = clamp(dot(n, l), 0.f, 1.f);
    float nov = clamp(dot(pixel.wi, n), 0.f, 1.f);

    return f * d * g / clamp(4.0f * (nol * nov + 0.05f), 0.f, 1.f);
}

__kernel void shade(__global const struct PixelElement* pixels, __global const struct LightElement* lights, 
    __global struct OutputElement *output, int num_pixels, float min_dist){
    int i = get_global_id(0);
    if(i >= num_pixels){
        return;
    }

    float coeff = (lights[i].coeff + 1.0f) / 2.0f;

    float3 diffuse = lambertian(pixels[i].diff_col, pixels[i].p, pixels[i].n, lights[i].p, lights[i].col, 
        lights[i].n);

    float3 specular = pixels[i].roughness < 1.0001f ? cookTorrance(pixels[i], lights[i]) : (float3)(0.0f);

    float3 color = diffuse + specular;

    float3 light_col = lights[i].col;

    float3 dir = pixels[i].p - lights[i].p;
    if(lights[i].type != 0){
        float dist = fmax(min_dist, length(dir));
        light_col /= (dist * dist);
    }

    if(lights[i].type == 2){
        float lnodir = fmax(0.0f, dot(normalize(dir), lights[i].n));
        light_col *= (float3)(lnodir) / (float3)(PI);
    }
    
    color *= light_col;

    //output[i].col = output[i].col + color * (float3)(coeff);

    output[i].col = output[i].col + color * (float3)(coeff);
}