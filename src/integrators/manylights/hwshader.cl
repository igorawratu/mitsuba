struct PixelElement{
    float3 diff_col;
    float3 spec_col;
    float3 p;
    float3 n;
    float roughness;
    float3 eta;
    float3 k;
    float3 wo;
    int type;
    int slice_id;
    int intersected;
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
    float rad;
    int type;
    int light_surface_type;
};

struct OutputElement{
    float3 col;
};

struct CoeffElement{
    float coeff;
};

#define PI 3.1415926

float rand(float3 p3)
{
    float3 floor;
	p3 = fract(p3 * 0.1031f, &floor);
    p3 += dot(p3, p3.yzx + 33.33f);
    float out;
    return fract((p3.x + p3.y) * p3.z, &out);
}

float3 rotate_to_frame(float3 old_frame, float3 new_frame, float3 v){
    if(length(old_frame - new_frame) < 0.00001f){
        return v;
    }

    float3 axis = cross((float4)(old_frame, 0.0f), (float4)(new_frame, 0.0f)).xyz;
    float c = dot(old_frame, new_frame);
    float s = length(axis);
    axis = normalize(axis);

    float invc = 1.0f - c;

    float3 mat[3] = {(float3)(axis.x * axis.x * invc + c, axis.x * axis.y * invc - s * axis.z, axis.x * axis.z * invc + s * axis.y),
                    (float3)(axis.x * axis.y * invc + s * axis.z, axis.y * axis.y * invc + c, axis.y * axis.z * invc - s * axis.x),
                    (float3)(axis.x * axis.z * invc - s * axis.y, axis.y * axis.z * invc + s * axis.x, axis.z * axis.z * invc + c)};

    return (float3)(mat[0].x * v.x + mat[0].y * v.y + mat[0].z * v.z,
                    mat[1].x * v.x + mat[1].y * v.y + mat[1].z * v.z,
                    mat[2].x * v.x + mat[2].y * v.y + mat[2].z * v.z);
}

float lambertian(float3 n, float3 wi){
    float wi_dot_n = clamp(dot(wi, n), 0.0f, 1.0f);

    return wi_dot_n / PI;
}

float d_ggx(float3 h, float3 n, float alpha){
    float noh = clamp(dot(h, n), 0.0f, 1.0f);
    float alpha2 = alpha * alpha;
    float noh2 = noh * noh;
    float denom = 1.0f + noh2 * (alpha2 - 1.0f);

    return alpha2 / (float)(PI * denom * denom);
}

float g1_ggx_smith(float alpha, float3 n, float3 v){
    float alpha2 = alpha * alpha;
    float nov = clamp(dot(n, v), 0.00001f, 1.0f);
    float tanv = (1.0f - nov) / nov;
    
    return 2.0f / (1.0f + sqrt(1.0f + alpha2 * tanv * tanv));
}

float g_ggx_smith(float alpha, float3 n, float3 v, float3 l){
    return g1_ggx_smith(alpha, n, v) * g1_ggx_smith(alpha, n, l);
}

float3 fresnel(float costheta, float3 eta, float3 k) {
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

float cookTorrance(float3 h, float3 wi, float3 wo, float3 n, float roughness){
    float d = d_ggx(h, n, roughness);
    float g = g_ggx_smith(roughness, n, wo, wi);

    float nol = clamp(dot(h, wi), 0.f, 1.f);
    float nov = clamp(dot(wo, n), 0.f, 1.f);

    return d * g / max(4.0f * nov, 0.0001f);
}

float3 evalBSDF(float3 n, float3 wi, float3 wo, float3 eta, float3 k, float roughness, float3 spec_col,
    float3 diff_col){
    float3 h = normalize(wi + wo);

    float3 f = fresnel(clamp(dot(h, wo), 0.0f, 1.0f), eta, k);

    float3 specular = roughness < 1.0001f ? spec_col * f * cookTorrance(h, wi, wo, n, roughness) 
        : (float3)(0.0f);

    float3 diffuse = diff_col * lambertian(n, wi);
    if(roughness < 1.0001f){
        //float3 inveta2 = (float3)(1.0f) / (eta * eta);
        //only once since mitsuba does the other 1-f multiplication when looking up diffuse reflectance
        diffuse *= (1.0f - f) * (1.0f - f)/* * inveta2*/;
    }

    return specular + diffuse;
}

__kernel void shade(__global const struct PixelElement* pixels, 
    __global const struct CoeffElement* coefficients,
    __global const struct LightElement* lights, 
    __global struct OutputElement *output, 
    int num_pixels, float min_dist, int clusters_per_slice, int curr_pass){
    int i = get_global_id(0);
    if(i >= num_pixels){
        return;
    }

    if(pixels[i].intersected == 0){
        return;
    }

    float coeff = coefficients[i].coeff;// > 0 ? 1.0f : 0.0f;
    int lidx = pixels[i].slice_id * clusters_per_slice + curr_pass;

    float3 dir = lights[lidx].type == 0 ? -lights[lidx].n : lights[lidx].p - pixels[i].p;
    float3 wi = normalize(dir);

    float3 bsdf_col = evalBSDF(pixels[i].n, wi, pixels[i].wo, pixels[i].eta, pixels[i].k, pixels[i].roughness,
        pixels[i].spec_col, pixels[i].diff_col);

    float3 light_col;
    if(lights[lidx].type == 2){
        if(lights[lidx].light_surface_type == 0){
            light_col = clamp(dot(-wi, lights[lidx].n), 0.0f, 1.0f) / (float)PI;
        }
        else{
            light_col = evalBSDF(lights[lidx].n, -wi, lights[lidx].wi, 
                lights[lidx].eta, lights[lidx].k, lights[lidx].roughness,
                lights[lidx].spec_col, lights[lidx].diff_col);
        }

        light_col *= lights[lidx].power;
    }
    else{
        light_col = lights[lidx].power;
    }

    if(lights[lidx].type != 0){
        float dist = fmax(min_dist, length(dir));
        light_col /= (dist * dist);
    }

    output[i].col = output[i].col + bsdf_col * light_col * (float3)(coeff);
}

float3 sampleCone(struct PixelElement pixel, struct LightElement light, float solid_angle, float3 seed){
    float r = rand(seed) * light.rad;
    float theta = rand(seed * 2) * 2.0f * PI;
    float sqrt_r = sqrt(r);
    float3 dir = light.p - pixel.p;
    float d = length(dir);
    dir = normalize(dir);
    float3 wi = (float3)(cos(theta) * sqrt_r, sin(theta) * sqrt_r, d);
    wi = normalize(rotate_to_frame((float3)(0.0f, 0.0f, 1.0f), dir, wi));

    float prough = max(0.001f, pixel.roughness);
    float lrough = max(0.001f, light.roughness);

    float3 bsdf_col = evalBSDF(pixel.n, wi, pixel.wo, pixel.eta, pixel.k, prough, pixel.spec_col, pixel.diff_col);
    float3 light_col;
    if(light.type == 2){
        if(light.light_surface_type == 0){
            light_col = clamp(dot(-wi, light.n), 0.0f, 1.0f) / (float)(PI) * light.diff_col;
        }
        else{
            light_col = evalBSDF(light.n, -wi, light.wi, 
                light.eta, light.k, light.roughness,
                light.spec_col, light.diff_col);
        }

        light_col *= light.power;
    }
    else{
        light_col = light.power;
    }

    float3 h = normalize(wi + pixel.wo);
    float costheta = clamp(dot(h, wi), 0.0f, 1.0f);
    
    float spec_mag = prough < 1.000001f ? length(pixel.spec_col) : 0.0f;
    float dif_mag = length(pixel.diff_col);
    float mag_denom = spec_mag + dif_mag;

    if(mag_denom < 0.0000001f){
        return 0.0f;
    }

    float dif_rat = dif_mag / mag_denom;

    float nol = clamp(dot(h, pixel.wo), 0.0f, 1.0f);
    float spec_prob = pixel.roughness < 1.0001f ? 
        (1.0f - dif_rat) * g1_ggx_smith(pixel.roughness, pixel.wo, h) * d_ggx(h, pixel.n, prough) / max(4.0f * nol, 0.0001f) : 0.0f;
    float diff_prob = lambertian(pixel.n, wi) * dif_rat;

    return bsdf_col * light_col / (spec_prob + diff_prob + 1.0f / solid_angle);
}

float3 sampleBSDF(struct PixelElement pixel, struct LightElement light, float solid_angle, float ca_ctheta, float3 seed){
    float2 v = (float2)(rand(seed), rand(seed * 2));
    float spec_mag = pixel.roughness < 1.000001f ? length(pixel.spec_col) : 0.0f;
    float dif_mag = length(pixel.diff_col);
    float mag_denom = spec_mag + dif_mag;

    if(mag_denom < 0.0000001f){
        return 0.0f;
    }

    float diff_rat = dif_mag / mag_denom;

    float3 dir = light.p - pixel.p;
    float d = length(dir);
    dir = normalize(dir);

    float3 wi;

    float prough = max(0.001f, pixel.roughness);
    float lrough = max(0.001f, light.roughness);
    
    if(v.x < diff_rat){
        v.x /= diff_rat;
        float theta = 2.0f * PI * v.y;
        float r = sqrt(v.x);
        wi = normalize((float3)(r * cos(theta), r * sin(theta), 1.0f - v.x));
        wi = rotate_to_frame((float3)(0.0f, 0.0f, 1.0f), pixel.n, wi);
    }
    else{
        v.x = (v.x - diff_rat) / (1.0f - diff_rat);
        float numer = prough * sqrt(v.x);
        float denom = sqrt(1.0f - v.x);
        float phi = atan2(numer, denom);
        float theta = 2.0f * PI * v.y;

        float sin_theta = sin(theta);
        float sin_phi = sin(phi);
        float cos_phi = cos(phi);
        float cos_theta = cos(theta);
        float3 wm = normalize((float3)(cos_theta * sin_phi, sin_theta * sin_phi, cos_phi));
        wm = rotate_to_frame((float3)(0.0f, 0.0f, 1.0f), pixel.n, wm);

        wi = 2.0f * dot(pixel.wo, wm) * wm - pixel.wo;
    }

    
    float dp = dot(dir, wi);

    if(dp > ca_ctheta)
    {
        float3 bsdf_col = evalBSDF(pixel.n, wi, pixel.wo, pixel.eta, pixel.k, prough, pixel.spec_col, pixel.diff_col);
        float3 light_col;
        if(light.type == 2){
            if(light.light_surface_type == 0){
                light_col = clamp(dot(-wi, light.n), 0.0f, 1.0f) / (float)(PI) * light.diff_col;
            }
            else{
                light_col = evalBSDF(light.n, -wi, light.wi, 
                    light.eta, light.k, light.roughness,
                    light.spec_col, light.diff_col);
            }

            light_col *= light.power;
        }
        else{
            light_col = light.power;
        }

        float3 h = normalize(wi + pixel.wo);
        float nol = clamp(dot(h, wi), 0.0f, 1.0f);
        float spec_prob = prough < 1.0001f ? 
            (1.0f - diff_rat) * g1_ggx_smith(pixel.roughness, pixel.wo, h) * d_ggx(h, pixel.n, prough) / max(4.0f * nol, 0.0001f) : 0.0f;
        float diff_prob = lambertian(pixel.n, wi) * diff_rat;

        return bsdf_col * light_col / (spec_prob + diff_prob + 1.0f / solid_angle);
    }

    return 0.0f;
}

__kernel void shadeVSL(__global const struct PixelElement* pixels, 
    __global const struct CoeffElement* coefficients,
    __global const struct LightElement* lights, 
    __global struct OutputElement *output, 
    int num_pixels, float min_dist, int clusters_per_slice, int curr_pass, int max_samples_per_vsl){
    int i = get_global_id(0);
    if(i >= num_pixels){
        return;
    }

    if(pixels[i].intersected == 0){
        return;
    }

    float3 seed = /*pixels[i].p / min_dist * 100.f * */(float)(curr_pass + 1.0f);

    float coeff = coefficients[i].coeff;// > 0 ? 1.0f : 0.0f;
    int curr_light_idx = pixels[i].slice_id * clusters_per_slice + curr_pass;

    float dp = dot(normalize(lights[curr_light_idx].p - pixels[i].p), pixels[i].n);
    if(dp < 0.0001f){
        return;
    }

    float central_disc_area = PI * lights[curr_light_idx].rad * lights[curr_light_idx].rad;
    float d = length(pixels[i].p - lights[curr_light_idx].p);
    float hypot_length = sqrt(d * d + lights[curr_light_idx].rad * lights[curr_light_idx].rad);
    float cos_theta = d / hypot_length;
    //float cos_theta = cos(asin(min(lights[curr_light_idx].rad / d, 1.0f)));
    float solid_angle = 2.0f * PI * (1.0f - cos_theta);
    int num_samples = max(1, (int)(sqrt(1.0f - cos_theta) * max_samples_per_vsl));

    float3 final_color = 0.0f;

    for(int sample = 0; sample < num_samples; ++sample){
        final_color += sampleCone(pixels[i], lights[curr_light_idx], solid_angle, seed * sample);
        final_color += sampleBSDF(pixels[i], lights[curr_light_idx], solid_angle, cos_theta, seed * sample);
    }

    output[i].col = output[i].col + final_color * coeff / (num_samples * central_disc_area);
}