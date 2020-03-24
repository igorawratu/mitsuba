#ifndef hwshader_OCL
#define hwshader_OCL
const char *hwshader_ocl =
"struct PixelElement{\n"
"    float3 diff_col;\n"
"    float3 spec_col;\n"
"    float3 p;\n"
"    float3 n;\n"
"    float roughness;\n"
"    float3 eta;\n"
"    float3 k;\n"
"    float3 wo;\n"
"    int type;\n"
"    int slice_id;\n"
"    int intersected;\n"
"};\n"
"\n"
"struct LightElement{\n"
"    float3 power;\n"
"    float3 diff_col;\n"
"    float3 spec_col;\n"
"    float3 p;\n"
"    float3 n;\n"
"    float roughness;\n"
"    float3 eta;\n"
"    float3 k;\n"
"    float3 wi;\n"
"    float rad;\n"
"    int type;\n"
"    int light_surface_type;\n"
"};\n"
"\n"
"struct OutputElement{\n"
"    float3 col;\n"
"};\n"
"\n"
"struct CoeffElement{\n"
"    float coeff;\n"
"};\n"
"\n"
"#define PI 3.1415926\n"
"\n"
"float rand(float3 p3)\n"
"{\n"
"    float3 floor;\n"
"	p3 = fract(p3 * 0.1031f, &floor);\n"
"    p3 += dot(p3, p3.yzx + 33.33f);\n"
"    float out;\n"
"    return fract((p3.x + p3.y) * p3.z, &out);\n"
"}\n"
"\n"
"float3 rotate_to_frame(float3 old_frame, float3 new_frame, float3 v){\n"
"    if(length(old_frame - new_frame) < 0.00001f){\n"
"        return v;\n"
"    }\n"
"\n"
"    float3 axis = cross((float4)(old_frame, 0.0f), (float4)(new_frame, 0.0f)).xyz;\n"
"    float c = dot(old_frame, new_frame);\n"
"    float s = length(axis);\n"
"    axis = normalize(axis);\n"
"\n"
"    float invc = 1.0f - c;\n"
"\n"
"    float3 mat[3] = {(float3)(axis.x * axis.x * invc + c, axis.x * axis.y * invc - s * axis.z, axis.x * axis.z * invc + s * axis.y),\n"
"                    (float3)(axis.x * axis.y * invc + s * axis.z, axis.y * axis.y * invc + c, axis.y * axis.z * invc - s * axis.x),\n"
"                    (float3)(axis.x * axis.z * invc - s * axis.y, axis.y * axis.z * invc + s * axis.x, axis.z * axis.z * invc + c)};\n"
"\n"
"    return (float3)(mat[0].x * v.x + mat[0].y * v.y + mat[0].z * v.z,\n"
"                    mat[1].x * v.x + mat[1].y * v.y + mat[1].z * v.z,\n"
"                    mat[2].x * v.x + mat[2].y * v.y + mat[2].z * v.z);\n"
"}\n"
"\n"
"float lambertian(float3 n, float3 wi){\n"
"    float wi_dot_n = clamp(dot(wi, n), 0.0f, 1.0f);\n"
"\n"
"    return wi_dot_n / PI;\n"
"}\n"
"\n"
"float d_ggx(float3 h, float3 n, float alpha){\n"
"    float noh = clamp(dot(h, n), 0.0f, 1.0f);\n"
"    float alpha2 = alpha * alpha;\n"
"    float noh2 = noh * noh;\n"
"    float denom = 1.0f + noh2 * (alpha2 - 1.0f);\n"
"\n"
"    return alpha2 / (float)(PI * denom * denom);\n"
"}\n"
"\n"
"float g1_ggx_smith(float alpha, float3 n, float3 v){\n"
"    float alpha2 = alpha * alpha;\n"
"    float nov = clamp(dot(n, v), 0.00001f, 1.0f);\n"
"    float tanv = (1.0f - nov) / nov;\n"
"    \n"
"    return 2.0f / (1.0f + sqrt(1.0f + alpha2 * tanv * tanv));\n"
"}\n"
"\n"
"float g_ggx_smith(float alpha, float3 n, float3 v, float3 l){\n"
"    return g1_ggx_smith(alpha, n, v) * g1_ggx_smith(alpha, n, l);\n"
"}\n"
"\n"
"float3 fresnel(float costheta, float3 eta, float3 k) {\n"
"    float cos2theta = costheta * costheta;\n"
"    float sin2theta = 1.0f - cos2theta;\n"
"    float sin4theta = sin2theta * sin2theta;\n"
"\n"
"    float3 temp1 =  eta * eta - k * k - (float3)(sin2theta);\n"
"    float3 a2pb2 = sqrt(temp1 * temp1 + k * k * eta * eta * 4.0f);\n"
"    float3 a = sqrt((a2pb2 + temp1) * 0.5f);\n"
"\n"
"    float3 term1 = a2pb2 + (float3)(cos2theta);\n"
"    float3 term2 = a * (2.0f * costheta);\n"
"\n"
"    float3 rs2 = (term1 - term2) / (term1 + term2);\n"
"\n"
"    float3 term3 = a2pb2 * cos2theta + (float3)(sin4theta);\n"
"    float3 term4 = term2 * sin2theta;\n"
"\n"
"    float3 rp2 = rs2 * (term3 - term4) / (term3 + term4);\n"
"\n"
"    return 0.5f * (rp2 + rs2);\n"
"}\n"
"\n"
"float cookTorrance(float3 h, float3 wi, float3 wo, float3 n, float roughness){\n"
"    float d = d_ggx(h, n, roughness);\n"
"    float g = g_ggx_smith(roughness, n, wo, wi);\n"
"\n"
"    float nol = clamp(dot(h, wi), 0.f, 1.f);\n"
"    float nov = clamp(dot(wo, n), 0.f, 1.f);\n"
"\n"
"    return d * g / max(4.0f * nov, 0.0001f);\n"
"}\n"
"\n"
"float3 evalBSDF(float3 n, float3 wi, float3 wo, float3 eta, float3 k, float roughness, float3 spec_col,\n"
"    float3 diff_col){\n"
"    float3 h = normalize(wi + wo);\n"
"\n"
"    float3 f = fresnel(clamp(dot(h, wo), 0.0f, 1.0f), eta, k);\n"
"\n"
"    float3 specular = roughness < 1.0001f ? spec_col * f * cookTorrance(h, wi, wo, n, roughness) \n"
"        : (float3)(0.0f);\n"
"\n"
"    float3 diffuse = diff_col * lambertian(n, wi);\n"
"    if(roughness < 1.0001f){\n"
"        //float3 inveta2 = (float3)(1.0f) / (eta * eta);\n"
"        //only once since mitsuba does the other 1-f multiplication when looking up diffuse reflectance\n"
"        diffuse *= (1.0f - f) * (1.0f - f)/* * inveta2*/;\n"
"    }\n"
"\n"
"    return specular + diffuse;\n"
"}\n"
"\n"
"__kernel void shade(__global const struct PixelElement* pixels, \n"
"    __global const struct CoeffElement* coefficients,\n"
"    __global const struct LightElement* lights, \n"
"    __global struct OutputElement *output, \n"
"    int num_pixels, float min_dist, int clusters_per_slice, int curr_pass){\n"
"    int i = get_global_id(0);\n"
"    if(i >= num_pixels){\n"
"        return;\n"
"    }\n"
"\n"
"    if(pixels[i].intersected == 0){\n"
"        return;\n"
"    }\n"
"\n"
"    float coeff = 1.0f;//coefficients[i].coeff;\n"
"    int lidx = pixels[i].slice_id * clusters_per_slice + curr_pass;\n"
"\n"
"    float3 dir = lights[lidx].type == 0 ? -lights[lidx].n : lights[lidx].p - pixels[i].p;\n"
"    float3 wi = normalize(dir);\n"
"\n"
"    float3 bsdf_col = evalBSDF(pixels[i].n, wi, pixels[i].wo, pixels[i].eta, pixels[i].k, pixels[i].roughness,\n"
"        pixels[i].spec_col, pixels[i].diff_col);\n"
"\n"
"    float3 light_col;\n"
"    if(lights[lidx].type == 2){\n"
"        if(lights[lidx].light_surface_type == 0){\n"
"            light_col = clamp(dot(-wi, lights[lidx].n), 0.0f, 1.0f) / (float)PI;\n"
"        }\n"
"        else{\n"
"            light_col = evalBSDF(lights[lidx].n, -wi, lights[lidx].wi, \n"
"                lights[lidx].eta, lights[lidx].k, lights[lidx].roughness,\n"
"                lights[lidx].spec_col, lights[lidx].diff_col);\n"
"        }\n"
"\n"
"        light_col *= lights[lidx].power;\n"
"    }\n"
"    else{\n"
"        light_col = lights[lidx].power;\n"
"    }\n"
"\n"
"    if(lights[lidx].type != 0){\n"
"        float dist = fmax(min_dist, length(dir));\n"
"        light_col /= (dist * dist);\n"
"    }\n"
"\n"
"    output[i].col = output[i].col + bsdf_col * light_col * (float3)(coeff);\n"
"}\n"
"\n"
"float3 sampleCone(struct PixelElement pixel, struct LightElement light, float solid_angle, float3 seed){\n"
"    float r = rand(seed) * light.rad;\n"
"    float theta = rand(seed * 2) * 2.0f * PI;\n"
"    float sqrt_r = sqrt(r);\n"
"    float3 dir = light.p - pixel.p;\n"
"    float d = length(dir);\n"
"    dir = normalize(dir);\n"
"    float3 wi = (float3)(cos(theta) * sqrt_r, sin(theta) * sqrt_r, d);\n"
"    wi = normalize(rotate_to_frame((float3)(0.0f, 0.0f, 1.0f), dir, wi));\n"
"\n"
"    float prough = max(0.001f, pixel.roughness);\n"
"    float lrough = max(0.001f, light.roughness);\n"
"\n"
"    float3 bsdf_col = evalBSDF(pixel.n, wi, pixel.wo, pixel.eta, pixel.k, prough, pixel.spec_col, pixel.diff_col);\n"
"    float3 light_col;\n"
"    if(light.type == 2){\n"
"        if(light.light_surface_type == 0){\n"
"            light_col = clamp(dot(-wi, light.n), 0.0f, 1.0f) / (float)(PI) * light.diff_col;\n"
"        }\n"
"        else{\n"
"            light_col = evalBSDF(light.n, -wi, light.wi, \n"
"                light.eta, light.k, light.roughness,\n"
"                light.spec_col, light.diff_col);\n"
"        }\n"
"\n"
"        light_col *= light.power;\n"
"    }\n"
"    else{\n"
"        light_col = light.power;\n"
"    }\n"
"\n"
"    float3 h = normalize(wi + pixel.wo);\n"
"    float costheta = clamp(dot(h, wi), 0.0f, 1.0f);\n"
"    \n"
"    float spec_mag = prough < 1.000001f ? length(pixel.spec_col) : 0.0f;\n"
"    float dif_mag = length(pixel.diff_col);\n"
"    float mag_denom = spec_mag + dif_mag;\n"
"\n"
"    if(mag_denom < 0.0000001f){\n"
"        return 0.0f;\n"
"    }\n"
"\n"
"    float dif_rat = dif_mag / mag_denom;\n"
"\n"
"    float nol = clamp(dot(h, pixel.wo), 0.0f, 1.0f);\n"
"    float spec_prob = pixel.roughness < 1.0001f ? \n"
"        (1.0f - dif_rat) * g1_ggx_smith(pixel.roughness, pixel.wo, h) * d_ggx(h, pixel.n, prough) / max(4.0f * nol, 0.0001f) : 0.0f;\n"
"    float diff_prob = lambertian(pixel.n, wi) * dif_rat;\n"
"\n"
"    return bsdf_col * light_col / (spec_prob + diff_prob + 1.0f / solid_angle);\n"
"}\n"
"\n"
"float3 sampleBSDF(struct PixelElement pixel, struct LightElement light, float solid_angle, float ca_ctheta, float3 seed){\n"
"    float2 v = (float2)(rand(seed), rand(seed * 2));\n"
"    float spec_mag = pixel.roughness < 1.000001f ? length(pixel.spec_col) : 0.0f;\n"
"    float dif_mag = length(pixel.diff_col);\n"
"    float mag_denom = spec_mag + dif_mag;\n"
"\n"
"    if(mag_denom < 0.0000001f){\n"
"        return 0.0f;\n"
"    }\n"
"\n"
"    float diff_rat = dif_mag / mag_denom;\n"
"\n"
"    float3 dir = light.p - pixel.p;\n"
"    float d = length(dir);\n"
"    dir = normalize(dir);\n"
"\n"
"    float3 wi;\n"
"\n"
"    float prough = max(0.001f, pixel.roughness);\n"
"    float lrough = max(0.001f, light.roughness);\n"
"    \n"
"    if(v.x < diff_rat){\n"
"        v.x /= diff_rat;\n"
"        float theta = 2.0f * PI * v.y;\n"
"        float r = sqrt(v.x);\n"
"        wi = normalize((float3)(r * cos(theta), r * sin(theta), 1.0f - v.x));\n"
"        wi = rotate_to_frame((float3)(0.0f, 0.0f, 1.0f), pixel.n, wi);\n"
"    }\n"
"    else{\n"
"        v.x = (v.x - diff_rat) / (1.0f - diff_rat);\n"
"        float numer = prough * sqrt(v.x);\n"
"        float denom = sqrt(1.0f - v.x);\n"
"        float phi = atan2(numer, denom);\n"
"        float theta = 2.0f * PI * v.y;\n"
"\n"
"        float sin_theta = sin(theta);\n"
"        float sin_phi = sin(phi);\n"
"        float cos_phi = cos(phi);\n"
"        float cos_theta = cos(theta);\n"
"        float3 wm = normalize((float3)(cos_theta * sin_phi, sin_theta * sin_phi, cos_phi));\n"
"        wm = rotate_to_frame((float3)(0.0f, 0.0f, 1.0f), pixel.n, wm);\n"
"\n"
"        wi = 2.0f * dot(pixel.wo, wm) * wm - pixel.wo;\n"
"    }\n"
"\n"
"    \n"
"    float dp = dot(dir, wi);\n"
"\n"
"    if(dp > ca_ctheta)\n"
"    {\n"
"        float3 bsdf_col = evalBSDF(pixel.n, wi, pixel.wo, pixel.eta, pixel.k, prough, pixel.spec_col, pixel.diff_col);\n"
"        float3 light_col;\n"
"        if(light.type == 2){\n"
"            if(light.light_surface_type == 0){\n"
"                light_col = clamp(dot(-wi, light.n), 0.0f, 1.0f) / (float)(PI) * light.diff_col;\n"
"            }\n"
"            else{\n"
"                light_col = evalBSDF(light.n, -wi, light.wi, \n"
"                    light.eta, light.k, light.roughness,\n"
"                    light.spec_col, light.diff_col);\n"
"            }\n"
"\n"
"            light_col *= light.power;\n"
"        }\n"
"        else{\n"
"            light_col = light.power;\n"
"        }\n"
"\n"
"        float3 h = normalize(wi + pixel.wo);\n"
"        float nol = clamp(dot(h, wi), 0.0f, 1.0f);\n"
"        float spec_prob = prough < 1.0001f ? \n"
"            (1.0f - diff_rat) * g1_ggx_smith(pixel.roughness, pixel.wo, h) * d_ggx(h, pixel.n, prough) / max(4.0f * nol, 0.0001f) : 0.0f;\n"
"        float diff_prob = lambertian(pixel.n, wi) * diff_rat;\n"
"\n"
"        return bsdf_col * light_col / (spec_prob + diff_prob + 1.0f / solid_angle);\n"
"    }\n"
"\n"
"    return 0.0f;\n"
"}\n"
"\n"
"__kernel void shadeVSL(__global const struct PixelElement* pixels, \n"
"    __global const struct CoeffElement* coefficients,\n"
"    __global const struct LightElement* lights, \n"
"    __global struct OutputElement *output, \n"
"    int num_pixels, float min_dist, int clusters_per_slice, int curr_pass, int max_samples_per_vsl){\n"
"    int i = get_global_id(0);\n"
"    if(i >= num_pixels){\n"
"        return;\n"
"    }\n"
"\n"
"    if(pixels[i].intersected == 0){\n"
"        return;\n"
"    }\n"
"\n"
"    float3 seed = /*pixels[i].p / min_dist * 100.f * */(float)(curr_pass + 1.0f);\n"
"\n"
"    float coeff = coefficients[i].coeff;// > 0 ? 1.0f : 0.0f;\n"
"    int curr_light_idx = pixels[i].slice_id * clusters_per_slice + curr_pass;\n"
"\n"
"    float dp = dot(normalize(lights[curr_light_idx].p - pixels[i].p), pixels[i].n);\n"
"    if(dp < 0.0001f){\n"
"        return;\n"
"    }\n"
"\n"
"    float central_disc_area = PI * lights[curr_light_idx].rad * lights[curr_light_idx].rad;\n"
"    float d = length(pixels[i].p - lights[curr_light_idx].p);\n"
"    float hypot_length = sqrt(d * d + lights[curr_light_idx].rad * lights[curr_light_idx].rad);\n"
"    float cos_theta = d / hypot_length;\n"
"    //float cos_theta = cos(asin(min(lights[curr_light_idx].rad / d, 1.0f)));\n"
"    float solid_angle = 2.0f * PI * (1.0f - cos_theta);\n"
"    int num_samples = max(1, (int)(sqrt(1.0f - cos_theta) * max_samples_per_vsl));\n"
"\n"
"    float3 final_color = 0.0f;\n"
"\n"
"    for(int sample = 0; sample < num_samples; ++sample){\n"
"        final_color += sampleCone(pixels[i], lights[curr_light_idx], solid_angle, seed * sample);\n"
"        final_color += sampleBSDF(pixels[i], lights[curr_light_idx], solid_angle, cos_theta, seed * sample);\n"
"    }\n"
"\n"
"    output[i].col = output[i].col + final_color * coeff / (num_samples * central_disc_area);\n"
"}\n";
#endif
