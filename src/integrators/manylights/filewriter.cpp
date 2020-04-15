#include "filewriter.h"

#include "common.h"
#include <fstream>
#include <iostream>

MTS_NAMESPACE_BEGIN

ref<Film> createFilm(std::uint32_t width, std::uint32_t height, bool hdr){
    Properties props = hdr ? Properties("hdrfilm") : Properties("ldrfilm");
    props.setInteger("width", width);
    props.setInteger("height", height);
    props.setFloat("gamma", 2.2);
    props.setBoolean("banner", false);

    ref<Film> film = static_cast<Film*> (PluginManager::getInstance()->createObject(MTS_CLASS(Film), props));

    return film;
}

void writeOutputImage(Scene* scene, std::string filename, std::uint32_t width, std::uint32_t height, bool hdr, const std::vector<Spectrum>& data){
    assert(width * height == data.size());

    ref<Film> film = createFilm(width, height, hdr);

    Bitmap::EComponentFormat cf = hdr ? Bitmap::EComponentFormat::EUInt16 : Bitmap::EComponentFormat::EUInt8;
    Vector2i size; size.x = width; size.y = height;
    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, cf, size);

    for(std::uint32_t i = 0; i < data.size(); ++i){
        Point2i curr_pixel;
        curr_pixel.x = i % width;
        curr_pixel.y = i / width;

        output_bitmap->setPixel(curr_pixel, data[i]);
    }

    film->setBitmap(output_bitmap);
    fs::path scene_path = scene->getDestinationFile();
    film->setDestinationFile(scene_path.parent_path() / filename, 0);
    film->develop(scene, 0.f);
}

float writeOutputErrorImage(Scene* scene, std::string filename, std::uint32_t width, std::uint32_t height, bool hdr, const std::vector<Spectrum>& data1,
    const std::vector<Spectrum>& data2, float max_error){
    assert(data1.size() == data2.size());

    std::vector<Spectrum> error_col(data1.size(), Spectrum(0.f));
    float total_err = 0.f;

    for(std::uint32_t i = 0; i < error_col.size(); ++i){
        Spectrum dif = data1[i] - data2[i];
        float r, g, b;
        dif.toLinearRGB(r, g, b);
        float curr_err = std::abs(r) + std::abs(g) + std::abs(b);
        total_err += curr_err;
        float error = std::min(1.f, curr_err / max_error);
        std::tie(r, g, b) = floatToRGB(error);
        error_col[i].fromLinearRGB(r, g, b);
    }

    writeOutputImage(scene, filename, width, height, hdr, error_col);

    return total_err;
}

void writeOutputData(Scene* scene, std::string filename, bool new_file, const std::vector<float>& data, char delimiter){
    std::ofstream file;
    auto mode = new_file ? std::ios::out : std::ios::out | std::ios::app;
    fs::path scene_path = scene->getDestinationFile();
    std::string filename_and_path = (scene_path.parent_path() / filename).string();
    file.open(filename_and_path, mode);

    if(!new_file){
        file << std::endl;
    }

    for(std::uint32_t i = 0; i < data.size(); ++i){
        file << data[i];
        if(i < data.size() - 1){
            file << delimiter;
        }
    }
    
    file.close();
}

MTS_NAMESPACE_END