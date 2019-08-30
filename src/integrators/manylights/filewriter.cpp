#include "filewriter.h"
#include "common.h"
#include <fstream>
#include <iostream>

MTS_NAMESPACE_BEGIN

ref<Film> createFilm(std::uint32_t width, std::uint32_t height, bool hdr){
    Properties props;
    props.setInteger("width", width);
    props.setInteger("height", height);
    props.setFloat("gamma", 2.2);

    ref<Film> film(hdr ? new HDRFilm(props) : new LDRFilm(props));

    return film;
}

void writeOutputImage(Scene* scene, std::string filename, std::uint32_t width, std::uint32_t height, bool hdr, const std::vector<Spectrum>& data){
    assert(width * height == data.size());

    ref<Film> film = createFilm(width, height, hdr);

    Bitmap::EComponentFormat cf = hdr ? EUint16 : EUint8;
    Vector2i size; size.x = width; size.y = height;
    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, cf, size);

    for(std::uint32_t i = 0; i < data.size(); ++i){
        Vector2i curr_pixel;
        curr_pixel.x = i % width;
        curr_pixel.y = i / width;

        output_bitmap.setPixel(curr_pixel, data[i]);
    }

    film->setBitmap(output_bitmap);
    film->setDestinationFile(fs::path(filename), 0);
    film->develop(scene, 0.f);
}

void writeOutputErrorImage(Scene* scene, std::string filename, std::uint32_t width, std::uint32_t height, bool hdr, const std::vector<Spectrum>& data1,
    const std::vector<Spectrum>& data2, float max_error){
    assert(data1.size() == data2.size());

    std::vector<Spectrum> error_col(data1.size());

    for(std::uint32_t i = 0; i < error_col.size(); ++i){
        Spectrum dif = data1[i] - data2[i];
        float r, g, b;
        dif.toLinearRGB(r, g, b);
        float error = std::max(1.f, (std::abs(r) + std::abs(g) + std::abs(b)) / max_error);
        std::tie(r, g, b) = floatToRGB(error);
        error_col[i].fromLinearRGB(r, g, b);
    }

    writeOutputImage(scene, filename, width, height, hdr, error_col);
}

void writeOutputData(std::string filename, bool new_file, const std::vector<float>& data, char delimiter){
    ofstream file;
    auto mode = new_file ? ios::out : ios::out | ios::app;
    file.open(filename, mode);

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