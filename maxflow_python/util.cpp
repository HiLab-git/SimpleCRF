#include "util.h"
#include <iostream>
#include <cmath>

// for 2d images
template<typename T>
T get_pixel(const T * data, int height, int width, int h, int w)
{
    return data[h * width + w];
}

template<typename T>
std::vector<T> get_pixel_vector(const T * data, int height, int width, int channel, int h, int w)
{
    std::vector<T> pixel_vector(channel);
    for (int c = 0; c < channel; c++){
        pixel_vector[c]= data[h * width * channel + w * channel + c];
    }
    return pixel_vector;
}

template<typename T>
void set_pixel(T * data,  int height, int width, int h, int w, T value)
{
    data[h * width + w] = value;
}

template<typename T>
float get_l2_distance(std::vector<T> p1, std::vector<T> p2)
{
    T sq_sum = 0.0;
    for(int d = 0; d < p1.size(); d++)
    {
        sq_sum = sq_sum + (p1[d] - p2[d]) * (p1[d] - p2[d]);
    }
    float dis = sqrt(sq_sum);
    return dis;
}

template
float get_pixel<float>(const float * data, int height, int width, int h, int w);

template
int get_pixel<int>(const int * data, int height, int width, int h, int w);

template
std::vector<float> get_pixel_vector<float>(const float * data, int height, int width, int channel, int h, int w);

template
unsigned char get_pixel<unsigned char>(const unsigned char * data, int height, int width, int h, int w);


template
void set_pixel<float>(float * data, int height, int width, int h, int w, float value);

template
void set_pixel<int>(int * data, int height, int width, int h, int w, int value);

template
void set_pixel<unsigned char>(unsigned char * data, int height, int width, int h, int w, unsigned char value);

template
float get_l2_distance(std::vector<float> p1, std::vector<float> p2);


// for 3d images
template<typename T>
T get_pixel(const T * data, int depth, int height, int width, int d, int h, int w)
{
    return data[(d*height + h) * width + w];
}

template<typename T>
std::vector<T> get_pixel_vector(const T * data, int depth,  int height, int width, int channel, int d, int h, int w)
{
    std::vector<T> pixel_vector(channel);
    for (int c = 0; c < channel; c++){
        pixel_vector[c]= data[d*height*width*channel +  h * width * channel + w * channel + c];
    }
    return pixel_vector;
}

template<typename T>
void set_pixel(T * data, int depth, int height, int width, int d, int h, int w, T value)
{
    data[(d*height + h) * width + w] = value;
}

template
float get_pixel<float>(const float * data, int depth, int height, int width, int d, int h, int w);

template
std::vector<float> get_pixel_vector<float>(const float * data, int depth, int height, int width, int channel, int d, int h, int w);

template
int get_pixel<int>(const int * data, int depth, int height, int width, int d, int h, int w);

template
unsigned char get_pixel<unsigned char>(const unsigned char * data,
                                       int depth, int height, int width,
                                       int d, int h, int w);


template
void set_pixel<float>(float * data, int depth, int height, int width, int d, int h, int w, float value);

template
void set_pixel<int>(int * data, int depth, int height, int width, int d, int h, int w, int value);

template
void set_pixel<unsigned char>(unsigned char * data,
                              int depth, int height, int width,
                              int d, int h, int w, unsigned char value);

