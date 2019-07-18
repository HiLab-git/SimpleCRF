add_executable( dense3DCrfInferenceOnNiis dense3DCrfInferenceOnNiis.cpp )
target_link_libraries( dense3DCrfInferenceOnNiis densecrf ${ITK_LIBRARIES})
