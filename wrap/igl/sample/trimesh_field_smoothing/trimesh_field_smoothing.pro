DEPENDPATH += . ../../../../..
INCLUDEPATH += . ../../../.. ../../../../eigenlib ../../../../include

CONFIG += console c++11
TEMPLATE = app
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle

QMAKE_CXXFLAGS += -std=c++11
TARGET = trimesh_field_smoothing
SOURCES += cube_style_precomputation.cpp \
    cube_style_single_iteration.cpp \
    energy_computing.cpp \
    fit_rotations_l1.cpp \
    main_test.cpp \
    normalize_unitbox.cpp \
    orthogonal_procrustes.cpp \
    shrinkage.cpp \
    ../../../../wrap/ply/plylib.cpp

HEADERS += \
    conversionMeshes.h \
    cube_style_data.h \
    cube_style_precomputation.h \
    cube_style_single_iteration.h \
    cubic_stylizing.h \
    edge_flip_cubization.h \
    energy_computing.h \
    exporter_cubic.h \
    fit_rotations_l1.h \
    normalize_unitbox.h \
    orthogonal_procrustes.h \
    shrinkage.h
