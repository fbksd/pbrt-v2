
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_SERVERRENDERER_H
#define PBRT_RENDERERS_SERVERRRENDERER_H

#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"


class PBRTServer;
class SceneInfo;
class SampleLayout;
class SampleQuerySplitBuffers;
class CropWindow;
class SamplesPipe;


class ServerRenderer : public Renderer {
public:
    ServerRenderer(Sampler *s, Camera *c, SurfaceIntegrator *si, VolumeIntegrator *vi, bool visIds);

    ~ServerRenderer();

    void Render(const Scene *scene);

    Spectrum Li(const Scene *scene, const RayDifferential &ray,
                const Sample *sample, RNG &rng, MemoryArena &arena,
                Intersection *isect = NULL, Spectrum *T = NULL, SampleBuffer* sampleBuffer = nullptr) const;

    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
                           const Sample *sample, RNG &rng, MemoryArena &arena) const;


    void getSceneInfo(SceneInfo *scene);
    void evaluateSamples(bool isSPP, int numSamples, int* resultSize);
    void evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow& crop, int* resultSize);
    void evaluateSamplesPDF(bool isSPP, int numSamples, const float* pdf, int* resultSize);

    void run(Sampler* sampler, Sample* origSample, SamplesPipe& pipe);

private:
    bool visualizeObjectIds;
    Sampler *sampler;
    Camera *camera;
    Sample *mainSample;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
    const Scene* scene;
    int w, h, maxSPP;
};



/**
 *  Discret Probability Density Function sampler.
 *
 *  Permite gerar amostras de acordo com uma pdf por partes dada como entrada no construtor.
 *  Cada valor retornado em sampler() está no intervalo [ 0, fv.size() ).
 *  O custo de sample() é \f$ O(log(n)) \f$ no tamanho de fv e o construtor tem custo linear.
 */
class DPDF
{
    friend class DPDF2;
public:
    DPDF(const float* fv, unsigned int n);
    DPDF(const std::vector<float>& fv);

    //! Gera um índice aleatório.
    unsigned int sample(const RNG& rng, float* pdf = NULL);

private:
    std::vector<float> m_pdf;
    std::vector<float> m_cpdf;
    float m_funcInt;
};


class DPDF2
{
public:
    DPDF2(const float* fv, int nx, int ny);
    ~DPDF2();

    //! Gera um índice aleatório.
    void sample(const RNG& rng, int* x, int* y, float* pdf = NULL);

private:
    int m_nx, m_ny;
    std::vector<DPDF*> m_dpdfX;
    DPDF* m_marginal;
};



#endif // PBRT_RENDERERS_SAMPLERRENDERER_H
