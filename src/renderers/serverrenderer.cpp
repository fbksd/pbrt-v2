
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


#include "renderers/serverrenderer.h"
#include "Benchmark/RenderingServer/RenderingServer.h"
#include "stdafx.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"

#include <QEventLoop>


// Global sample parameters variable. Declared as 'extern' in pbrt.h and filled through the system.
extern string g_samplerName;


static uint32_t hash(char *key, uint32_t len)
{
    uint32_t hash = 0, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}


static bool isPerfectSquare(int n)
{
    if (n < 0)
        return false;
    int root(round(sqrt(n)));
    return n == root * root;
}


//======================================================================
//                    AUXILIAR SAMPLERS
//======================================================================

#include "samplers/adaptive.h"
#include "samplers/bestcandidate.h"
#include "samplers/halton.h"
#include "samplers/lowdiscrepancy.h"
#include "samplers/random.h"
#include "samplers/stratified.h"

class SparseSampler: public Sampler
{
public:
    SparseSampler(int xstart, int xend, int ystart, int yend, int ns, float sopen, float sclose):
        Sampler(xstart, xend, ystart, yend, ns, sopen, sclose)
    {
        m_numSamples = ns;
        m_sampler = new RandomSampler(0, 1, 0, 1, m_numSamples, sopen, sclose);
    }

    ~SparseSampler()
    { delete m_sampler; }

    int MaximumSampleCount() { return 1; }

    int GetMoreSamples(Sample *sample, RNG &rng)
    {
        int result = m_sampler->GetMoreSamples(sample, rng);
        if(result)
        {
            sample->imageX = Lerp(sample->imageX, xPixelStart, xPixelEnd);
            sample->imageY = Lerp(sample->imageY, yPixelStart, yPixelEnd);
        }

        return result;
    }

    int RoundSize(int sz) const { return sz; }

    Sampler *GetSubSampler(int num, int count)
    {
        if(count < m_numSamples)
            return NULL;

        // Share the number of samples to be generated to the `count` subsamplers.
        int x0, x1, y0, y1;
        ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
        int spc = m_numSamples / count;
        int r = m_numSamples % count;
        // FIXME: if m_numSamples is not divisible by count, the first subsampler gets the remaining,
        // this strategy is biased.
        int n = spc + (num == 0 ? r : 0);
        return new SparseSampler(x0, x1, y0, y1, n, shutterOpen, shutterClose);
    }

private:
    Sampler* m_sampler;
    int m_numSamples; // this is the plain number of samples, not in spp.
};


class MixSampler: public Sampler
{
public:
    MixSampler(Sampler* originalSampler, int xstart, int xend, int ystart, int yend, float sopen, float sclose, int n, bool isSPP):
        // NOTE: n is passed to Sampler here, but Sampler::samplesPerPixel is not used.
        Sampler(xstart, xend, ystart, yend, n, sopen, sclose)
    {
        m_spp = m_partSpp = m_rNumSamples = 0;
        m_mainSampler = m_fallbackSampler = m_sparseSampler = nullptr;

        if(isSPP)
        {
            SPP_CASE:
            int spp =  n;
            if(spp == originalSampler->RoundSize(spp))
            {
                // just use the main sampler
                m_spp = spp;
                m_partSpp = 0;
                m_stage = MAIN;
            }
            else
            {
                int roundSpp = spp - 1;
                while((roundSpp != originalSampler->RoundSize(roundSpp)) && (roundSpp > 0))
                    --roundSpp;
                if(roundSpp > 0)
                {
                    // use main sampler with roundSpp and fallback sampler with (spp - roundSpp)
                    m_spp = roundSpp;
                    m_partSpp = spp - roundSpp;
                    m_stage = MAIN;
                }
                else
                {
                    // use fallback sampler with spp
                    m_spp = 0;
                    m_partSpp = spp;
                    m_stage = FALLBACK;
                }
            }

            if(m_spp)
                m_mainSampler = createMainSampler();
            if(m_partSpp)
                m_fallbackSampler = new RandomSampler(xstart, xend, ystart, yend, m_partSpp, sopen, sclose);
        }
        else
        {
            int w = xend - xstart;
            int h = yend - ystart;
            int spp = n / (float)(w*h);
            m_rNumSamples = n % (w*h);
            if(m_rNumSamples)
            {
                m_sparseSampler = new SparseSampler(xstart, xend, ystart, yend, m_rNumSamples, sopen, sclose);
                m_stage = SPARSE;
            }

            if(spp > 0)
            {
                n = spp;
                goto SPP_CASE;
            }

        }
    }

    MixSampler(int xstart, int xend, int ystart, int yend, float sopen, float sclose, int spp, int partSpp, int rNumSamples):
        Sampler(xstart, xend, ystart, yend, spp, sopen, sclose)
    {
        m_mainSampler = m_fallbackSampler = m_sparseSampler = nullptr;
        m_spp = spp;
        m_partSpp = partSpp;
        m_rNumSamples = rNumSamples;

        if(m_rNumSamples)
        {
            m_sparseSampler = new SparseSampler(xstart, xend, ystart, yend, m_rNumSamples, sopen, sclose);
            m_stage = SPARSE;
        }
        if(m_partSpp)
        {
            m_fallbackSampler = new RandomSampler(xstart, xend, ystart, yend, m_partSpp, sopen, sclose);
            m_stage = FALLBACK;
        }
        if(m_spp)
        {
            m_mainSampler = createMainSampler();
            m_stage = MAIN;
        }
    }

    virtual ~MixSampler()
    {
        delete m_mainSampler;
        delete m_fallbackSampler;
        delete m_sparseSampler;
    }

    virtual int MaximumSampleCount()
    {
        int result = 0;
        if(m_mainSampler)
            result = std::max(result, m_mainSampler->MaximumSampleCount());
        if(m_fallbackSampler)
            result = std::max(result, m_fallbackSampler->MaximumSampleCount());
        if(m_sparseSampler)
            result = std::max(result, m_sparseSampler->MaximumSampleCount());
        return result;
    }

    virtual int GetMoreSamples(Sample *sample, RNG &rng)
    {
        int count = 0;

        switch(m_stage)
        {
        case MAIN:
            if((count = m_mainSampler->GetMoreSamples(sample, rng)) > 0)
                return count;
            else if(m_fallbackSampler)
                m_stage = FALLBACK;
            else if(m_sparseSampler)
            {
                m_stage = SPARSE;
                goto SPARSE_CASE;
            }
            else
                return 0;

        case FALLBACK:
            if((count = m_fallbackSampler->GetMoreSamples(sample, rng)) > 0)
                return count;
            else if(m_sparseSampler)
                m_stage = SPARSE;
            else
                return 0;

        case SPARSE:
            SPARSE_CASE:
            return m_sparseSampler->GetMoreSamples(sample, rng);
        }

        return count; //just to avoid warning
    }

    virtual Sampler *GetSubSampler(int num, int count)
    {
        int x0, x1, y0, y1;
        ComputeSubWindow(num, count, &x0, &x1, &y0, &y1);
        if (x0 == x1 || y0 == y1) return NULL;
        return new MixSampler(x0, x1, y0, y1, shutterOpen, shutterClose, m_spp, m_partSpp, m_rNumSamples);
    }

    virtual int RoundSize(int size) const
    { return size; }

private:
    enum Stage
    {
        MAIN,
        FALLBACK,
        SPARSE
    };

    Sampler* createMainSampler()
    {
        string name = g_samplerName;
        Sampler* sampler = nullptr;
        if (name == "bestcandidate")
            sampler = new BestCandidateSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, m_spp, shutterOpen, shutterClose);
        else if (name == "halton")
            sampler = new HaltonSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, m_spp, shutterOpen, shutterClose);
        else if (name == "lowdiscrepancy")
            sampler = new LDSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, m_spp, shutterOpen, shutterClose);
        else if (name == "random")
            sampler = new RandomSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, m_spp, shutterOpen, shutterClose);
        else if (name == "stratified")
        {
            int xsamp, ysamp = 0;
            if(isPerfectSquare(m_spp))
                xsamp = ysamp = sqrt(m_spp);
            else
            {
                xsamp = m_spp;
                ysamp = 1;
            }
            sampler = new StratifiedSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, xsamp, ysamp, true, shutterOpen, shutterClose);
        }

        return sampler;
    }

    Sampler* m_mainSampler;     // used when spp == m_mainSampler->RoundSize(spp)
    Sampler* m_fallbackSampler; // used when (spp != m_mainSampler->RoundSize(spp)) && (spp != m_fallbackSampler->RoundSize(spp))
    SparseSampler* m_sparseSampler;   // used for the remaining samples (m_sparseSampler < 1 spp)
    int m_spp;        // round spp count for the main sampler
    int m_partSpp;    // not round spp count. used in the fallback sampler
    int m_rNumSamples; // remaining samples for the sparse sampler
    Stage m_stage;
};


class AdaptiveMixSampler: public MixSampler
{
public:
    AdaptiveMixSampler(Sampler* originalSampler, int xstart, int xend, int ystart, int yend, float sopen, float sclose, int n, bool isSPP, const float* pdf):
        MixSampler(originalSampler, xstart, xend, ystart, yend, sopen, sclose, n, isSPP),
        dpdf(pdf, xend - xstart, yend - ystart)
    {}

    virtual int GetMoreSamples(Sample *samples, RNG &rng)
    {
        int count = MixSampler::GetMoreSamples(samples, rng);
        for(int i = 0; i < count; ++i)
        {
            int x = 0, y = 0;
            dpdf.sample(rng, &x, &y);
            samples[i].imageX = x + xPixelStart;
            samples[i].imageY = y + yPixelStart;
        }
        return count;
    }

private:
    DPDF2 dpdf;
};

//======================================================================


//======================================================================
//                           ServerRenderer
//======================================================================


// SamplerRenderer Method Definitions
ServerRenderer::ServerRenderer(Sampler *s, Camera *c,SurfaceIntegrator *si, VolumeIntegrator *vi, bool visIds) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    visualizeObjectIds = visIds;
}


ServerRenderer::~ServerRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


void ServerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();

    this->scene = scene;
    this->w = camera->film->xResolution;
    this->h = camera->film->yResolution;
    mainSample = new Sample(sampler, surfaceIntegrator, volumeIntegrator, scene);

    QEventLoop eventLoop;
    std::unique_ptr<RenderingServer> server(new RenderingServer);
    QObject::connect(server.get(), &RenderingServer::getSceneInfo,
        [this](SceneInfo* scene)
        { this->getSceneInfo(scene); }
    );
    QObject::connect(server.get(), &RenderingServer::evaluateSamples,
        [this](bool isSPP, int numSamples, int* resultSize)
        { this->evaluateSamples(isSPP, numSamples, resultSize); }
    );
    QObject::connect(server.get(), &RenderingServer::evaluateSamplesCrop,
        [this](bool isSPP, int numSamples, const CropWindow& crop, int* resultSize)
        { this->evaluateSamplesCrop(isSPP, numSamples, crop, resultSize); }
    );
    QObject::connect(server.get(), &RenderingServer::evaluateSamplesPDF,
        [this](bool isSPP, int numSamples, const float* pdf, int* resultSize)
        { this->evaluateSamplesPDF(isSPP, numSamples, pdf, resultSize); }
    );
    QObject::connect(server.get(), &RenderingServer::finishRender, &eventLoop, &QEventLoop::quit);
    server->startServer(2227);
    eventLoop.exec();

    PBRT_FINISHED_RENDERING();
}


Spectrum ServerRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T, SampleBuffer* sampleBuffer) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect))
        Li = surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena, sampleBuffer);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
}


Spectrum ServerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}

void ServerRenderer::getSceneInfo(SceneInfo *desc)
{
    desc->set<int>("width", camera->film->xResolution);
    desc->set<int>("height", camera->film->yResolution);
    desc->set<int>("max_spp", sampler->samplesPerPixel);
    desc->set<int>("max_samples", sampler->samplesPerPixel * camera->film->xResolution * camera->film->yResolution);

    // NOTE: In PBRT, shutterOpen and shutterClose defaults are 0 and 1, respectively,
    // which makes them not good to decide if a scene as motion blur.
    // This info (`has_motion_blur`) is being set in the benchmark configuration file.
    desc->set<float>("shutter_open", camera->shutterOpen);
    desc->set<float>("shutter_close", camera->shutterClose);

    bool hasAreaLights = false;
    for(Light* light : scene->lights)
    {
        if(light->IsDeltaLight() == false)
        {
            hasAreaLights = true;
            break;
        }
    }
    desc->set<bool>("has_area_lights", hasAreaLights);
}

void ServerRenderer::evaluateSamples(bool isSPP, int numSamples, int* resultSize)
{
    int totalNumSamples = isSPP ? w * h * numSamples : numSamples;
    *resultSize = totalNumSamples;

    MixSampler mixSampler(sampler,
                          0,
                          w,
                          0,
                          h,
                          sampler->shutterOpen,
                          sampler->shutterClose,
                          numSamples,
                          isSPP);

    SamplesPipe pipe;
    run(&mixSampler, mainSample, pipe);
}

void ServerRenderer::evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow &crop, int *resultSize)
{
    int totalNumSamples = isSPP ? w * h * numSamples : numSamples;
    *resultSize = totalNumSamples;

    MixSampler mixSampler(sampler,
                          crop.beginX,
                          crop.endX,
                          crop.beginY,
                          crop.endY,
                          sampler->shutterOpen,
                          sampler->shutterClose,
                          numSamples,
                          isSPP);

    SamplesPipe pipe;
    run(&mixSampler, mainSample, pipe);
}

void ServerRenderer::evaluateSamplesPDF(bool isSPP, int numSamples, const float *pdf, int *resultSize)
{
    int totalNumSamples = isSPP ? w * h * numSamples : numSamples;
    *resultSize = totalNumSamples;

    AdaptiveMixSampler mixSampler(sampler,
                                  0,
                                  w,
                                  0,
                                  h,
                                  sampler->shutterOpen,
                                  sampler->shutterClose,
                                  numSamples,
                                  isSPP,
                                  pdf);

    SamplesPipe pipe;
    run(&mixSampler, mainSample, pipe);
}

void ServerRenderer::run(Sampler* sampler, Sample* origSample, SamplesPipe& pipe)
{
    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(0);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            SampleBuffer sampleBuffer = pipe.getBuffer();

            // Replace sample positions if they are provided as INPUT
            samples[i].imageX = sampleBuffer.set(IMAGE_X, samples[i].imageX);
            samples[i].imageY = sampleBuffer.set(IMAGE_Y, samples[i].imageY);
            samples[i].lensU = sampleBuffer.set(LENS_U, samples[i].lensU);
            samples[i].lensV = sampleBuffer.set(LENS_V, samples[i].lensV);
            samples[i].time = sampleBuffer.set(TIME, samples[i].time);

            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (visualizeObjectIds) {
                if (rayWeight > 0.f && scene->Intersect(rays[i], &isects[i])) {
                    // random shading based on shape id...
                    uint32_t ids[2] = { isects[i].shapeId, isects[i].primitiveId };
                    uint32_t h = hash((char *)ids, sizeof(ids));
                    float rgb[3] = { float(h & 0xff), float((h >> 8) & 0xff),
                                     float((h >> 16) & 0xff) };
                    Ls[i] = Spectrum::FromRGB(rgb);
                    Ls[i] /= 255.f;
                }
                else
                    Ls[i] = 0.f;
            }
            else {
            if (rayWeight > 0.f)
                Ls[i] = rayWeight * Li(scene, rays[i], &samples[i], rng, arena, &isects[i], &Ts[i], &sampleBuffer);
            else {
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }

            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned "
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (std::isinf(Ls[i].y())) {
                Error("Infinite luminance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            }

            // Save features
            float rgb[3]; Ls[i].ToRGB(rgb);
            sampleBuffer.set(COLOR_R, rgb[0]);
            sampleBuffer.set(COLOR_G, rgb[1]);
            sampleBuffer.set(COLOR_B, rgb[2]);

            pipe << sampleBuffer;

            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _SamplerRendererTask_ is done with its image region
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
}



DPDF::DPDF(const float *fv, unsigned int n) :
    m_pdf(n),
    m_cpdf(n)
{
    float sum = std::accumulate(fv, fv + n, 0.f);
    m_funcInt = sum ;

    float nf = 1.f / sum;
    sum = 0.f;

    for(unsigned int i=0; i < n; ++i)
    {
        m_pdf[i] = fv[i] * nf;
        m_cpdf[i] = sum + m_pdf[i];
        sum += m_pdf[i];
    }

    //Garante a última entrada de m_cpdf será 1.f. Pode não ser devido a erro numérico.
    m_cpdf.back() = 1.f;
}

DPDF::DPDF(const vector<float>& fv) :
    m_pdf(fv.size()),
    m_cpdf(fv.size())
{
    float sum = accumulate(fv.begin(), fv.end(), 0.f);
    m_funcInt = sum ;

    float nf = 1.f / sum;
    sum = 0.f;

    for(unsigned int i=0; i < fv.size(); ++i)
    {
        m_pdf[i] = fv[i] * nf;
        m_cpdf[i] = sum + m_pdf[i];
        sum += m_pdf[i];
    }

    //Garante a última entrada de m_cpdf será 1.f. Pode não ser devido a erro numérico.
    m_cpdf.back() = 1.f;
}

unsigned int DPDF::sample(const RNG& rng, float* pdf)
{
    float *p = std::lower_bound(&(*m_cpdf.begin()), &(*m_cpdf.end()), rng.RandomFloat());
    unsigned int i = p - &(*m_cpdf.begin());

    if(pdf && m_pdf.size()) *pdf = m_pdf[i];

    return i;
}




DPDF2::DPDF2(const float* fv, int nx, int ny)
{
    m_dpdfX.reserve(nx);

    for(int x = 0; x < nx; ++x)
        m_dpdfX.push_back(new DPDF(&fv[x*ny], ny));

    vector<float> marginalFunc;
    marginalFunc.reserve(nx);

    for(int x = 0; x < nx; ++x)
        marginalFunc.push_back(m_dpdfX[x]->m_funcInt);

    m_marginal = new DPDF(marginalFunc);
}

DPDF2::~DPDF2()
{
    for(unsigned int i=0; i< m_dpdfX.size(); ++i)
        delete m_dpdfX[i];

    delete m_marginal;
}

void DPDF2::sample(const RNG& rng, int* x, int* y, float* pdf)
{
    float pdfs[2];

    *y = m_marginal->sample(rng, &pdfs[1]);
    *x = m_dpdfX[*y]->sample(rng, &pdfs[0]);

    if(pdf) *pdf = pdfs[0] * pdfs[1];
}
