#include "renderers/serverrenderer.h"
#include "stdafx.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "samplers/stratified.h"
#include "samplers/random.h"
#include "samplers/lowdiscrepancy.h"
#include "Benchmark/RenderingServer/RenderingServer.h"
#include <QEventLoop>

namespace
{

uint32_t hash(char *key, uint32_t len)
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

bool isPerfectSquare(int n)
{
    if(n < 0)
        return false;
    int root(std::round(std::sqrt(n)));
    return n == root * root;
}


// Sampler used when the # of samples is not in SPP.
class SparseSampler: public Sampler
{
public:
    SparseSampler(int xstart, int xend, int ystart, int yend, int ns, float sopen, float sclose):
        Sampler(xstart, xend, ystart, yend, ns, sopen, sclose)
    {
        m_numSamples = ns;
        if(isPerfectSquare(m_numSamples))
        {
            int ss = std::sqrt(m_numSamples);
            m_sampler = new StratifiedSampler(0, 1, 0, 1, ss, ss, true, sopen, sclose);
        }
        else
        {
            int xs = m_numSamples;
            m_sampler = new StratifiedSampler(0, 1, 0, 1, xs, 1, true, sopen, sclose);
        }
    }

    ~SparseSampler()
    { delete m_sampler; }

    int MaximumSampleCount() { return m_sampler->MaximumSampleCount(); }

    int GetMoreSamples(Sample *sample, RNG &rng)
    {
        int result = m_sampler->GetMoreSamples(sample, rng);
        for(size_t i = 0; i < result; ++i)
        {
            sample[i].imageX = Lerp(sample[i].imageX, xPixelStart, xPixelEnd);
            sample[i].imageY = Lerp(sample[i].imageY, yPixelStart, yPixelEnd);
        }

        return result;
    }

    int RoundSize(int sz) const { return sz; }

    Sampler *GetSubSampler(int num, int count)
    {
        int sps = getSubSamplerSize(num, count);
        if(sps)
            return new SparseSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, sps, shutterOpen, shutterClose);
        else
            return nullptr;
    }

    int getSubSamplerSize(int num, int count) const
    {
        int sps = m_numSamples / count;
        int rest = m_numSamples % count;
        sps += num == 0 ? rest : 0;
        return sps;
    }

private:
    Sampler* m_sampler;
    int m_numSamples; // this is the plain number of samples, not in spp.
};

}


// ServerRendererTask Definitions
void ServerRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _ServerRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    SamplesPipe pipe;
    pipe.seek(m_pipeOffset);

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    // Get samples from _Sampler_ and update image
    int sampleCount;
    int s = 0;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        if(m_seekPipeByPixel && s == 0)
        {
            //WARNING: Assumes that sampleCount <= spp and that all sampleCount samples
            // are for the same pixel.
            int x = samples[0].imageX;
            int y = samples[0].imageY;
            pipe.seek(x, y, sampler->samplesPerPixel, camera->film->xResolution);
        }

        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            SampleBuffer sampleBuffer = pipe.getBuffer();
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
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i], &sampleBuffer);
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
            if(std::isfinite(rays[i].maxt))
                sampleBuffer.set(DEPTH, rays[i].maxt);
            else
                sampleBuffer.set(DEPTH, 0.f);
            pipe << sampleBuffer;

            if(++s == sampler->samplesPerPixel)
                s = 0;

            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }

    // Clean up after _ServerRendererTask_ is done with its image region
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// ServerRenderer Method Definitions
ServerRenderer::ServerRenderer(Sampler *s, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 bool visIds) {
    m_sampler = s;
    m_camera = c;
    m_surfaceIntegrator = si;
    m_volumeIntegrator = vi;
    m_visualizeObjectIds = visIds;
}


ServerRenderer::~ServerRenderer() {
    delete m_sampler;
    delete m_camera;
    delete m_surfaceIntegrator;
    delete m_volumeIntegrator;
}


void ServerRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    m_surfaceIntegrator->Preprocess(scene, m_camera, this);
    m_volumeIntegrator->Preprocess(scene, m_camera, this);
    m_scene = scene;
    m_width = m_camera->film->xResolution;
    m_height = m_camera->film->yResolution;
    m_sample = new Sample(m_sampler, m_surfaceIntegrator, m_volumeIntegrator, scene);
    PBRT_FINISHED_PREPROCESSING();

    QEventLoop eventLoop;
    std::unique_ptr<RenderingServer> server(new RenderingServer);
    QObject::connect(server.get(), &RenderingServer::setParameters,
        [this](int maxSPP, const SampleLayout& layout, float* samplesBuffer, float* pdfBuffer)
        { m_layout = layout; }
    );
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
        Li = m_surfaceIntegrator->Li(scene, this, ray, *isect, sample,
                                   rng, arena, sampleBuffer);
    else {
        // Handle ray that doesn't intersect any geometry
        for (uint32_t i = 0; i < scene->lights.size(); ++i)
           Li += scene->lights[i]->Le(ray);
    }
    Spectrum Lvi = m_volumeIntegrator->Li(scene, this, ray, m_sample, rng,
                                        T, arena);
    return *T * Li + Lvi;
}


Spectrum ServerRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return m_volumeIntegrator->Transmittance(scene, this, ray, m_sample,
                                           rng, arena);
}

void ServerRenderer::getSceneInfo(SceneInfo *info)
{
    info->set<int>("width", m_width);
    info->set<int>("height", m_height);
    info->set<int>("max_spp", m_sampler->samplesPerPixel);
    info->set<int>("max_samples", m_sampler->samplesPerPixel * m_width * m_height);

    // NOTE: In PBRT, shutterOpen and shutterClose defaults are 0 and 1, respectively,
    // which makes them not good to decide if a scene as motion blur.
    // This info (`has_motion_blur`) is being set in the benchmark configuration file.
    info->set<float>("shutter_open", m_camera->shutterOpen);
    info->set<float>("shutter_close", m_camera->shutterClose);

    bool hasAreaLights = false;
    for(Light* light : m_scene->lights)
    {
        if(light->IsDeltaLight() == false)
        {
            hasAreaLights = true;
            break;
        }
    }
    info->set<bool>("has_area_lights", hasAreaLights);
}

void ServerRenderer::evaluateSamples(bool isSPP, int numSamples, int *resultSize)
{
    PBRT_STARTED_RENDERING();
    int numPixels = m_width * m_height;
    int totalNumSamples = isSPP ? numPixels * numSamples : numSamples;
    *resultSize = totalNumSamples;
    int spp = isSPP ? numSamples : totalNumSamples / (float)(numPixels);
    int rest = totalNumSamples % numPixels;
    int nSppTasks = max(32 * NumSystemCores(), numPixels / (16*16));
    nSppTasks = RoundUpPow2(nSppTasks);
    int nRestTasks = max(1, rest / (16*16));
    nRestTasks = RoundUpPow2(nRestTasks);
    bool hasInputPos = m_layout.hasInput("IMAGE_X") || m_layout.hasInput("IMAGE_Y");

    std::unique_ptr<Sampler> sppSampler;
    std::unique_ptr<SparseSampler> sparseSampler;
    if(spp)
        sppSampler.reset(new RandomSampler(0, m_width, 0, m_height, spp,
                                           m_sampler->shutterOpen,
                                           m_sampler->shutterClose));
    else
        nSppTasks = 0;

    if(rest)
        sparseSampler.reset(new SparseSampler(0, m_width, 0, m_height, rest,
                                              m_sampler->shutterOpen,
                                              m_sampler->shutterClose));
    else
        nRestTasks = 0;

    ProgressReporter reporter(nSppTasks + nRestTasks, "Rendering");
    vector<Task*> renderTasks = createTasks(sppSampler.get(), sparseSampler.get(), hasInputPos, nSppTasks, nRestTasks, numPixels, spp, reporter);
    EnqueueTasks(renderTasks);
    WaitForAllTasks();
    for (size_t i = 0; i < renderTasks.size(); ++i)
        delete renderTasks[i];
    reporter.Done();
    PBRT_FINISHED_RENDERING();
}

void ServerRenderer::evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow &crop, int *resultSize)
{

}

void ServerRenderer::evaluateSamplesPDF(bool isSPP, int numSamples, const float *pdf, int *resultSize)
{

}

std::vector<Task*> ServerRenderer::createTasks(Sampler* sppSampler,
                                               Sampler* sparseSampler,
                                               bool hasInputPos,
                                               int nSppTasks,
                                               int nSparseTasks,
                                               int numPixels,
                                               int spp,
                                               ProgressReporter& reporter)
{
    std::vector<Task*> renderTasks;
    size_t sppOffset = 0;
    if(sppSampler)
    {
        if(hasInputPos)
        {
            size_t currentOffset = 0;
            for (int i = 0; i < nSppTasks; ++i)
            {
                renderTasks.push_back(new ServerRendererTask(m_scene, this, m_camera,
                                                             reporter, sppSampler, m_sample,
                                                             m_visualizeObjectIds,
                                                             nSppTasks-1-i, nSppTasks, currentOffset, false));
                currentOffset += size_t(sppSampler->getSubSamplerSize(nSppTasks-1-i, nSppTasks)) * m_layout.getSampleSize() * spp;
            }
        }
        else
        {
            for (int i = 0; i < nSppTasks; ++i)
                renderTasks.push_back(new ServerRendererTask(m_scene, this, m_camera,
                                                             reporter, sppSampler, m_sample,
                                                             m_visualizeObjectIds,
                                                             nSppTasks-1-i, nSppTasks, 0, true));
        }
        sppOffset = size_t(numPixels) * spp * m_layout.getSampleSize();
    }

    if(sparseSampler)
    {
        size_t currentOffset = sppOffset;
        for (int i = 0; i < nSparseTasks; ++i)
        {
            renderTasks.push_back(new ServerRendererTask(m_scene, this, m_camera,
                                                         reporter, sparseSampler, m_sample,
                                                         m_visualizeObjectIds,
                                                         nSparseTasks-1-i, nSparseTasks, currentOffset, false));
            currentOffset += sparseSampler->getSubSamplerSize(nSparseTasks-1-i, nSparseTasks) * m_layout.getSampleSize();
        }
    }

    return renderTasks;
}
