#ifndef SERVERRENDERER_H
#define SERVERRENDERER_H

#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
#include <fbksd/core/SampleLayout.h>

class SceneInfo;
class CropWindow;

// SamplerRenderer Declarations
class ServerRenderer : public Renderer {
public:
    // SamplerRenderer Public Methods
    ServerRenderer(Sampler *s, Camera *c, SurfaceIntegrator *si,
                    VolumeIntegrator *vi, bool visIds);
    ~ServerRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL, SampleBuffer* sampleBuffer = nullptr) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;

    void getSceneInfo(SceneInfo *scene);
    bool evaluateSamples(int64_t spp, int64_t remainingCount);

private:
    std::vector<Task*> createTasks(Sampler* sppSampler,
                                   Sampler* sparseSampler,
                                   bool hasInputPos,
                                   int nSppTasks,
                                   int nSparseTasks,
                                   int numPixels,
                                   int spp,
                                   ProgressReporter& reporter);
    bool m_visualizeObjectIds;
    Sampler *m_sampler;
    Camera *m_camera;
    SurfaceIntegrator *m_surfaceIntegrator;
    VolumeIntegrator *m_volumeIntegrator;
    const Scene* m_scene;
    Sample* m_sample;
    SampleLayout m_layout;
    int m_width;
    int m_height;
};



// SamplerRendererTask Declarations
class ServerRendererTask : public Task {
public:
    // SamplerRendererTask Public Methods
    ServerRendererTask(const Scene *sc, Renderer *ren, Camera *c,
                        ProgressReporter &pr, Sampler *ms, Sample *sam,
                        bool visIds, int tn, int tc, size_t pipeOffset = 0, bool seekPipeByPixel = true)
      : reporter(pr)
    {
        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; visualizeObjectIds = visIds; taskNum = tn; taskCount = tc;
#ifdef PBRT2_RANDOM_SEEDING
        seed = rand();
#else
        seed = taskNum;
#endif
        m_pipeOffset = pipeOffset;
        m_seekPipeByPixel = seekPipeByPixel;
    }
    void Run();
private:
    const Scene *scene;
    const Renderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    bool visualizeObjectIds;
    int taskNum, taskCount;
    int seed;
    size_t m_pipeOffset;
    bool m_seekPipeByPixel;
};



#endif
