#include <iostream>
using namespace std;

#include "samplers/adaptive.h"
#include "samplers/bestcandidate.h"
#include "samplers/halton.h"
#include "samplers/lowdiscrepancy.h"
#include "samplers/random.h"
#include "samplers/stratified.h"

extern string g_samplerName;


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

    ~MixSampler()
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
            if(IsPowerOf2(m_spp))
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


int main()
{
    g_samplerName = "lowdiscrepancy";
    Sampler* origSampler = new LDSampler(0, 5, 0, 5, 3, 0.f, 1.f);
    Sampler* mixSampler = new MixSampler(origSampler, 0, 5, 0, 5, 0.f, 1.f, 1, true);
    Sample origSample(mixSampler, NULL, NULL, NULL);
    int maxSamples = mixSampler->MaximumSampleCount();
    Sample *samples = origSample.Duplicate(maxSamples);
    RNG rng;

//    int sampleCount = 0;
//    size_t totalSamples = 0;
//    while ((sampleCount = mixSampler->GetMoreSamples(samples, rng)) > 0)
//    {
//        totalSamples += sampleCount;
//        if(totalSamples == 100)
//            totalSamples = 100;

//        for(int s = 0; s < sampleCount; ++s)
//            cout << samples[s].imageX << " " << samples[s].imageY << endl;
//    }

    int nTasks = 2;
    std::vector<Sampler*> samplers;
    for (int i = 0; i < nTasks; ++i)
        samplers.push_back(mixSampler->GetSubSampler(i, nTasks));

    for (int i = 0; i < nTasks; ++i)
    {
        int sampleCount = 0;
        size_t totalSamples = 0;
        while ((sampleCount = samplers[i]->GetMoreSamples(samples, rng)) > 0)
        {
            totalSamples += sampleCount;
            if(totalSamples == 100)
                totalSamples = 100;

            for(int s = 0; s < sampleCount; ++s)
                cout << samples[s].imageX << " " << samples[s].imageY << endl;
        }
    }

    return 0;
}
