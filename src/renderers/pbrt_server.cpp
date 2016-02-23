
#include "pbrt_server.h"

#include "Benchmark/RenderingServer/RenderingServer.h"
#include "serverrenderer.h"


PBRTServer::PBRTServer(ServerRenderer* renderer):
    QObject(0)
{
    this->renderer = renderer;

    server = new RenderingServer;
    QObject::connect(server, &RenderingServer::getSceneInfo, this, &PBRTServer::getSceneInfo);
    QObject::connect(server, &RenderingServer::setMaxSPP, this, &PBRTServer::setMaxSPP);
    QObject::connect(server, &RenderingServer::setSampleLayout, this, &PBRTServer::setSampleLayout);
    QObject::connect(server, &RenderingServer::evaluateSamples, this, &PBRTServer::evaluateSamples);
    QObject::connect(server, &RenderingServer::evaluateSamplesCrop, this, &PBRTServer::evaluateSamplesCrop);
    QObject::connect(server, &RenderingServer::evaluateSamplesPDF, this, &PBRTServer::evaluateSamplesPDF);
    QObject::connect(server, &RenderingServer::finishRender, this, &PBRTServer::finishRender);
    server->startServer(2227);
}

PBRTServer::~PBRTServer()
{
    delete server;
}

void PBRTServer::setSampleBuffers(float *input, float *output)
{
    server->setSampleBuffers(input, output);
}

void PBRTServer::getSceneInfo(SceneInfo *scene)
{
    renderer->getSceneInfo(scene);
}

void PBRTServer::setMaxSPP(int maxSPP)
{
    renderer->setMaxSPP(maxSPP);
}

void PBRTServer::setSampleLayout(const SampleLayout& layout)
{
    renderer->setSampleLayout(layout);
}

void PBRTServer::evaluateSamples(bool isSPP, int numSamples, int* resultSize)
{
    renderer->evaluateSamples(isSPP, numSamples, resultSize);
}

void PBRTServer::evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow &crop, int *resultSize)
{
    renderer->evaluateSamples(isSPP, numSamples, crop, resultSize);
}

void PBRTServer::evaluateSamplesPDF(bool isSPP, int numSamples, const float *pdf, int *resultSize)
{
    renderer->evaluateSamples(isSPP, numSamples, pdf, resultSize);
}
