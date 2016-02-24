#ifndef PBRT_SERVER_H
#define PBRT_SERVER_H


#include <QObject>

class SceneInfo;
class RenderingServer;
class ServerRenderer;
class SampleLayout;
class CropWindow;


class PBRTServer : public QObject
{
    Q_OBJECT
public:
    PBRTServer(ServerRenderer *renderer);
    ~PBRTServer();

signals:
    void finishRender();

private slots:
    void getSceneInfo(SceneInfo* scene);
    void evaluateSamples(bool isSPP, int numSamples, int* resultSize);
    void evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow& crop, int* resultSize);
    void evaluateSamplesPDF(bool isSPP, int numSamples, const float* pdf, int* resultSize);

private:
    RenderingServer* server;
    ServerRenderer* renderer;
};



#endif // PBRT_SERVER_H
