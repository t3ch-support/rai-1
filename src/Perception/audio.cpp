#include "audio.h"

#include <Core/util.h>

#ifdef HAVE_LIBAV
extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
}
#include "avutil.h"
extern Mutex libav_open_mutex; // is defined in videoEncoder.cpp right now

#define DEFAULT_CONTAINER "wav"

class sAudioWriter_libav {
private:
    AVCodec* codec;
    AVFormatContext *oc;
    AVStream *s;

public:
    sAudioWriter_libav(const char* filename, unsigned int sample_rate=44100, unsigned int num_channels=2) : codec(NULL), oc(NULL), s(NULL) {
        Lock lock(libav_open_mutex);
        avcodec_register_all();
        oc = avformat_alloc_context();
        if(!oc) {
            HALT("Could not allocate format context");
        }
        oc->oformat = mt_guess_format(filename, DEFAULT_CONTAINER);
        oc->audio_codec_id = CODEC_ID_PCM_S16LE;       
        codec = avcodec_find_encoder(oc->audio_codec_id);
        if(!codec)
            HALT("Audion codec not found");

        snprintf(oc->filename, sizeof(oc->filename), "%s", filename);
        s= avformat_new_stream(oc, codec);
        if(!s) {
            HALT("Could not allocate stream structure");
        }
        s->codec->sample_fmt = AV_SAMPLE_FMT_S16;
        s->codec->sample_rate = sample_rate;
        s->codec->channels = num_channels;

        // some formats want stream headers to be separate
        if(oc->oformat->flags & AVFMT_GLOBALHEADER)
            oc->flags |= CODEC_FLAG_GLOBAL_HEADER;

        /* open it */
        if (avcodec_open2(s->codec, codec, NULL) < 0)
            HALT("Encoder failed to open");

        if (avio_open(&(oc->pb), filename, URL_WRONLY) < 0) {
            HALT("Could not open " << filename);
        }
        avformat_write_header(oc, NULL);
    }
    ~sAudioWriter_libav() {
        avcodec_close(s->codec);
        avformat_free_context(oc);
    }

    void write(const byteA& audio_samples) {
        AVPacket pkt;
        av_init_packet(&pkt);
        pkt.data = NULL;

        AVFrame frame;
        frame.data[0] = audio_samples.p;
        frame.linesize[0] = audio_samples.d0;
        frame.nb_samples = audio_samples.d0 / 4; // 2 channels, 2 bytes per sample
        frame.pts = AV_NOPTS_VALUE;

        // note: allocates data, so not most efficient, but should do
        int got_packet, ret;

        if((ret = avcodec_encode_audio2(s->codec, &pkt, &frame, &got_packet)) != 0) {
            HALT("Could not encode audio frame: " << ret);
            return;
        }
        if(got_packet) {
            if((ret = av_write_frame(oc, &pkt)) != 0) {
                HALT("Error while writing audio frame: " << ret);
            }
        }
        av_free_packet(&pkt);
    }
};
#else
class sAudioRecorder_libav {};
#endif

AudioWriter_libav::AudioWriter_libav(const char* filename)
#ifdef HAVE_LIBAV
    : s(new sAudioWriter_libav(filename))
#endif
{
}
AudioWriter_libav::~AudioWriter_libav() {
#if HAVE_LIBAV
    delete(s);
#endif
}

void AudioWriter_libav::writeSamples_R44100_2C_S16_NE(const byteA &samples) {
#ifdef HAVE_LIBAV
    s->write(samples);
#endif
}

// PULSEAUDIO-based audio grabbing implementation

#ifdef HAVE_PULSEAUDIO
#include <pulse/simple.h>
#include <pulse/error.h>

class sAudioPoller_PA {
private:
    pa_simple *pa;
public:
    sAudioPoller_PA(const char* appname, const char* dev) {
        static const pa_sample_spec ss = {
            .format = PA_SAMPLE_S16NE,
            .rate = 44100,
            .channels = 2
        };
        int error;
        if (!(pa = pa_simple_new(NULL, appname, PA_STREAM_RECORD, dev, "record", &ss, NULL, NULL, &error))) {
            HALT(": pa_simple_new() failed: " << pa_strerror(error));
        }
    }
    ~sAudioPoller_PA() {
        pa_simple_free(pa);
    }
    int read(byteA& buf) {
        int ret, error;
        if((ret = pa_simple_read(pa, buf.p, buf.d0, &error)) < 0) {
            MT_MSG(pa_strerror(error));
            return -1;
        }
        return ret;
    }
};

#endif

AudioPoller_PA::AudioPoller_PA(const char* appname, const char* dev)
#ifdef HAVE_PULSEAUDIO
    : s(new sAudioPoller_PA(appname, dev))
#endif
{
}

AudioPoller_PA::~AudioPoller_PA() {
#ifdef HAVE_PULSEAUDIO
    delete s;
#endif
}

bool AudioPoller_PA::read(byteA& buf) {
#ifdef HAVE_PULSEAUDIO
    return(s->read(buf) >= 0);
#else
    MT_MSG("AudioPoller_PA::read not available, libpulse missing");
#endif
}

