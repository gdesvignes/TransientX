/*
 * rfi.cpp
 *
 *  Created on: May 6, 2020
 *      Author: ypmen
 */

#include "string.h"

#include "rfi.h"
#include "kdtree.h"
#include "dedisperse.h"

using namespace std;

RFI::RFI(){}

RFI::RFI(const RFI &rfi) : DataBuffer<float>(rfi)
{
    weights = rfi.weights;
}

RFI & RFI::operator=(const RFI &rfi)
{
    DataBuffer<float>::operator=(rfi);

    weights = rfi.weights;

    return *this;  
}

RFI::~RFI(){}

void RFI::prepare(const DataBuffer<float> &databuffer)
{
    nsamples = databuffer.nsamples;
    nchans = databuffer.nchans;

    resize(nsamples, nchans);

    tsamp = databuffer.tsamp;
    frequencies = databuffer.frequencies;

    weights.resize(nchans, 1);
}

void RFI::zap(DataBuffer<float> &databuffer, const vector<pair<double, double>> &zaplist)
{
    for (long int j=0; j<nchans; j++)
    {
        for (auto k=zaplist.begin(); k!=zaplist.end(); ++k)
        {
            if (frequencies[j]>=(*k).first and frequencies[j]<=(*k).second)
            {
                weights[j] = 0;
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples; i++)
    {   
        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]*weights[j];
        }
    }

    equalized = databuffer.equalized;
    counter += nsamples;
}

void RFI::zdot(DataBuffer<float> &databuffer)
{
    vector<float> xe(nchans, 0.);
    vector<float> xs(nchans, 0.);
    vector<float> alpha(nchans, 0.);
    vector<float> beta(nchans, 0.);
    vector<float> s(nsamples, 0.);
    float se = 0.;
    float ss = 0.;

    for (long int i=0; i<nsamples; i++)
    {
        float temp = 0.;
        for (long int j=0; j<nchans; j++)
        {
            temp += databuffer.buffer[i*nchans+j];
        }

        for (long int j=0; j<nchans; j++)
        {
            xe[j] += databuffer.buffer[i*nchans+j];
        }

        temp /= nchans;
        se += temp;
        ss += temp*temp;
        
        for (long int j=0; j<nchans; j++)
        {
            xs[j] += databuffer.buffer[i*nchans+j]*temp;
        }

        s[i] = temp;
    }

    float tmp = se*se-ss*nsamples;
    for (long int j=0; j<nchans; j++)
    {
        alpha[j] = (xe[j]*se-xs[j]*nsamples)/tmp;
        beta[j] = (xs[j]*se-xe[j]*ss)/tmp;
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-alpha[j]*s[i]-beta[j];
        }
    }

    equalized = false;
}

void RFI::zero(DataBuffer<float> &databuffer)
{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples; i++)
    {
        float s = 0;
        for (long int j=0; j<nchans; j++)
        {
            s += databuffer.buffer[i*nchans+j];
        }
        s /= nchans;

        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = databuffer.buffer[i*nchans+j]-s;
        }
    }

    equalized = false;
}

bool RFI::mask(DataBuffer<float> &databuffer, float threRFI2, int td, int fd)
{
    long int nsamples_ds = nsamples/td;
    long int nchans_ds = nchans/fd;

    vector<float> buffer_ds(nsamples_ds*nchans_ds, 0.);

    //downsample
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples_ds; i++)
    {
        for (long int n=0; n<td; n++)
        {
            for (long int k=0; k<fd; k++)
            {
                for (long int j=0; j<nchans_ds; j++)
                {
                    buffer_ds[i*nchans_ds+j] += databuffer.buffer[(i*td+n)*nchans+j*fd+k];
                }
            }
        }
    }

    buffer = databuffer.buffer;

    vector<float> buffer_dscopy = buffer_ds;
	std::nth_element(buffer_dscopy.begin(), buffer_dscopy.begin()+nsamples_ds*nchans_ds/4, buffer_dscopy.end(), std::less<float>());
    float Q1 = buffer_dscopy[nsamples_ds*nchans_ds/4];
	std::nth_element(buffer_dscopy.begin(), buffer_dscopy.begin()+nsamples_ds*nchans_ds/2, buffer_dscopy.end(), std::less<float>());
    float Q2 = buffer_dscopy[nsamples_ds*nchans_ds/2];
    std::nth_element(buffer_dscopy.begin(), buffer_dscopy.begin()+nsamples_ds*nchans_ds/4, buffer_dscopy.end(), std::greater<float>());
    float Q3 = buffer_dscopy[nsamples_ds*nchans_ds/4];

	float mean = Q2;
	float var = ((Q3-Q1)/1.349)*((Q3-Q1)/1.349);
    float thre = threRFI2*var;

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples_ds; i++)
    {
        for (long int n=0; n<td; n++)
        {
            for (long int j=0; j<nchans_ds; j++)
            {
                for (long int k=0; k<fd; k++)
                {
                    if ((buffer_ds[i*nchans_ds+j]-mean)*(buffer_ds[i*nchans_ds+j]-mean)>thre)
                        buffer[(i*td+n)*nchans+j*fd+k] = mean;
                }
            }
        }
    }

    equalized = databuffer.equalized;

    return true;
}

bool RFI::kadaneF(DataBuffer<float> &databuffer, float threRFI2, double widthlimit, int td, int fd)
{
    if (!databuffer.equalized)
    {
        cerr<<"Error: data is not equalize"<<endl;
        return false;
    }

    vector<float> bufferT(nchans*nsamples, 0.);

    transpose_pad<float>(&bufferT[0], &databuffer.buffer[0], nsamples, nchans);

    long int nsamples_ds = nsamples/td;
    long int nchans_ds = nchans/fd;

    vector<float> bufferT_ds(nchans_ds*nsamples_ds, 0.);

    //downsample
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int j=0; j<nchans_ds; j++)
    {
        for (long int k=0; k<fd; k++)
        {
            for (long int n=0; n<td; n++)
            {
                for (long int i=0; i<nsamples_ds; i++)       
                {
                    bufferT_ds[j*nsamples_ds+i] += bufferT[(j*fd+k)*nsamples+(i*td+n)];
                }
            }
        }
    }

#ifdef _OPENMP
    float *chdata_t = new float [num_threads*nsamples_ds];
    memset(chdata_t, 0, sizeof(float)*num_threads*nsamples_ds);
#else
    float *chdata_t = new float [nsamples_ds];
    memset(chdata_t, 0, sizeof(float)*nsamples_ds);
#endif

    int wnlimit = widthlimit/tsamp/td;

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int j=0; j<nchans_ds; j++)
    {
#ifdef _OPENMP
        float *chdata = chdata_t+omp_get_thread_num()*nsamples_ds;
#else
        float *chdata = chdata_t;
#endif

        for (long int i=0; i<nsamples_ds; i++)
        {
            chdata[i] = bufferT_ds[j*nsamples_ds+i];
        }

        long int start, end;
        float snr2=0.;
        float boxsum = kadane<float>(chdata, nsamples_ds, &start, &end);
        long int wn = end-start+1;

        float mean = 0.;
        float var = td*fd;    
        if (wn > wnlimit)
        {
            snr2 = boxsum*boxsum/(wn*var);
        }
        else
        {
            snr2 = 0.;
        }

        start *= td;
        end += 1;
        end *= td;
        if (snr2 > threRFI2)
        {
            for (long int k=0; k<fd; k++)
            {
                for (long int i=start; i<end; i++)
                {
                    bufferT[(j*fd+k)*nsamples+i] = 0.;
                }
            }
        }

        //<0
        for (long int i=0; i<nsamples_ds; i++)
        {
            chdata[i] = -chdata[i];
        }

        snr2=0.;
        boxsum = kadane<float>(chdata, nsamples_ds, &start, &end);
        wn = end-start+1;

        if (wn > wnlimit)
        {
            snr2 = boxsum*boxsum/(wn*var);
        }
        else
        {
            snr2 = 0.;
        }

        start *= td;
        end += 1;
        end *= td;
        if (snr2 > threRFI2)
        {
            for (long int k=0; k<fd; k++)
            {
                for (long int i=start; i<end; i++)
                {
                    bufferT[(j*fd+k)*nsamples+i] = 0.;
                }
            }
        }
        
    }

    transpose_pad<float>(&buffer[0], &bufferT[0], nchans, nsamples);

    equalized = databuffer.equalized;

    delete [] chdata_t;

    return true;
}

bool RFI::kadaneT(DataBuffer<float> &databuffer, float threRFI2, double bandlimit, int td, int fd)
{
    if (!databuffer.equalized)
    {
        cerr<<"Error: data is not equalize"<<endl;
        return false;
    }

    long int nsamples_ds = nsamples/td;
    long int nchans_ds = nchans/fd;

    vector<float> buffer_ds(nsamples_ds*nchans_ds, 0.);

    //downsample
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples_ds; i++)
    {
        for (long int n=0; n<td; n++)
        {
            for (long int k=0; k<fd; k++)
            {
                for (long int j=0; j<nchans_ds; j++)
                {
                    buffer_ds[i*nchans_ds+j] += databuffer.buffer[(i*td+n)*nchans+j*fd+k];
                }
            }
        }
    }

    buffer = databuffer.buffer;

#ifdef _OPENMP
    float *tsdata_t = new float [num_threads*nchans_ds];
    memset(tsdata_t, 0, sizeof(float)*num_threads*nchans_ds);
#else
    float *tsdata_t = new float [nchans_ds];
    memset(tsdata_t, 0, sizeof(float)*nchans_ds);
#endif

    int chnlimit = abs(bandlimit/(frequencies[1]-frequencies[0])/fd);

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
    for (long int i=0; i<nsamples_ds; i++)
    {
#ifdef _OPENMP
        float *tsdata = tsdata_t+omp_get_thread_num()*nchans_ds;
#else
        float *tsdata = tsdata_t;
#endif

        for (long int j=0; j<nchans_ds; j++)
        {
            tsdata[j] = buffer_ds[i*nchans_ds+j];
        }

        long int start, end;
        float snr2=0.;
        float boxsum = kadane<float>(tsdata, nchans_ds, &start, &end);
        long int chn = end-start+1;

        float mean = 0.;
        float var = td*fd;    

        if (chn > chnlimit)
        {
            snr2 = boxsum*boxsum/(chn*var);
        }
        else
        {
            snr2 = 0.;
        }

        start *= fd;
        end += 1;
        end *= fd;
        if (chn > nchans_ds*0.8)
        {
            start = 0;
            end = nchans-1;
        }

        if (snr2 > threRFI2)
        {
            for (long int k=0; k<td; k++)
            {
                for (long int j=start; j<end; j++)
                {
                    buffer[(i*td+k)*nchans+j] = 0.;
                }
            }
        }

        //<0
        for (long int j=0; j<nchans_ds; j++)
        {
            tsdata[j] = -tsdata[j];
        }

        snr2=0.;
        boxsum = kadane<float>(tsdata, nchans_ds, &start, &end);
        chn = end-start+1;

        if (chn > chnlimit)
        {
            snr2 = boxsum*boxsum/(chn*var);
        }
        else
        {
            snr2 = 0.;
        }

        start *= fd;
        end += 1;
        end *= fd;
        
        if (chn > nchans_ds*0.8)
        {
            start = 0;
            end = nchans-1;
        }

        if (snr2 > threRFI2)
        {
            for (long int k=0; k<td; k++)
            {
                for (long int j=start; j<end; j++)
                {
                    buffer[(i*td+k)*nchans+j] = 0.;
                }
            }
        }
    }

    equalized = databuffer.equalized;

    delete [] tsdata_t;

    return true;
}
